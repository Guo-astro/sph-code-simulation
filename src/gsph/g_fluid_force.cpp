#include "defines.hpp"
#include "particle.hpp"
#include "periodic.hpp"
#include "simulation.hpp"
#include "bhtree.hpp"
#include "kernel/kernel_function.hpp"
#include "gsph/g_fluid_force.hpp"

#ifdef EXHAUSTIVE_SEARCH
#include "exhaustive_search.hpp"
#endif

namespace sph
{
namespace gsph
{

void FluidForce::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::FluidForce::initialize(param);
    m_is_2nd_order = param->gsph.is_2nd_order;
    m_gamma = param->physics.gamma;

    // Select Riemann solver based on configuration
    if(param->gsph.riemann_solver == RiemannSolverType::ITERATIVE) {
        iterative_solver();
    } else if(param->gsph.riemann_solver == RiemannSolverType::KITAJIMA) {
        kitajima_solver();
    } else {
        hll_solver();
    }
}

// van Leer (1979) limiter
inline real limiter(const real dq1, const real dq2)
{
    const real dq1dq2 = dq1 * dq2;
    if(dq1dq2 <= 0) {
        return 0.0;
    } else {
        return 2.0 * dq1dq2 / (dq1 + dq2);
    }
}

// Cha & Whitworth (2003)
void FluidForce::calculation(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();
    auto * tree = sim->get_tree().get();
    const real dt = sim->get_dt();

    // for MUSCL
    auto & grad_d = sim->get_vector_array("grad_density");
    auto & grad_p = sim->get_vector_array("grad_pressure");
    vec_t * grad_v[DIM] = {
        sim->get_vector_array("grad_velocity_0").data(),
#if DIM == 2
        sim->get_vector_array("grad_velocity_1").data(),
#elif DIM == 3
        sim->get_vector_array("grad_velocity_1").data(),
        sim->get_vector_array("grad_velocity_2").data(),
#endif
    };

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);
        
        // neighbor search
#ifdef EXHAUSTIVE_SEARCH
        int const n_neighbor = exhaustive_search(p_i, p_i.sml, particles, num, neighbor_list, m_neighbor_number * neighbor_list_size, periodic, true);
#else
        int const n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, true);
#endif

        // fluid force
        const vec_t & r_i = p_i.pos;
        const vec_t & v_i = p_i.vel;
        const real h_i = p_i.sml;
        const real rho2_inv_i = 1.0 / sqr(p_i.dens);

        vec_t acc(0.0);
        real dene = 0.0;

        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
            const real r = std::abs(r_ij);

            if(r >= std::max(h_i, p_j.sml) || r == 0.0) {
                continue;
            }

            const real r_inv = 1.0 / r;
            const vec_t e_ij = r_ij * r_inv;
            const real ve_i = inner_product(v_i, e_ij);
            const real ve_j = inner_product(p_j.vel, e_ij);
            real vstar, pstar;

            if(m_is_2nd_order) {
                // Murante et al. (2011)

                real right[4], left[4];
                const real delta_i = 0.5 * (1.0 - p_i.sound * dt * r_inv);
                const real delta_j = 0.5 * (1.0 - p_j.sound * dt * r_inv);

                // velocity
                const real dv_ij = ve_i - ve_j;
                vec_t dv_i, dv_j;
                for(int k = 0; k < DIM; ++k) {
                    dv_i[k] = inner_product(grad_v[k][i], e_ij);
                    dv_j[k] = inner_product(grad_v[k][j], e_ij);
                }
                const real dve_i = inner_product(dv_i, e_ij) * r;
                const real dve_j = inner_product(dv_j, e_ij) * r;
                right[0] = ve_i - limiter(dv_ij, dve_i) * delta_i;
                left[0] = ve_j + limiter(dv_ij, dve_j) * delta_j;

                // density
                const real dd_ij = p_i.dens - p_j.dens;
                const real dd_i = inner_product(grad_d[i], e_ij) * r;
                const real dd_j = inner_product(grad_d[j], e_ij) * r;
                right[1] = p_i.dens - limiter(dd_ij, dd_i) * delta_i;
                left[1] = p_j.dens + limiter(dd_ij, dd_j) * delta_j;

                // pressure
                const real dp_ij = p_i.pres - p_j.pres;
                const real dp_i = inner_product(grad_p[i], e_ij) * r;
                const real dp_j = inner_product(grad_p[j], e_ij) * r;
                right[2] = p_i.pres - limiter(dp_ij, dp_i) * delta_i;
                left[2] = p_j.pres + limiter(dp_ij, dp_j) * delta_j;

                // sound speed
                right[3] = std::sqrt(m_gamma * right[2] / right[1]);
                left[3] = std::sqrt(m_gamma * left[2] / left[1]);

                m_solver(left, right, pstar, vstar);
            } else {
                const real right[4] = {
                    ve_i,
                    p_i.dens,
                    p_i.pres,
                    p_i.sound,
                };
                const real left[4] = {
                    ve_j,
                    p_j.dens,
                    p_j.pres,
                    p_j.sound,
                };

                m_solver(left, right, pstar, vstar);
            }

            const vec_t dw_i = kernel->dw(r_ij, r, h_i);
            const vec_t dw_j = kernel->dw(r_ij, r, p_j.sml);
            const vec_t v_ij = e_ij * vstar;
            const real rho2_inv_j = 1.0 / sqr(p_j.dens);
            const vec_t f = dw_i * (p_j.mass * pstar * rho2_inv_i) + dw_j * (p_j.mass * pstar * rho2_inv_j);

            acc -= f;
            dene -= inner_product(f, v_ij - v_i);
        }

        p_i.acc = acc;
        p_i.dene = dene;
    }
}

void FluidForce::hll_solver()
{
    m_solver = [&](const real left[], const real right[], real & pstar, real & vstar) {
        const real u_l   = left[0];
        const real rho_l = left[1];
        const real p_l   = left[2];
        const real c_l   = left[3];

        const real u_r   = right[0];
        const real rho_r = right[1];
        const real p_r   = right[2];
        const real c_r   = right[3];

        const real roe_l = std::sqrt(rho_l);
        const real roe_r = std::sqrt(rho_r);
        const real roe_inv = 1.0 / (roe_l + roe_r);

        const real u_t = (roe_l * u_l + roe_r * u_r) * roe_inv;
        const real c_t = (roe_l * c_l + roe_r * c_r) * roe_inv;
        const real s_l = std::min(u_l - c_l, u_t - c_t);
        const real s_r = std::max(u_r + c_r, u_t + c_t);

        const real c1 = rho_l * (s_l - u_l);
        const real c2 = rho_r * (s_r - u_r);
        const real c3 = 1.0 / (c1 - c2);
        const real c4 = p_l - u_l * c1;
        const real c5 = p_r - u_r * c2;
        
        vstar = (c5 - c4) * c3;
        pstar = (c1 * c5 - c2 * c4) * c3;
    };
}

// Iterative Riemann solver from van Leer (1997) 
// Implementation based on PySPH's van_leer function
// Reference: https://github.com/pypr/pysph/blob/main/pysph/sph/gas_dynamics/riemann_solver.py
void FluidForce::iterative_solver()
{
    m_solver = [&](const real left[], const real right[], real & pstar, real & vstar) {
        // Left and right states (a = left, b = right)
        const real u_l   = left[0];   // velocity
        const real rho_l = left[1];   // density
        const real p_l   = left[2];   // pressure
        // left[3] is sound speed from state, but we compute Lagrangian sound speed below
        
        const real u_r   = right[0];
        const real rho_r = right[1];
        const real p_r   = right[2];
        // right[3] is sound speed from state, but we compute Lagrangian sound speed below

        // Check for negative values
        if (rho_l < 0.0 || rho_r < 0.0 || p_l < 0.0 || p_r < 0.0) {
            pstar = 0.0;
            vstar = 0.0;
            return;
        }

        // Constants
        const real gamma2 = 1.0 + m_gamma;
        const real gamma1 = 0.5 * gamma2 / m_gamma;
        const real smallp = 1e-25;
        
        // Specific volumes
        const real V_l = 1.0 / rho_l;
        const real V_r = 1.0 / rho_r;
        
        // Lagrangian sound speeds
        const real cl = std::sqrt(m_gamma * p_l * rho_l);
        const real cr = std::sqrt(m_gamma * p_r * rho_r);
        
        // Initial guess for pstar
        real pstar_guess = p_r - p_l - cr * (u_r - u_l);
        pstar_guess = p_l + pstar_guess * cl / (cl + cr);
        pstar_guess = std::max(pstar_guess, smallp);
        
        // Newton-Raphson iteration
        const int max_iter = 20;
        const real tol = 1e-6;
        real pstar_iter = pstar_guess;
        bool converged = false;
        
        for (int iter = 0; iter < max_iter; ++iter) {
            const real pstar_old = pstar_iter;
            
            // Left wave impedance
            real w_l = 1.0 + gamma1 * (pstar_iter - p_l) / p_l;
            w_l = cl * std::sqrt(w_l);
            
            // Right wave impedance  
            real w_r = 1.0 + gamma1 * (pstar_iter - p_r) / p_r;
            w_r = cr * std::sqrt(w_r);
            
            // Left derivative
            real z_l = 4.0 * V_l * w_l * w_l;
            z_l = -z_l * w_l / (z_l - gamma2 * (pstar_iter - p_l));
            
            // Right derivative
            real z_r = 4.0 * V_r * w_r * w_r;
            z_r = z_r * w_r / (z_r - gamma2 * (pstar_iter - p_r));
            
            // Intermediate velocities
            const real ustar_l = u_l - (pstar_iter - p_l) / w_l;
            const real ustar_r = u_r + (pstar_iter - p_r) / w_r;
            
            // Newton-Raphson update
            pstar_iter = pstar_iter + (ustar_r - ustar_l) * (z_l * z_r) / (z_r - z_l);
            pstar_iter = std::max(smallp, pstar_iter);
            
            // Check convergence
            converged = (std::abs(pstar_iter - pstar_old) / pstar_iter < tol);
            if (converged) {
                break;
            }
        }
        
        // Recalculate wave impedances with final pstar
        real w_l = 1.0 + gamma1 * (pstar_iter - p_l) / p_l;
        w_l = cl * std::sqrt(w_l);
        
        real w_r = 1.0 + gamma1 * (pstar_iter - p_r) / p_r;
        w_r = cr * std::sqrt(w_r);
        
        // Calculate averaged ustar
        const real ustar_l = u_l - (pstar_iter - p_l) / w_l;
        const real ustar_r = u_r + (pstar_iter - p_r) / w_r;
        const real ustar = 0.5 * (ustar_l + ustar_r);
        
        // Return results
        pstar = pstar_iter;
        vstar = ustar;
    };
}

// Kitajima-style iterative Riemann solver
// Non-relativistic adaptation inspired by Kitajima et al. (2025) methodology
// Uses Newton-Raphson iteration with explicit shock/rarefaction wave treatment
void FluidForce::kitajima_solver()
{
    m_solver = [&](const real left[], const real right[], real & pstar, real & vstar) {
        const real u_l   = left[0];   // velocity
        const real rho_l = left[1];   // density
        const real p_l   = left[2];   // pressure
        const real c_l   = left[3];   // sound speed
        
        const real u_r   = right[0];
        const real rho_r = right[1];
        const real p_r   = right[2];
        const real c_r   = right[3];

        // Check for negative values
        if (rho_l < 0.0 || rho_r < 0.0 || p_l < 0.0 || p_r < 0.0) {
            pstar = 0.0;
            vstar = 0.0;
            return;
        }

        const real smallp = 1e-25;
        const real tol = 1e-6;
        const int max_iter = 50;
        
        // Compute initial pressure guess using two-shock approximation
        // Similar to PVRS (Primitive Variable Riemann Solver)
        const real rho_avg = 0.5 * (rho_l + rho_r);
        const real c_avg = 0.5 * (c_l + c_r);
        real pstar_guess = 0.5 * (p_l + p_r) - 0.5 * (u_r - u_l) * rho_avg * c_avg;
        pstar_guess = std::max(pstar_guess, smallp);
        
        // Helper lambda for computing post-wave state
        auto get_velocity_jump = [&](const real p, const real p_k, const real rho_k, 
                                      const real c_k, const real u_k, bool is_left) -> std::pair<real, real> {
            const real sign = is_left ? -1.0 : 1.0;
            real du_dp;  // derivative of velocity jump w.r.t. pressure
            real du;     // velocity jump
            
            if (p > p_k) {
                // Shock wave
                // Rankine-Hugoniot relations
                const real A = 2.0 / ((m_gamma + 1.0) * rho_k);
                const real B = (m_gamma - 1.0) / (m_gamma + 1.0) * p_k;
                const real sqrt_term = std::sqrt(A / (p + B));
                
                du = (p - p_k) * sqrt_term;
                du_dp = sqrt_term * (1.0 - 0.5 * (p - p_k) / (p + B));
            } else {
                // Rarefaction wave
                // Isentropic relations
                const real gamma_exp = (m_gamma - 1.0) / (2.0 * m_gamma);
                const real p_ratio = p / p_k;
                const real factor = std::pow(p_ratio, gamma_exp) - 1.0;
                
                du = 2.0 * c_k / (m_gamma - 1.0) * factor;
                du_dp = (1.0 / (m_gamma * p)) * std::pow(p_ratio, gamma_exp) * 2.0 * c_k / (m_gamma - 1.0);
            }
            
            return {sign * du, sign * du_dp};
        };
        
        // Newton-Raphson iteration to find pstar
        real pstar_iter = pstar_guess;
        bool converged = false;
        
        for (int iter = 0; iter < max_iter; ++iter) {
            // Left wave contribution
            auto [du_l, du_dp_l] = get_velocity_jump(pstar_iter, p_l, rho_l, c_l, u_l, true);
            
            // Right wave contribution  
            auto [du_r, du_dp_r] = get_velocity_jump(pstar_iter, p_r, rho_r, c_r, u_r, false);
            
            // Function f(p) = u_l + du_l - (u_r + du_r) = 0
            const real f = u_l + du_l - u_r - du_r;
            const real df_dp = du_dp_l - du_dp_r;
            
            // Newton-Raphson update
            const real pstar_new = pstar_iter - f / df_dp;
            pstar_iter = std::max(pstar_new, smallp);
            
            // Check convergence
            const real rel_change = std::abs(pstar_new - pstar_iter) / (pstar_iter + smallp);
            if (rel_change < tol) {
                converged = true;
                break;
            }
        }
        
        // Compute star velocity
        auto [du_l_final, _] = get_velocity_jump(pstar_iter, p_l, rho_l, c_l, u_l, true);
        auto [du_r_final, __] = get_velocity_jump(pstar_iter, p_r, rho_r, c_r, u_r, false);
        
        const real vstar_l = u_l + du_l_final;
        const real vstar_r = u_r + du_r_final;
        
        // Return results
        pstar = pstar_iter;
        vstar = 0.5 * (vstar_l + vstar_r);  // Average the two estimates
    };
}

}
}