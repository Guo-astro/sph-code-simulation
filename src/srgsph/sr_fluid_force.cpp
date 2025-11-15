#include "defines.hpp"
#include "particle.hpp"
#include "periodic.hpp"
#include "simulation.hpp"
#include "bhtree.hpp"
#include "kernel/kernel_function.hpp"
#include "srgsph/sr_fluid_force.hpp"
#include "srgsph/sr_primitive_recovery.hpp"

#ifdef EXHAUSTIVE_SEARCH
#include "exhaustive_search.hpp"
#endif

namespace sph
{
namespace srgsph
{

void FluidForce::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::FluidForce::initialize(param);
    m_is_2nd_order = param->srgsph.is_2nd_order;
    m_gamma = param->physics.gamma;
    m_c_speed = param->srgsph.c_speed;
    m_c_shock = param->srgsph.c_shock;
    m_c_cd = param->srgsph.c_cd;

    exact_riemann_solver();
}

// van Leer (1979) limiter - same as non-relativistic GSPH
inline real limiter(const real dq1, const real dq2)
{
    const real dq1dq2 = dq1 * dq2;
    if(dq1dq2 <= 0) {
        return 0.0;
    } else {
        return 2.0 * dq1dq2 / (dq1 + dq2);
    }
}

/**
 * Main SR-GSPH calculation
 * Computes dS/dt (canonical momentum) and de/dt (canonical energy)
 * Based on Eqs. 64-65 in paper
 */
void FluidForce::calculation(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();
    auto * tree = sim->get_tree().get();
    const real dt = sim->get_dt();

    // For MUSCL reconstruction
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
        
        // Neighbor search
#ifdef EXHAUSTIVE_SEARCH
        int const n_neighbor = exhaustive_search(p_i, p_i.sml, particles, num,
                                                neighbor_list, m_neighbor_number * neighbor_list_size,
                                                periodic, true);
#else
        int const n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, true);
#endif

        // Initialize forces
        const vec_t & r_i = p_i.pos;
        const vec_t & v_i = p_i.vel;
        const real h_i = p_i.sml;
        const real nu_i = p_i.nu;

        vec_t dS(0.0);    // dS/dt (canonical momentum)
        real de = 0.0;    // de/dt (canonical energy)

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
            
            // Project velocities onto line connecting particles
            const real ve_i = inner_product(v_i, e_ij);
            const real ve_j = inner_product(p_j.vel, e_ij);
            
            real vstar, pstar;

            if(m_is_2nd_order) {
                // MUSCL reconstruction for 2nd order accuracy
                // Murante et al. (2011) adapted for SR

                real right[4], left[4];  // [velocity, density, pressure, sound]
                const real delta_i = 0.5 * (1.0 - p_i.sound * dt * r_inv);
                const real delta_j = 0.5 * (1.0 - p_j.sound * dt * r_inv);

                // Velocity reconstruction
                const real dv_ij = ve_i - ve_j;
                vec_t dv_i, dv_j;
                for(int k = 0; k < DIM; ++k) {
                    dv_i[k] = inner_product(grad_v[k][i], e_ij);
                    dv_j[k] = inner_product(grad_v[k][j], e_ij);
                }
                const real dve_i = inner_product(dv_i, e_ij) * r;
                const real dve_j = inner_product(dv_j, e_ij) * r;
                
                // Monotonicity constraint (Eq. 66)
                // Switch to 1st order if:
                // 1. Strong shock: C_shock * e_ij·(v_i - v_j) > min(c_s,i, c_s,j)
                // 2. Large pressure jump: |log10(P_i/P_j)| > C_c.d.
                const real min_cs = std::min(p_i.sound, p_j.sound);
                const real shock_indicator = m_c_shock * std::abs(ve_i - ve_j);
                const real pressure_jump = std::abs(std::log10(p_i.pres / p_j.pres));
                const bool is_shock = (shock_indicator > min_cs) || (pressure_jump > m_c_cd);

                if(is_shock) {
                    // Fall back to 1st order at shocks
                    right[0] = ve_i;
                    left[0] = ve_j;
                    right[1] = p_i.dens;
                    left[1] = p_j.dens;
                    right[2] = p_i.pres;
                    left[2] = p_j.pres;
                    right[3] = p_i.sound;
                    left[3] = p_j.sound;
                } else {
                    // 2nd order MUSCL
                    right[0] = ve_i - limiter(dv_ij, dve_i) * delta_i;
                    left[0] = ve_j + limiter(dv_ij, dve_j) * delta_j;

                    // Density
                    const real dd_ij = p_i.dens - p_j.dens;
                    const real dd_i = inner_product(grad_d[i], e_ij) * r;
                    const real dd_j = inner_product(grad_d[j], e_ij) * r;
                    right[1] = p_i.dens - limiter(dd_ij, dd_i) * delta_i;
                    left[1] = p_j.dens + limiter(dd_ij, dd_j) * delta_j;

                    // Pressure
                    const real dp_ij = p_i.pres - p_j.pres;
                    const real dp_i = inner_product(grad_p[i], e_ij) * r;
                    const real dp_j = inner_product(grad_p[j], e_ij) * r;
                    right[2] = p_i.pres - limiter(dp_ij, dp_i) * delta_i;
                    left[2] = p_j.pres + limiter(dp_ij, dp_j) * delta_j;

                    // Sound speed (recompute from reconstructed state)
                    right[3] = std::sqrt(m_gamma * right[2] / right[1]);
                    left[3] = std::sqrt(m_gamma * left[2] / left[1]);
                }

                m_solver(left, right, pstar, vstar);
            } else {
                // 1st order Riemann solver
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

            // Kernel gradients (variable smoothing length)
            // Eq. 64-65: uses ∇_i W(x_i-x_j, 2h_i) - ∇_j W(x_i-x_j, 2h_j)
            const vec_t dw_i = kernel->dw(r_ij, r, 2.0 * h_i);
            const vec_t dw_j = kernel->dw(r_ij, r, 2.0 * p_j.sml);
            const vec_t dw_ij = dw_i - dw_j;

            // Volume weighting V²_ij (simplified - could use interpolation)
            // For now, use average of particle volumes
            const real V_i = nu_i / p_i.N;
            const real V_j = p_j.nu / p_j.N;
            const real V2_ij = 0.5 * (V_i * V_i + V_j * V_j);

            // Equation of motion (Eq. 64):
            // ⟨ν_i Ṡ_i⟩ = -Σ_j P*_ij V²_ij [∇_i W - ∇_j W]
            dS -= dw_ij * (pstar * V2_ij);

            // Energy equation (Eq. 65):
            // ⟨ν_i ė_i⟩ = -Σ_j P*_ij v*_ij · V²_ij [∇_i W - ∇_j W]
            const vec_t v_ij = e_ij * vstar;
            de -= inner_product(v_ij, dw_ij) * (pstar * V2_ij);
        }

        // Store time derivatives
        // Note: These are ν_i * dS/dt and ν_i * de/dt
        // Divide by ν_i to get actual derivatives
        p_i.dS = dS / nu_i;
        p_i.de = de / nu_i;
    }
}

/**
 * HLL Riemann solver for special relativity
 * Simplified version - full exact solver would use Pons et al. (2000)
 * 
 * For now, use relativistic HLL which is robust and simple
 */
void FluidForce::exact_riemann_solver()
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

        // Relativistic HLL (simplified)
        // Wave speeds in SR are more complex, but for moderate velocities
        // this approximation works reasonably well

        // Roe-averaged state
        const real roe_l = std::sqrt(rho_l);
        const real roe_r = std::sqrt(rho_r);
        const real roe_inv = 1.0 / (roe_l + roe_r);

        const real u_t = (roe_l * u_l + roe_r * u_r) * roe_inv;
        const real c_t = (roe_l * c_l + roe_r * c_r) * roe_inv;

        // Wave speeds (need relativistic correction for high velocity)
        const real c_light = m_c_speed;
        const real gamma_l = 1.0 / std::sqrt(1.0 - u_l * u_l / (c_light * c_light));
        const real gamma_r = 1.0 / std::sqrt(1.0 - u_r * u_r / (c_light * c_light));

        // Approximate SR wave speeds
        const real s_l = std::min(u_l - c_l * gamma_l, u_t - c_t);
        const real s_r = std::max(u_r + c_r * gamma_r, u_t + c_t);

        // HLL fluxes
        const real c1 = rho_l * (s_l - u_l);
        const real c2 = rho_r * (s_r - u_r);
        const real c3 = 1.0 / (c1 - c2);
        const real c4 = p_l - u_l * c1;
        const real c5 = p_r - u_r * c2;
        
        vstar = (c5 - c4) * c3;
        pstar = (c1 * c5 - c2 * c4) * c3;
    };
}

}
}
