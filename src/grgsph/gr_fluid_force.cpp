/**
 * GR-GSPH Fluid Force Implementation
 *
 * Godunov SPH for General Relativity using exact Riemann solver.
 * Based on:
 * - Kitajima et al. (2025) SR-GSPH
 * - Liptai & Price (2019) GRSPH
 */

#include "grgsph/gr_fluid_force.hpp"
#include "grgsph/gr_riemann.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "periodic.hpp"
#include "openmp.hpp"
#include "kernel/kernel_function.hpp"
#include "bhtree.hpp"
#include "logger.hpp"
#include "defines.hpp"
#include "exception.hpp"
#include <cmath>
#include <algorithm>

#ifdef EXHAUSTIVE_SEARCH
#include "exhaustive_search.hpp"
#endif

namespace sph {
namespace grgsph {

void GRFluidForce::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::FluidForce::initialize(param);

    m_gamma = param->physics.gamma;
    m_c_speed = param->srgsph.c_speed;
    m_c_shock = (param->srgsph.c_shock > 0.0) ? param->srgsph.c_shock : 3.0;
    m_c_cd = (param->srgsph.c_cd > 0.0) ? param->srgsph.c_cd : 0.2;
    m_use_muscl = param->srgsph.is_2nd_order;
    m_riemann_solver = param->srgsph.riemann_solver;

    // Default to Minkowski if no metric specified
    if (!m_metric) {
        m_metric = std::make_unique<MinkowskiMetric>();
    }
}

/**
 * Solve interface state using Riemann solver
 */
void GRFluidForce::solve_interface_state(
    const real left_state[5],
    const real right_state[5],
    real& P_star,
    real& v_x_star,
    real& v_t_star) const
{
    // Build GRRiemannState from input arrays
    // Array format: [v^x, n, P, c_s, v^t]
    GRRiemannState left, right;

    left.V_x = left_state[0];
    left.n = left_state[1];
    left.P = left_state[2];
    left.c_s = left_state[3];
    left.V_t = left_state[4];
    left.h = 1.0 + m_gamma / (m_gamma - 1.0) * left.P / std::max(left.n, 1e-15);

    const real v2_L = left.V_x * left.V_x + left.V_t * left.V_t;
    left.W = 1.0 / std::sqrt(std::max(1e-10, 1.0 - v2_L));

    right.V_x = right_state[0];
    right.n = right_state[1];
    right.P = right_state[2];
    right.c_s = right_state[3];
    right.V_t = right_state[4];
    right.h = 1.0 + m_gamma / (m_gamma - 1.0) * right.P / std::max(right.n, 1e-15);

    const real v2_R = right.V_x * right.V_x + right.V_t * right.V_t;
    right.W = 1.0 / std::sqrt(std::max(1e-10, 1.0 - v2_R));

    // Select solver type (use exact for EXACT, KITAJIMA, ITERATIVE; use HLLC for HLLC, HLL)
    const bool use_exact = (m_riemann_solver == RiemannSolverType::EXACT ||
                            m_riemann_solver == RiemannSolverType::KITAJIMA ||
                            m_riemann_solver == RiemannSolverType::ITERATIVE);

    // Solve Riemann problem
    GRRiemannResult result = GRRiemannSolver::solve(left, right, m_gamma, use_exact);

    P_star = result.P_star;
    v_x_star = result.v_x_star;
    v_t_star = result.v_t_star;
}

/**
 * Compute gravitational source term
 *
 * f_i = (√-g / 2ρ*) T^μν ∂g_μν/∂x^i
 *
 * The stress-energy tensor for ideal fluid:
 * T^μν = ρ h u^μ u^ν + P g^μν
 *
 * In 3+1 form:
 * T^00 = ρ h Γ² / α² - P / α²
 * T^0i = ρ h Γ² V^i / α
 * T^ij = ρ h Γ² V^i V^j + P γ^ij
 */
vec_t GRFluidForce::compute_metric_source(
    const GRPrimitive& prim,
    const Metric31& metric,
    const MetricDerivatives& derivs,
    real rho_star) const
{
    vec_t f_source(0.0);

    if (rho_star < 1e-15) return f_source;

    const real n = prim.n;
    const real P = prim.P;
    const real h = prim.h;
    const real Gamma = prim.Gamma;
    const real Gamma2 = Gamma * Gamma;
    const real alpha = metric.alpha;
    const real alpha2 = alpha * alpha;

    // Stress-energy tensor components
    const real rho_h = n * h;

    // T^00
    const real T00 = rho_h * Gamma2 / alpha2 - P / alpha2;

    // T^0i
    real T0i[3];
    for (int i = 0; i < 3; ++i) {
        T0i[i] = rho_h * Gamma2 * prim.V[i] / alpha;
    }

    // T^ij
    real Tij[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Tij[i][j] = rho_h * Gamma2 * prim.V[i] * prim.V[j] + P * metric.gamma_inv[i][j];
        }
    }

    // Compute source term f_k = (√-g / 2ρ*) T^μν ∂g_μν/∂x^k
    const real prefactor = metric.sqrt_neg_g / (2.0 * rho_star);

    for (int k = 0; k < DIM; ++k) {
        // Contribution from g_tt derivative
        real contribution = T00 * derivs.dg_tt[k];

        // Contribution from g_ti derivatives (factor of 2 for symmetry)
        for (int i = 0; i < 3; ++i) {
            contribution += 2.0 * T0i[i] * derivs.dg_ti[i][k];
        }

        // Contribution from γ_ij derivatives
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                contribution += Tij[i][j] * derivs.dgamma_ij[i][j][k];
            }
        }

        f_source[k] = prefactor * contribution;
    }

    return f_source;
}

/**
 * Normalize derivatives (same as SR-GSPH)
 * Scale derivatives by baryon number so dS and de are per-baryon quantities.
 */
void GRFluidForce::normalize_derivatives(SPHParticle& p)
{
    if (p.nu <= 0.0) {
        THROW_ERROR("GR-GSPH particle baryon number must be positive.");
    }

    const real inv_nu = 1.0 / p.nu;
    p.dS *= inv_nu;
    p.de *= inv_nu;
    p.dS_t *= inv_nu;  // Normalize tangent momentum derivative for 1D tests
    p.acc = p.dS;      // Keep legacy aliases in sync
    p.dene = p.de;
}

/**
 * Main force calculation loop
 *
 * For each particle pair:
 * 1. Compute local metrics at both positions
 * 2. Transform velocities to Eulerian frame
 * 3. Solve Riemann problem for interface state
 * 4. Compute pressure gradient force
 * 5. Add gravitational source term
 */
void GRFluidForce::calculation(std::shared_ptr<Simulation> sim)
{
    auto& particles = sim->get_particles();
    auto* periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto* kernel = sim->get_kernel().get();
    auto* tree = sim->get_tree().get();

    // Reset derivatives
#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        particles[i].dS = vec_t(0.0);
        particles[i].de = 0.0;
    }

    const real sqrt_two = std::sqrt(2.0);

    // Main force loop
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < num; ++i) {
        auto& p_i = particles[i];

        if (p_i.is_ghost) continue;

        // Compute metric at particle i position
        Metric31 metric_i;
        MetricDerivatives derivs_i;
        if (m_metric) {
            m_metric->compute_all(p_i.pos, metric_i, derivs_i);
        }

        // Find neighbors
        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

#ifdef EXHAUSTIVE_SEARCH
        const int n_neighbor = exhaustive_search(p_i, p_i.sml * 6.0, particles, num,
                                                  neighbor_list, m_neighbor_number * neighbor_list_size,
                                                  periodic, false);
#else
        const int n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, false);
#endif

        // Local accumulator for forces
        vec_t dS_acc(0.0);
        real de_acc = 0.0;

        for (int n = 0; n < n_neighbor; ++n) {
            const int j = neighbor_list[n];
            if (j == i) continue;

            auto& p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(p_i.pos, p_j.pos);
            const real r = std::abs(r_ij);

            if (r > 3.0 * sqrt_two * std::max(p_i.sml, p_j.sml)) continue;

            const vec_t n_ij = r_ij * (-1.0 / r);

            // Compute metric at particle j position (if different from flat)
            Metric31 metric_j;
            if (m_metric) {
                m_metric->compute(p_j.pos, metric_j);
            }

            // For now, use simple primitive variables from particles
            // In full implementation, would use GRParticle structure
            real p_L = p_i.pres;
            real rho_L = p_i.dens;
            real vx_L = p_i.vel[0];
#if DIM >= 2
            real vy_L = p_i.vel[1];
#endif
#if DIM == 3
            real vz_L = p_i.vel[2];
#endif
            real v_t_L = p_i.vel_t;

            real p_R = p_j.pres;
            real rho_R = p_j.dens;
            real vx_R = p_j.vel[0];
#if DIM >= 2
            real vy_R = p_j.vel[1];
#endif
#if DIM == 3
            real vz_R = p_j.vel[2];
#endif
            real v_t_R = p_j.vel_t;

            // Build velocity vectors
#if DIM == 1
            const vec_t v_L(vx_L);
            const vec_t v_R(vx_R);
#elif DIM == 2
            const vec_t v_L(vx_L, vy_L);
            const vec_t v_R(vx_R, vy_R);
#else
            const vec_t v_L(vx_L, vy_L, vz_L);
            const vec_t v_R(vx_R, vy_R, vz_R);
#endif

            // Project velocities onto line of sight
            real v_x_L = inner_product(v_L, n_ij);
            real v_x_R = inner_product(v_R, n_ij);

            // Riemann solution
            real P_star = 0.0;
            real v_x_star = 0.0;
            real v_t_star = 0.0;

            if (p_j.is_ghost) {
                // Wall boundary condition
                P_star = p_L;
                v_x_star = 0.0;
                v_t_star = v_t_L;
            } else {
                // State arrays: [v^x, n, P, c_s, v^t]
                real left_state[5] = {v_x_L, rho_L, p_L, p_i.sound, v_t_L};
                real right_state[5] = {v_x_R, rho_R, p_R, p_j.sound, v_t_R};

                solve_interface_state(left_state, right_state, P_star, v_x_star, v_t_star);
            }

            // Reconstruct star velocity vector
            vec_t v_star_vec = n_ij * v_x_star;
#if DIM >= 2
            // Add tangent component (upwind from contact wave direction)
            const vec_t v_t_vec_L = v_L - n_ij * v_x_L;
            const vec_t v_t_vec_R = v_R - n_ij * v_x_R;
            const real v_t_L_mag = std::abs(v_t_vec_L);
            const real v_t_R_mag = std::abs(v_t_vec_R);

            if (v_t_star > 1e-10) {
                if (v_x_star > 0.0 && v_t_L_mag > 1e-10) {
                    v_star_vec = v_star_vec + v_t_vec_L * (v_t_star / v_t_L_mag);
                } else if (v_x_star <= 0.0 && v_t_R_mag > 1e-10) {
                    v_star_vec = v_star_vec + v_t_vec_R * (v_t_star / v_t_R_mag);
                }
            }
#endif

            // Volume elements
            const real V_i = p_i.nu / p_i.N;
            const real V_j = p_j.nu / p_j.N;

            // Kernel gradients
            const vec_t grad_W_i = kernel->dw(r_ij, r, sqrt_two * p_i.sml);
            const vec_t grad_W_j = kernel->dw(r_ij, r, sqrt_two * p_j.sml);

            // Grad-h corrections
            const real omega_i = p_i.gradh;
            const real omega_j = p_j.gradh;

            // Force: F_i = -P* (V_i² Ω_i ∇W_i + V_j² Ω_j ∇W_j)
            const vec_t force = grad_W_i * (-P_star * V_i * V_i * omega_i)
                              + grad_W_j * (-P_star * V_j * V_j * omega_j);

            // Power: dE/dt = v* · F
            const real power = inner_product(v_star_vec, force);

            // Accumulate
            for (int d = 0; d < DIM; ++d) {
                dS_acc[d] += force[d];
            }
            de_acc += power;
        }

        // Add gravitational source term
        if (m_metric) {
            // Build primitive state for source computation
            GRPrimitive prim;
            prim.n = p_i.dens;
            prim.P = p_i.pres;
            prim.h = 1.0 + m_gamma / (m_gamma - 1.0) * prim.P / std::max(prim.n, 1e-15);

            const real v2 = inner_product(p_i.vel, p_i.vel) + p_i.vel_t * p_i.vel_t;
            prim.Gamma = 1.0 / std::sqrt(std::max(1e-10, 1.0 - v2));

            for (int d = 0; d < DIM; ++d) {
                prim.V[d] = p_i.vel[d];
            }

            // Compute ρ* = n Γ (for Minkowski √γ = 1)
            const real rho_star = prim.n * prim.Gamma * metric_i.sqrt_gamma;

            // Gravitational source
            vec_t f_grav = compute_metric_source(prim, metric_i, derivs_i, rho_star);

            for (int d = 0; d < DIM; ++d) {
                dS_acc[d] += f_grav[d];
            }
        }

        // Store results
        p_i.dS = dS_acc;
        p_i.de = de_acc;
    }

    // Normalize derivatives
#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        if (particles[i].is_ghost) continue;
        normalize_derivatives(particles[i]);
    }
}

} // namespace grgsph
} // namespace sph
