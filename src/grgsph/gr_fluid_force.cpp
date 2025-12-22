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

    // Geodesic mode: skip pairwise SPH forces (for testing geodesic motion)
    m_geodesic_mode = param->grgsph.geodesic_mode;
    m_bh_mass = param->grgsph.bh_mass;

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
 * Compute gravitational momentum source term (GRGSPH Formulation §5.3)
 *
 * From the formulation, the gravitational force is:
 *   f_i = (√-g / 2N*) T^μν ∂g_μν/∂x^i
 *
 * Implemented using the equivalent Valencia formulation (L&P 2019 Eq. 15):
 *   S_i^{source} = α√γ Γ^α_{iβ} T_α^β
 *
 * where:
 *   Γ^α_{iβ} = g^{αν} Γ_{νiβ}  (Christoffel symbol with raised first index)
 *   T_α^β = g_{αμ} T^{μβ}      (mixed stress-energy tensor)
 *
 * For Schwarzschild (zero shift, static):
 *   - Γ^0_{i0} = (1/α) ∂α/∂x^i    (gravitational redshift)
 *   - Γ^j_{ik} = spatial Christoffel symbols (velocity-dependent terms)
 *
 * The source is converted to per-baryon acceleration for primitive recovery.
 */
vec_t GRFluidForce::compute_metric_source(
    const GRPrimitive& prim,
    const Metric31& metric,
    const MetricDerivatives& derivs,
    real N_star) const
{
    vec_t f_source(0.0);

    if (N_star < 1e-15) return f_source;

    const real n = prim.n;
    const real P = prim.P;
    const real h = prim.h;
    const real W = prim.Gamma;  // Lorentz factor
    const real W2 = W * W;
    const real alpha = metric.alpha;
    const real alpha2 = alpha * alpha;

    // =========================================================================
    // STEP 1: Compute contravariant stress-energy tensor T^μν
    // =========================================================================
    const real rho_h = n * h;

    // T^0i = ρh W² V^i / α
    real T_up_0i[3];
    for (int i = 0; i < 3; ++i) {
        T_up_0i[i] = rho_h * W2 * prim.V[i] / alpha;
    }

    // T^ij = ρh W² V^i V^j + P γ^ij
    real T_up_ij[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            T_up_ij[i][j] = rho_h * W2 * prim.V[i] * prim.V[j] + P * metric.gamma_inv[i][j];
        }
    }

    // =========================================================================
    // STEP 2: Compute mixed stress-energy tensor T_α^β = g_{αμ} T^{μβ}
    // =========================================================================
    // For zero shift: g_00 = -α², g_0i = 0, g_ij = γ_ij

    // T_0^0 = g_{0μ} T^{μ0} = g_00 T^00 = -α² × (ρhW²/α² - P/α²) = -ρhW² + P
    const real T_mix_0_0 = -rho_h * W2 + P;

    // T_0^k = g_{0μ} T^{μk} = g_00 T^0k = -α² × (ρhW² V^k/α) = -α ρhW² V^k
    real T_mix_0_k[3];
    for (int k = 0; k < 3; ++k) {
        T_mix_0_k[k] = -alpha * rho_h * W2 * prim.V[k];
    }

    // T_j^0 = g_{jμ} T^{μ0} = g_jk T^k0 = γ_jk × (ρhW² V^k/α)
    real T_mix_j_0[3];
    for (int j = 0; j < 3; ++j) {
        T_mix_j_0[j] = 0.0;
        for (int k = 0; k < 3; ++k) {
            T_mix_j_0[j] += metric.gamma_ij[j][k] * rho_h * W2 * prim.V[k] / alpha;
        }
    }

    // T_j^k = g_{jμ} T^{μk} = g_jl T^lk = γ_jl (ρhW² V^l V^k + P γ^lk)
    //       = ρhW² V_j V^k + P δ_j^k
    real T_mix_j_k[3][3];
    for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 3; ++k) {
            // V_j = γ_jl V^l (lowered velocity)
            real V_low_j = 0.0;
            for (int l = 0; l < 3; ++l) {
                V_low_j += metric.gamma_ij[j][l] * prim.V[l];
            }
            T_mix_j_k[j][k] = rho_h * W2 * V_low_j * prim.V[k];
            if (j == k) T_mix_j_k[j][k] += P;
        }
    }

    // =========================================================================
    // STEP 3: Compute Christoffel symbols Γ^α_{iβ}
    // =========================================================================
    // For static metric with zero shift:
    //   Γ^0_{i0} = (1/α) ∂α/∂x^i
    //   Γ^0_{ik} = 0
    //   Γ^j_{i0} = 0
    //   Γ^j_{ik} = (1/2) γ^{jl} (∂γ_{li}/∂x^k + ∂γ_{lk}/∂x^i - ∂γ_{ik}/∂x^l)

    // Γ^0_{i0}: computed from lapse gradient
    // ∂α/∂x^i can be obtained from dg_tt = ∂(-α²)/∂x^i = -2α ∂α/∂x^i
    real Gamma_0_i0[3];
    for (int i = 0; i < 3; ++i) {
        // dg_tt[i] = -2α ∂α/∂x^i, so ∂α/∂x^i = -dg_tt[i]/(2α)
        real dalpha_dxi = -derivs.dg_tt[i] / (2.0 * alpha);
        Gamma_0_i0[i] = dalpha_dxi / alpha;
    }

    // Γ^j_{ik}: spatial Christoffel symbols
    // Γ_{lik} = (1/2)(∂γ_{li}/∂x^k + ∂γ_{lk}/∂x^i - ∂γ_{ik}/∂x^l)
    // Γ^j_{ik} = γ^{jl} Γ_{lik}
    real Gamma_j_ik[3][3][3];  // Gamma_j_ik[j][i][k]
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            for (int k = 0; k < 3; ++k) {
                real Gamma_lik = 0.0;
                for (int l = 0; l < 3; ++l) {
                    // Γ_{lik} = (1/2)(∂γ_{li}/∂x^k + ∂γ_{lk}/∂x^i - ∂γ_{ik}/∂x^l)
                    real term = 0.5 * (derivs.dgamma_ij[l][i][k] +
                                       derivs.dgamma_ij[l][k][i] -
                                       derivs.dgamma_ij[i][k][l]);
                    Gamma_lik += metric.gamma_inv[j][l] * term;
                }
                Gamma_j_ik[j][i][k] = Gamma_lik;
            }
        }
    }

    // =========================================================================
    // STEP 4: Compute source term S_i = α√γ Γ^α_{iβ} T_α^β
    // =========================================================================
    real source_i[3];
    for (int i = 0; i < 3; ++i) {
        real S = 0.0;

        // Contribution from Γ^0_{i0} T_0^0
        S += Gamma_0_i0[i] * T_mix_0_0;

        // Contribution from Γ^0_{ik} T_0^k = 0 (Γ^0_{ik} = 0 for static metric)

        // Contribution from Γ^j_{i0} T_j^0 = 0 (Γ^j_{i0} = 0 for static metric)

        // Contribution from Γ^j_{ik} T_j^k
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                S += Gamma_j_ik[j][i][k] * T_mix_j_k[j][k];
            }
        }

        // Multiply by α√γ
        source_i[i] = alpha * metric.sqrt_gamma * S;
    }

    // =========================================================================
    // STEP 5: Convert to per-baryon acceleration for velocity evolution
    // =========================================================================
    // From GRGSPH formulation §5.3:
    //   f_i = (√-g / 2N*) T^μν ∂g_μν/∂x^i
    //
    // The Valencia form gives momentum source S_i (= source_i computed above).
    // To convert to per-baryon acceleration:
    //   f^k = γ^{ki} source_i / (N* h W²)
    //
    // where the factor h W² comes from the momentum-velocity relation:
    //   p_i = n h W² v_i (in natural units)
    //
    // Note: N* = √γ n Γ = √γ n W

    for (int k = 0; k < DIM; ++k) {
        real f_k = 0.0;
        for (int i = 0; i < 3; ++i) {
            f_k += metric.gamma_inv[k][i] * source_i[i];
        }
        // Convert Valencia source to per-baryon acceleration
        // Factor α/(N* h W²) accounts for:
        // - α: coordinate time vs proper time
        // - N*: conserved baryon density
        // - h W²: relativistic momentum-velocity relation
        f_source[k] = f_k * alpha / (N_star * h * W2);
    }

    return f_source;
}

/**
 * Compute gravitational work (energy source) for dynamical spacetimes
 *
 * From GRGSPH Formulation §5.3:
 *   Λ = -(√-g / 2N*) T^μν ∂g_μν/∂t
 *
 * This term represents the energy exchange between the fluid and the
 * dynamically evolving gravitational field.
 *
 * Using the stress-energy tensor decomposition:
 *   T^μν = ρh W² u^μ u^ν + P g^μν
 *
 * where u^μ = (1, v^i)/α in coordinate basis, we expand:
 *
 *   T^μν ∂g_μν/∂t = T^00 ∂g_00/∂t + 2 T^0i ∂g_0i/∂t + T^ij ∂g_ij/∂t
 *
 * For the 3+1 decomposition with g_00 = -α², g_0i = β_i, g_ij = γ_ij:
 *   ∂g_00/∂t = -2α ∂α/∂t
 *   ∂g_0i/∂t = ∂β_i/∂t
 *   ∂g_ij/∂t = ∂γ_ij/∂t
 *
 * Examples of time-dependent metrics:
 *   - Cosmological FLRW: a(t) expansion, γ_ij = a²(t) δ_ij
 *   - Gravitational waves: γ_ij = δ_ij + h_ij(t,x)
 *   - Binary BH inspiral: time-dependent M, a, or numerical metric
 */
real GRFluidForce::compute_gravitational_work(
    const GRPrimitive& prim,
    const Metric31& metric,
    const MetricTimeDerivatives& time_derivs,
    real N_star) const
{
    // For stationary metrics or negligible density, no work
    if (N_star < 1e-15) return 0.0;

    const real n = prim.n;
    const real P = prim.P;
    const real h = prim.h;
    const real W = prim.Gamma;  // Lorentz factor
    const real W2 = W * W;
    const real alpha = metric.alpha;
    const real alpha2 = alpha * alpha;

    // =========================================================================
    // STEP 1: Compute contravariant stress-energy tensor T^μν
    // =========================================================================
    const real rho_h = n * h;

    // T^00 = ρh W²/α² - P/α²
    //      = (ρh W² - P)/α²
    const real T_up_00 = (rho_h * W2 - P) / alpha2;

    // T^0i = ρh W² V^i / α
    real T_up_0i[3];
    for (int i = 0; i < 3; ++i) {
        T_up_0i[i] = rho_h * W2 * prim.V[i] / alpha;
    }

    // T^ij = ρh W² V^i V^j + P γ^ij
    real T_up_ij[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            T_up_ij[i][j] = rho_h * W2 * prim.V[i] * prim.V[j] + P * metric.gamma_inv[i][j];
        }
    }

    // =========================================================================
    // STEP 2: Compute T^μν ∂g_μν/∂t
    // =========================================================================
    // Using g_00 = -α², so dg_tt_dt = ∂(-α²)/∂t
    real T_dg_dt = 0.0;

    // T^00 ∂g_00/∂t term
    T_dg_dt += T_up_00 * time_derivs.dg_tt_dt;

    // 2 T^0i ∂g_0i/∂t term (factor of 2 from symmetry)
    // g_0i = β_i (lowered shift), so dg_ti_dt[i] = ∂β_i/∂t
    for (int i = 0; i < 3; ++i) {
        T_dg_dt += 2.0 * T_up_0i[i] * time_derivs.dg_ti_dt[i];
    }

    // T^ij ∂γ_ij/∂t term
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            T_dg_dt += T_up_ij[i][j] * time_derivs.dgamma_ij_dt[i][j];
        }
    }

    // =========================================================================
    // STEP 3: Compute gravitational work Λ = -(√-g / 2N*) T^μν ∂g_μν/∂t
    // =========================================================================
    const real Lambda = -metric.sqrt_neg_g * T_dg_dt / (2.0 * N_star);

    return Lambda;
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
#ifndef EXHAUSTIVE_SEARCH
    auto* tree = sim->get_tree().get();
#endif

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

        // In geodesic mode, skip pairwise SPH forces (only use gravitational source)
        // This is for testing geodesic motion without hydrodynamic effects
        if (!m_geodesic_mode) {
        for (int n = 0; n < n_neighbor; ++n) {
            const int j = neighbor_list[n];
            if (j == i) continue;

            auto& p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(p_i.pos, p_j.pos);
            const real r = std::abs(r_ij);

            if (r > 3.0 * sqrt_two * std::max(p_i.sml, p_j.sml)) continue;

            // === GR-GSPH FORMULATION (Kitajima + Liptai & Price) ===
            // Following formulation document §8.1-8.3

            // Compute metric at particle j position
            Metric31 metric_j;
            if (m_metric) {
                m_metric->compute(p_j.pos, metric_j);
            }

            // Unit vector from j to i (line-of-sight direction)
            // Note: Use -r_ij/r to match SR-GSPH convention where n_ij points from j to i
            // This ensures consistent "left" (i) and "right" (j) state identification
            const vec_t n_ij = r_ij * (-1.0 / r);

            // Project coordinate velocities onto line-of-sight
            // v_n = v · n̂ (simple dot product in coordinate space)
            const real v_n_L = inner_product(p_i.vel, n_ij);
            const real v_n_R = inner_product(p_j.vel, n_ij);

            // Tangent velocity = total velocity minus normal component
            const vec_t v_t_L_vec = p_i.vel - n_ij * v_n_L;
            const vec_t v_t_R_vec = p_j.vel - n_ij * v_n_R;
            const real v_t_L = std::abs(v_t_L_vec);
            const real v_t_R = std::abs(v_t_R_vec);

            // Primitive variables
            const real p_L = p_i.pres;
            const real rho_L = p_i.dens;
            const real p_R = p_j.pres;
            const real rho_R = p_j.dens;

            // Solve Riemann problem with coordinate velocities
            // The Riemann solver handles the relativistic Lorentz factor internally
            real P_star = 0.0;
            real v_n_star = 0.0;
            real v_t_star = 0.0;

            if (p_j.is_ghost) {
                // Wall boundary condition
                P_star = p_L;
                v_n_star = 0.0;
                v_t_star = v_t_L;
            } else {
                // State arrays: [v_n, rho, P, c_s, v_t]
                real left_state[5] = {v_n_L, rho_L, p_L, p_i.sound, v_t_L};
                real right_state[5] = {v_n_R, rho_R, p_R, p_j.sound, v_t_R};

                solve_interface_state(left_state, right_state, P_star, v_n_star, v_t_star);
            }

            // Reconstruct star velocity in coordinate frame
            // Use upwind tangent direction
            vec_t v_star_vec;
            if (v_n_star > 0.0 && v_t_L > 1e-15) {
                // Flow from i to j: use left tangent direction
                v_star_vec = n_ij * v_n_star + v_t_L_vec * (v_t_star / v_t_L);
            } else if (v_n_star <= 0.0 && v_t_R > 1e-15) {
                // Flow from j to i: use right tangent direction
                v_star_vec = n_ij * v_n_star + v_t_R_vec * (v_t_star / v_t_R);
            } else {
                v_star_vec = n_ij * v_n_star;
            }

            // Volume elements
            const real V_i = p_i.nu / p_i.N;
            const real V_j = p_j.nu / p_j.N;

            // === GRGSPH Formulation §8.3: V²_ij interpolation ===
            // V²_ij = (1/2)(V²_ij(h_i) + V²_ij(h_j))
            // For SPH approximation: use symmetric average of volume products
            const real V2_ij_i = V_i * V_i;  // V² at h_i
            const real V2_ij_j = V_j * V_j;  // V² at h_j
            // Note: V2_ij_interp available for symmetric formulation
            (void)(0.5 * (V2_ij_i + V2_ij_j));  // Suppress unused warning

            // === GRGSPH Formulation §8.3: √(-g) averaging ===
            // √(-g)_ij = (1/2)(√(-g_i) + √(-g_j))
            const real sqrt_neg_g_avg = 0.5 * (metric_i.sqrt_neg_g + metric_j.sqrt_neg_g);

            // Kernel gradients with √2 h argument (formulation §8.1)
            const vec_t grad_W_i = kernel->dw(r_ij, r, sqrt_two * p_i.sml);
            const vec_t grad_W_j = kernel->dw(r_ij, r, sqrt_two * p_j.sml);

            // Grad-h corrections
            const real omega_i = p_i.gradh;
            const real omega_j = p_j.gradh;

            // === GRGSPH Momentum Equation (Formulation §8.1) ===
            // dS_k/dt = -P*_ij √(-g)_ij V²_ij [∇_k W_ij(√2 h_i) - ∇_k W_ij(√2 h_j)]
            // Note: grad_W_i and grad_W_j have opposite signs due to antisymmetry
            // Using grad-h corrected form with separate V² terms for stability
            const real coef_i = -P_star * V2_ij_i * omega_i * sqrt_neg_g_avg;
            const real coef_j = -P_star * V2_ij_j * omega_j * sqrt_neg_g_avg;
            const vec_t force = grad_W_i * coef_i + grad_W_j * coef_j;

            // === GRGSPH Energy Equation (Formulation §8.2) ===
            // de/dt = -P*_ij v*^k_ij √(-g)_ij V²_ij [∇_k W_ij(√2 h_i) - ∇_k W_ij(√2 h_j)]
            // Power: dE/dt = v* · F (with proper GR weighting already in force)
            const real power = inner_product(v_star_vec, force);

            // Accumulate
            for (int d = 0; d < DIM; ++d) {
                dS_acc[d] += force[d];
            }
            de_acc += power;
        }
        } // End of if (!m_geodesic_mode)

        // Add gravitational source terms using Valencia formulation
        // Momentum: f_i = (√-g / 2N*) T^μν ∂g_μν/∂x^i (formulation §5.3)
        // Energy:   Λ = -(√-g / 2N*) T^μν ∂g_μν/∂t   (formulation §5.3)
        if (m_metric) {
            GRPrimitive prim;
            prim.n = p_i.dens;
            prim.P = p_i.pres;
            prim.h = 1.0 + m_gamma / (m_gamma - 1.0) * prim.P / std::max(prim.n, 1e-15);

            // Transform coordinate velocity to Eulerian velocity
            // V^i = (v^i + β^i) / α
            real v_coord[3] = {0.0, 0.0, 0.0};
            for (int d = 0; d < DIM; ++d) {
                v_coord[d] = p_i.vel[d];
            }
            real V_eul[3];
            metric_i.coord_to_eulerian(v_coord, V_eul);

            // Copy to prim structure
            for (int d = 0; d < 3; ++d) {
                prim.V[d] = V_eul[d];
            }

            // Compute proper Lorentz factor using spatial metric
            prim.Gamma = metric_i.lorentz_factor(V_eul);

            // Compute N* = n Γ √γ (conserved baryon number density, formulation §4.1)
            const real N_star = prim.n * prim.Gamma * metric_i.sqrt_gamma;

            // Momentum source: f_i = (√-g / 2N*) T^μν ∂g_μν/∂x^i
            // Using Valencia formulation: S_i = α√γ Γ^α_{iβ} T_α^β
            vec_t f_grav = compute_metric_source(prim, metric_i, derivs_i, N_star);

            for (int d = 0; d < DIM; ++d) {
                dS_acc[d] += f_grav[d];
            }

            // Energy source: Λ = -(√-g / 2N*) T^μν ∂g_μν/∂t (formulation §5.3)
            // For stationary metrics (Schwarzschild, Kerr), ∂g_μν/∂t = 0, so Λ = 0
            // For time-dependent metrics, compute full gravitational work
            real Lambda = 0.0;
            if (!m_metric->is_stationary()) {
                MetricTimeDerivatives time_derivs;
                m_metric->compute_time_derivatives(p_i.pos, 0.0, time_derivs);  // TODO: pass actual time
                Lambda = compute_gravitational_work(prim, metric_i, time_derivs, N_star);
            }
            de_acc += Lambda;
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
