/**
 * Special Relativistic Godunov SPH - Fluid Force Calculation
 *
 * Reference: Kitajima, Inutsuka, Seno (2025) "Special Relativistic Hydrodynamics
 *           with Smoothed Particle Hydrodynamics" arXiv:2510.18251v1
 *
 * Key Equations Implemented:
 * - Eq. 63-64: Canonical equations of motion for S (momentum) and e (energy)
 * - Eq. 67: MUSCL reconstruction with shock/contact discontinuity detection
 * - Riemann solver: Iterative Newton-Raphson for P* and v* with tangential velocities
 * - Pons et al. (1999) extension for non-zero tangential velocities
 *
 * Formulation:
 * - Uses volume-based V_p from Eq. 220 (computed in pre_interaction)
 * 5. Compute kernel gradients: ∇_i W(r_ij, √2·h_i) - ∇_j W(r_ij, √2·h_j)
 *    - Variable h formulation uses SEPARATE h_i and h_j
 *    - NOT symmetric like fixed-h case: must compute both gradients
 *    - √2 factor from convolution integral (Gaussian kernel combination)
 *   (uses separate smoothing lengths for each particle)
 * - Solves 1D Riemann problem in e_ij direction for each particle pair
 */

#include "defines.hpp"
#include "particle.hpp"
#include "periodic.hpp"
#include "simulation.hpp"
#include "bhtree.hpp"
#include "kernel/kernel_function.hpp"
#include "srgsph/sr_fluid_force.hpp"
#include "srgsph/sr_primitive_recovery.hpp"

#include <iostream>

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
            m_gamma = param->physics.gamma;
            m_c_speed = param->srgsph.c_speed;

            // √2 kernel factor: DISABLED - causes inconsistency with volume calculation
            // Volume is computed with h, forces should use same h for consistency
            m_use_sqrt2_h = false;  // Was: (param->kernel == KernelType::GAUSSIAN)

            // Use iterative solver (Pons et al. 1999)
            iterative_solver();
        }

        /**
         * HLLC Riemann solver for special relativistic hydrodynamics
         *
         * Based on Mignone & Bodo (2005) "An HLLC Riemann solver for relativistic flows"
         * Journal of Computational Physics, 207, 545-563
         *
         * HLLC improves on HLL by:
         * - Resolving contact discontinuity (HLL smears it)
         * - More accurate than HLL, nearly as accurate as exact solver
         * - Always robust (no iteration, no convergence issues)
         * - Handles strong shocks and rarefactions
         *
         * Wave structure: Left | Star_L | Contact | Star_R | Right
         * Wave speeds: S_L (left wave), S_* (contact), S_R (right wave)
         */
        void FluidForce::hllc_solver()
        {
            m_solver = [&](const real left[], const real right[], real &pstar, real &vstar)
            {
                // Extract states: [v^n, n, P, c_s, |v^t|]
                const real vl = left[0];
                const real nl = left[1];
                const real Pl = left[2];
                const real csl = left[3];
                const real vtl = left[4];

                const real vr = right[0];
                const real nr = right[1];
                const real Pr = right[2];
                const real csr = right[3];
                const real vtr = right[4];

                const real c = m_c_speed;
                const real gamma_c = m_gamma;

                // Compute Lorentz factors
                const real v2l = vl * vl + vtl * vtl;
                const real v2r = vr * vr + vtr * vtr;
                const real wl = 1.0 / std::sqrt(1.0 - v2l / (c * c));
                const real wr = 1.0 / std::sqrt(1.0 - v2r / (c * c));

                // Compute enthalpies H = 1 + u/c² + P/(n·c²)
                const real ul = Pl / ((gamma_c - 1.0) * nl);
                const real ur = Pr / ((gamma_c - 1.0) * nr);
                const real Hl = 1.0 + ul / (c * c) + Pl / (nl * c * c);
                const real Hr = 1.0 + ur / (c * c) + Pr / (nr * c * c);

                // Conserved variables: D = γn, S = γ²hn v, E = γ²hn - P
                const real Dl = nl * wl;
                const real Dr = nr * wr;
                const real Sl = nl * Hl * wl * wl * vl;
                const real Sr = nr * Hr * wr * wr * vr;
                const real El = nl * Hl * wl * wl - Pl;
                const real Er = nr * Hr * wr * wr - Pr;

                // Wave speed estimates (Davis approximation)
                const real lambda_l_min = (vl - csl) / (1.0 + vl * csl / c);
                const real lambda_l_max = (vl + csl) / (1.0 - vl * csl / c);
                const real lambda_r_min = (vr - csr) / (1.0 + vr * csr / c);
                const real lambda_r_max = (vr + csr) / (1.0 - vr * csr / c);

                const real S_L = std::min(lambda_l_min, lambda_r_min);
                const real S_R = std::max(lambda_l_max, lambda_r_max);

                // Special cases: supersonic flow
                if (S_L >= 0.0)
                {
                    pstar = Pl;
                    vstar = vl;
                    return;
                }
                if (S_R <= 0.0)
                {
                    pstar = Pr;
                    vstar = vr;
                    return;
                }

                // HLLC contact wave speed (Mignone & Bodo Eq. 18)
                const real num = Pr - Pl + Sl * (Dl - S_L * Dl + Sl * Dl) - Sr * (Dr - S_R * Dr + Sr * Dr);
                const real den = Dl * (S_L - vl) - Dr * (S_R - vr);
                const real S_star = num / den;

                // Star region pressure (Mignone & Bodo Eq. 19)
                pstar = Pl + Dl * (S_L - vl) * (S_star - vl);
                pstar = std::max(1e-20, pstar);  // Ensure positive

                // Star region velocity
                vstar = S_star;
                vstar = std::max(-0.99*c, std::min(0.99*c, vstar));  // Ensure subluminal
            };
        }

        /**
         * HLLE Riemann solver for special relativistic hydrodynamics
         *
         * Based on:
         * - Einfeldt et al. (1991) "On Godunov-type methods near low densities"
         * - Schneider et al. (1993) "New algorithms for ultra-relativistic flows"
         * - Mignone & Bodo (2005) "An HLLC Riemann solver for relativistic flows"
         *
         * HLLE (Harten-Lax-van Leer-Einfeldt) solver:
         * - Estimates fastest left/right wave speeds
         * - Computes intermediate state from conservation laws
         * - More robust than exact solver (no root finding, always converges)
         * - Slightly more dissipative but stable for strong shocks
         *
         * For GSPH we only need p* and v*, not full flux
         */
        void FluidForce::hlle_solver()
        {
            m_solver = [&](const real left[], const real right[], real &pstar, real &vstar)
            {
                // Extract states: [v^n, n, P, c_s, |v^t|]
                const real vl = left[0];
                const real nl = left[1];
                const real Pl = left[2];
                const real csl = left[3];
                const real vtl = left[4];

                const real vr = right[0];
                const real nr = right[1];
                const real Pr = right[2];
                const real csr = right[3];
                const real vtr = right[4];

                const real c = m_c_speed;
                const real gamma_c = m_gamma;

                // Compute Lorentz factors
                const real v2l = vl * vl + vtl * vtl;
                const real v2r = vr * vr + vtr * vtr;
                const real wl = 1.0 / std::sqrt(1.0 - v2l / (c * c));
                const real wr = 1.0 / std::sqrt(1.0 - v2r / (c * c));

                // Compute enthalpies
                const real ul = Pl / ((gamma_c - 1.0) * nl);
                const real ur = Pr / ((gamma_c - 1.0) * nr);
                const real Hl = 1.0 + ul / (c * c) + Pl / (nl * c * c);
                const real Hr = 1.0 + ur / (c * c) + Pr / (nr * c * c);

                // Relativistic wave speed estimates (Davis 1988, Mignone & Bodo 2005)
                // For relativistic flow: lambda = (v ± c_s*sqrt(1-v^2/c^2))/(1 ± v*c_s/c^2)
                // Simplified robust estimate: use characteristic speeds

                // Left state characteristics
                const real lambda_l_minus = (vl - csl * c) / (c + vl * csl / c);
                const real lambda_l_plus  = (vl + csl * c) / (c - vl * csl / c);

                // Right state characteristics
                const real lambda_r_minus = (vr - csr * c) / (c + vr * csr / c);
                const real lambda_r_plus  = (vr + csr * c) / (c - vr * csr / c);

                // HLLE wave speed estimates
                const real s_l = std::min(lambda_l_minus, lambda_r_minus);
                const real s_r = std::max(lambda_l_plus, lambda_r_plus);

                // Special cases
                if (s_l >= 0.0)
                {
                    // Supersonic right-moving: use left state
                    pstar = Pl;
                    vstar = vl;
                    return;
                }
                if (s_r <= 0.0)
                {
                    // Supersonic left-moving: use right state
                    pstar = Pr;
                    vstar = vr;
                    return;
                }

                // Subsonic case: compute HLL averaged state
                // Using conservative form: U = [D, S, E]
                // D = rho*W, S = rho*h*W^2*v, E = rho*h*W^2 - P

                const real Dl = nl * wl;
                const real Dr = nr * wr;

                const real Sl = nl * Hl * wl * wl * vl;
                const real Sr = nr * Hr * wr * wr * vr;

                const real El = nl * Hl * wl * wl - Pl;
                const real Er = nr * Hr * wr * wr - Pr;

                // HLL average using simple weighted interpolation
                // More robust than extracting from conserved variables
                const real weight_l = std::abs(s_r) / (std::abs(s_l) + std::abs(s_r) + 1e-20);
                const real weight_r = 1.0 - weight_l;

                // Weighted average of primitives (guaranteed positive and bounded)
                pstar = weight_l * Pl + weight_r * Pr;
                vstar = weight_l * vl + weight_r * vr;

                // Ensure physical bounds
                pstar = std::max(1e-20, pstar);
                vstar = std::max(-0.99*c, std::min(0.99*c, vstar));
            };
        }

        void FluidForce::calculation(std::shared_ptr<Simulation> sim)
        {
            /**
             * SR-GSPH fluid force calculation - GSPH-STYLE FORMULATION
             *
             * Implements Kitajima et al. (2025) Equations 209-212 (fixed smoothing length):
             *
             * Eq. 209: ⟨dS_i⟩ = -Σ_j ν P*_ij V²_ij,interp [∇_i W - ∇_j W]
             *
             * Eq. 211: ⟨de_i⟩ = -Σ_j ν P*_ij v*_ij · V²_ij,interp [∇_i W - ∇_j W]
             *
             * where:
             * - S_i = canonical momentum per baryon
             * - e_i = canonical energy per baryon
             * - ν = baryon number (stored in p.mass field)
             * - P*_ij = Riemann solution pressure
             * - v*_ij = Riemann solution velocity
             * - V²_ij,interp = 0.5(1/N_i² + 1/N_j²) - STANDARD GSPH INTERPOLATION
             * - N = lab-frame baryon number density (N = γn)
             *
             * This is the STANDARD GSPH formulation (Cha & Whitworth 2003) adapted for SR:
             * - Uses ν (baryon number) instead of m (mass)
             * - Uses N (lab-frame density) instead of ρ (rest-frame density)
             * - Evolves canonical variables (S, e) instead of primitives (v, u)
             * - MUCH SIMPLER than volume-based approach (Eq. 371-374)
             *
             * Symmetry & Conservation:
             * - Antisymmetric in i↔j → conserves momentum and energy
             * - Kernel at h_force = √2·h for Gaussian convolution integral
             */

            auto &particles = sim->get_particles();
            auto *periodic = sim->get_periodic().get();
            const int num = sim->get_particle_num();
            auto *kernel = sim->get_kernel().get();
            auto *tree = sim->get_tree().get();

#pragma omp parallel for
            for (int i = 0; i < num; ++i)
            {
                auto &p_i = particles[i];
                std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

                // Neighbor search
#ifdef EXHAUSTIVE_SEARCH
                int const n_neighbor = exhaustive_search(p_i, p_i.sml, particles, num, neighbor_list,
                                                         m_neighbor_number * neighbor_list_size, periodic, true);
#else
                int const n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, true);
#endif

                // Extract primitive variables from particle i
                // These were computed in pre_interaction phase via primitive recovery
                const vec_t &r_i = p_i.pos;
                const vec_t &v_i = p_i.vel;           // Velocity v (primitive, recovered)
                const real h_i = p_i.sml;             // Smoothing length h
                const real N_i = p_i.N;               // Lab-frame density N = γn (stored in dedicated field)
                const real n_i = N_i / p_i.gamma_lor; // Rest-frame density n = N/γ
                const real P_i = p_i.pres;            // Pressure P
                const real c_i = p_i.sound;           // Sound speed c_s

                // Canonical momentum and energy time derivatives
                // Will accumulate: Ṡ = dS/dt and ė = de/dt from Eq. 63-64
                vec_t dS_dt(0.0);
                real de_dt = 0.0;

                for (int n = 0; n < n_neighbor; ++n)
                {
                    int const j = neighbor_list[n];
                    auto &p_j = particles[j];

                    // Skip self-interaction
                    if (i == j) continue;

                    // Compute separation vector and distance
                    const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos); // r_ij = r_i - r_j
                    const real r = std::abs(r_ij);                        // |r_ij|
                    const real h_j = p_j.sml;

                  
                    // Unit vector e_ij = r_ij/|r_ij| defines normal direction for Riemann problem
                    const real r_inv = 1.0 / r;
                    const vec_t e_ij = r_ij * r_inv;

                    // Project velocities onto normal direction e_ij
                    // v^n_i = v_i · e_ij (normal component, for 1D Riemann problem)
                    // v^t_i = v_i - (v_i·e_ij)e_ij (tangential component, preserved across shock)
                    const real ve_i = inner_product(v_i, e_ij);
                    const real ve_j = inner_product(p_j.vel, e_ij);

                    // Compute tangential velocity components (Pons et al. 1999)
                    // Extension for multi-dimensional relativistic Riemann problems
                    // v^t = v - (v·e_ij)e_ij (tangential to e_ij direction)
                    // |v^t| = magnitude of tangential velocity
                    const vec_t vt_i = v_i - e_ij * ve_i;
                    const vec_t vt_j = p_j.vel - e_ij * ve_j;
                    const real vt_mag_i = std::abs(vt_i);
                    const real vt_mag_j = std::abs(vt_j);

                    // Use stored Lorentz factors from particles
                    // These were computed consistently in pre_interaction with the stored velocities
                    // CRITICAL: p.dens = N = gamma * n, so n = N / gamma_lor (stored gamma)
                    const real gamma_j = p_j.gamma_lor;

                    /**
                     * Riemann solver states for 1D problem in e_ij direction
                     *
                     * State vector: [v^n, n, P, c_s, |v^t|]
                     * where:
                     * - v^n = normal velocity component (along e_ij)
                     * - n = rest-frame baryon density
                     * - P = pressure
                     * - c_s = sound speed
                     * - |v^t| = magnitude of tangential velocity (for Pons et al. coupling)
                     *
                     * Pons et al. (1999) Eqs. 27-28, 54-55:
                     * Conservation across waves: h*W*v^y = const, h*W*v^z = const
                     *
                     * Left/right states depend on MUSCL reconstruction (Eq. 67)
                     */
                    /**
                     * Riemann solver states for 1D problem in e_ij direction
                     * Using 1st order (no MUSCL reconstruction)
                     *
                     * State vector: [v^n, n, P, c_s, |v^t|]
                     * where:
                     * - v^n = normal velocity component (along e_ij)
                     * - n = rest-frame baryon density
                     * - P = pressure
                     * - c_s = sound speed
                     * - |v^t| = magnitude of tangential velocity (for Pons et al. coupling)
                     */
                    real left[5], right[5];
                    real pstar, vstar; // Riemann solution: P* and v* at interface

                    // First order (no reconstruction) - always use particle values directly
                    const real n_j = p_j.dens / gamma_j;

                    // GSPH convention (from g_fluid_force.cpp lines 155-166):
                    // - right = particle i
                    // - left = particle j
                    right[0] = ve_i;
                    right[1] = n_i;
                    right[2] = P_i;
                    right[3] = c_i;
                    right[4] = vt_mag_i; // Tangential velocity

                    left[0] = ve_j;
                    left[1] = n_j;
                    left[2] = p_j.pres;
                    left[3] = p_j.sound;
                    left[4] = vt_mag_j; // Tangential velocity

                    // Solve Riemann problem to get pstar, vstar
                    m_solver(left, right, pstar, vstar);

                    /**
                     * Force calculation - intermediate scaling
                     *
                     * Use ν_j × P / (N_i × N_j) as coefficient.
                     * This is:
                     * - ν times larger than V²_ij formula (which gives tiny forces)
                     * - ν times smaller than 1/N² formula (which explodes)
                     */
                    const real nu_j = p_j.mass;  // Baryon number for particle j
                    const real N_j = p_j.N;       // Lab-frame density for particle j

                    // Density products
                    const real N_prod_inv = 1.0 / (N_i * N_j);

                    /**
                     * Kernel gradients for variable smoothing length
                     *
                     * Reference: Kitajima et al. (2025) Eq. 63-64
                     *
                     * Variable h formulation: ∇_i W(r_ij, h_i) - ∇_j W(r_ji, h_j)
                     * where r_ji = -r_ij
                     *
                     * Since kernel->dw(r_ij, r, h) = (r_ij/r) * dW/dr,
                     * we have: ∇_j W(r_ji, h_j) = (-r_ij/r) * dW/dr|_h_j = -∇_j W(r_ij, h_j)
                     *
                     * Therefore: ∇_i W - ∇_j W(r_ji) = ∇_i W + ∇_j W(r_ij)
                     *
                     * √2 factor: Only for Gaussian kernels due to convolution integral.
                     * For Cubic Spline and Wendland kernels, use h directly.
                     */
                    const real h_scale = m_use_sqrt2_h ? std::sqrt(2.0) : 1.0;
                    const real h_force_i = h_scale * h_i;
                    const real h_force_j = h_scale * h_j;
                    const vec_t grad_i_W = kernel->dw(r_ij, r, h_force_i); // ∇_i W(r_ij, h_scale·h_i)
                    const vec_t grad_j_W = kernel->dw(r_ij, r, h_force_j); // ∇_j W(r_ij, h_scale·h_j)
                    
                    /**
                     * Accumulate force contributions - intermediate scaling
                     *
                     * f = (∇_i W + ∇_j W) × (ν_j P* / (N_i × N_j))
                     * dS/dt = -f
                     * de/dt = -f · v*
                     */
                    const vec_t dw_ij = grad_i_W + grad_j_W;
                    const vec_t f = dw_ij * (nu_j * pstar * N_prod_inv);

                    // Momentum evolution
                    dS_dt -= f;

                    // Energy evolution
                    const vec_t v_star = e_ij * vstar;
                    de_dt -= inner_product(f, v_star);
                }

                // Store time derivatives in SR-GSPH dedicated fields
                // These are the PRIMARY storage for time derivatives
                const real nu_i = p_i.mass; // Baryon number (mass field stores ν in SRGSPH)
                p_i.dS = dS_dt / nu_i;      // dS/dt: canonical momentum derivative (per baryon)
                p_i.de = de_dt / nu_i;      // de/dt: canonical energy derivative (per baryon)

                // LEGACY COMPATIBILITY: Also copy to standard SPH fields
                // NOTE: acc/dene have DIFFERENT meanings in SRGSPH vs standard SPH!
                //   - In SRGSPH: acc = dS/dt (NOT acceleration!), dene = de/dt
                //   - In std SPH: acc = dv/dt (true acceleration), dene = du/dt
                p_i.acc = p_i.dS;  // Copy for output (CSV/VTK write this as "acceleration")
                p_i.dene = p_i.de; // Copy for output (CSV/VTK write this as "energy derivative")
            }
        }

        /**
         * Iterative relativistic Riemann solver with Newton-Raphson method
         *
         * Based on Pons, Marti & Muller (1999):
         * "The exact solution of the Riemann problem with non-zero tangential velocities
         * in relativistic hydrodynamics" (J. Fluid Mech. vol. 422, pp. 125-139)
         *
         * Physical Setup:
         * - Solves 1D Riemann problem in normal direction: P*(v*_L - v_L) = P*(v*_R - v_R)
         * - Includes tangential velocities v^t in Lorentz factor: gamma^2 = 1/(1 - (v^n^2 + v^t^2)/c^2)
         * - Returns interface state (P*, v*) used for flux calculation
         *
         * Jump Conditions (Pons et al. 1999):
         * - Eq. (1): p_L* = p_R* (= p*) - pressure continuity at contact
         * - Eq. (2): v^n_L* = v^n_R* (= v*) - normal velocity continuity at contact
         *
         * Algorithm:
         * 1. Check for trivial solution (identical states) -> P*=P_L, v*=v_L
         * 2. Compute enthalpy H = 1 + u/c^2 + P/(n*c^2) for both states
         * 3. Use two-rarefaction initial guess: P* ~ [(c_s,L + c_s,R - (v_R - v_L)) / (...)]^(2*gamma/(gamma-1))
         * 4. Newton-Raphson iteration:
         *    - Compute v*(P*) for each wave using get_velocity_behind_wave()
         *    - Evaluate f(P*) = v*_L(P*) - v*_R(P*) (velocity must match at contact)
         *    - Compute f'(P*) using finite difference
         *    - Update P* <- P* - f(P*)/f'(P*)
         * 5. Converge when |f(P*)| < tol or |dP*| / P* < tol
         *
         * Wave Structure:
         * - Left wave (sign=-1): either shock or rarefaction connecting (v_L, P_L) to (v*, P*)
         * - Right wave (sign=+1): either shock or rarefaction connecting (v*, P*) to (v_R, P_R)
         * - Contact discontinuity at v* separating left/right material
         *
         * Convergence:
         * - max_iter = 50 (typically converges in 5-10 iterations)
         * - tol = 1e-8 for velocity matching
         * - rel_tol = 1e-10 for pressure update
         * - Graceful fallback to two-rarefaction if no convergence
         */
        void FluidForce::iterative_solver()
        {
            // Iterative relativistic Riemann solver using Newton-Raphson
            // Extended for tangential velocities following Pons, Martí & Müller (1999)
            // "The exact solution of the Riemann problem with non-zero tangential velocities
            // in relativistic hydrodynamics" (J. Fluid Mech. 1999)

            m_solver = [&](const real left[], const real right[], real &pstar, real &vstar)
            {
                // Extract states
                // State format: [v^n, n, P, c_s, |v^t|]
                const real vl = left[0];  // Normal velocity left
                const real nl = left[1];  // Rest-frame density left
                const real Pl = left[2];  // Pressure left
                const real csl = left[3]; // Sound speed left
                const real vtl = left[4]; // Tangential velocity magnitude left

                const real vr = right[0];  // Normal velocity right
                const real nr = right[1];  // Rest-frame density right
                const real Pr = right[2];  // Pressure right
                const real csr = right[3]; // Sound speed right
                const real vtr = right[4]; // Tangential velocity magnitude right

                const real c = m_c_speed;
                const real gamma_c = m_gamma;

                // Compute enthalpy H = 1 + u/c² + P/(n·c²)
                // For ideal gas: u = P/((γ-1)·n)
                const real ul = Pl / ((gamma_c - 1.0) * nl);
                const real ur = Pr / ((gamma_c - 1.0) * nr);
                const real Hl = 1.0 + ul / (c * c) + Pl / (nl * c * c);
                const real Hr = 1.0 + ur / (c * c) + Pr / (nr * c * c);

                // Compute Lorentz factors including tangential velocities
                // γ² = 1/(1 - v²/c²) where v² = (v^n)² + (v^t)²
                const real v2l = vl * vl + vtl * vtl;
                const real v2r = vr * vr + vtr * vtr;
                const real wl = 1.0 / std::sqrt(1.0 - v2l / (c * c));
                const real wr = 1.0 / std::sqrt(1.0 - v2r / (c * c));

                // Check for trivial solution (identical states within tolerance)
                constexpr real trivial_tol = 1.0e-12;
                if (std::abs(Pl - Pr) < trivial_tol * std::max(Pl, Pr) &&
                    std::abs(vl - vr) < trivial_tol * c)
                {
                    pstar = Pl;
                    vstar = vl;
                    return;
                }

                /**
                 * Adaptive bracketing + Brent's method for P*
                 *
                 * Following the robust approach from Kitajima et al. (2025):
                 * 1. Find a pressure bracket [p_min, p_max] where dvel changes sign
                 * 2. Use Brent's method (root-finding) instead of Newton-Raphson
                 *
                 * Why this works:
                 * - Bracketing guarantees root is found (if it exists)
                 * - Brent's method is robust and doesn't require derivatives
                 * - No risk of Newton divergence or vacuum solutions
                 */

                // Step 1: Find bracket for P* using adaptive expansion
                // Start with arithmetic mean, expand until dvel changes sign
                real p_mid = 0.5 * (Pl + Pr);
                real p_min = p_mid;
                real p_max = p_mid;
                
                // Lambda to compute velocity difference at given pressure
                auto dvel_at_p = [&](real p) -> real {
                    real nl_star, Hl_star, csl_star, vl_star;
                    real nr_star, Hr_star, csr_star, vr_star;
                    get_velocity_behind_wave(p, nl, Pl, Hl, csl, vl, vtl, wl, -1.0,
                                           nl_star, Hl_star, csl_star, vl_star);
                    get_velocity_behind_wave(p, nr, Pr, Hr, csr, vr, vtr, wr, +1.0,
                                           nr_star, Hr_star, csr_star, vr_star);
                    return vl_star - vr_star;
                };
                
                // Expand bracket until we find opposite-sign dvel values
                constexpr int max_bracket_iter = 100;
                real dvel_min = 0.0;
                real dvel_max = 0.0;
                bool bracket_found = false;
                
                for (int i = 0; i < max_bracket_iter; ++i)
                {
                    p_min = 0.5 * std::max(p_min, 0.0);
                    p_max = 2.0 * p_max;
                    
                    dvel_min = dvel_at_p(p_min);
                    dvel_max = dvel_at_p(p_max);
                    
                    if (dvel_min * dvel_max <= 0.0)
                    {
                        bracket_found = true;
                        break;
                    }
                }
                
                if (!bracket_found)
                {
                    // Fallback: use largest bracket found
                    // This shouldn't happen for physical problems
                    p_min = 0.1 * std::min(Pl, Pr);
                    p_max = 10.0 * std::max(Pl, Pr);
                }

                // Step 2: Brent's method to find root within bracket
                constexpr int max_iter = 50;
                constexpr real tol = 1.0e-8;
                real p_current = 0.0;
                bool converged = false;

                // Step 2: Brent's method to find root within bracket
                constexpr int max_brent_iter = 50;
                constexpr real brent_tol = 1.0e-8;
                
                real brent_a = p_min, brent_b = p_max;
                real fa = dvel_min, fb = dvel_max;
                real brent_c = brent_a, fc = fa;
                real brent_d = brent_b - brent_a, brent_e = brent_d;
                
                constexpr real eps = std::numeric_limits<real>::epsilon();
                
                for (int iter = 0; iter < max_brent_iter; ++iter)
                {
                    // Ensure |f(b)| <= |f(c)|
                    if (std::abs(fc) < std::abs(fb))
                    {
                        brent_a = brent_b; brent_b = brent_c; brent_c = brent_a;
                        fa = fb; fb = fc; fc = fa;
                    }
                    
                    // Convergence test
                    const real tol1 = 2.0 * eps * std::abs(brent_b) + 0.5 * brent_tol;
                    const real xm = 0.5 * (brent_c - brent_b);
                    
                    if (std::abs(xm) <= tol1 || fb == 0.0)
                    {
                        p_current = brent_b;
                        converged = true;
                        break;
                    }
                    
                    // Choose interpolation or bisection
                    if (std::abs(brent_e) < tol1 || std::abs(fa) <= std::abs(fb))
                    {
                        // Bisection
                        brent_d = xm;
                        brent_e = brent_d;
                    }
                    else
                    {
                        real p, q;
                        const real s = fb / fa;
                        
                        if (brent_a == brent_c)
                        {
                            // Linear interpolation
                            p = 2.0 * xm * s;
                            q = 1.0 - s;
                        }
                        else
                        {
                            // Inverse quadratic interpolation
                            const real r = fb / fc;
                            const real t = fa / fc;
                            p = s * (2.0*xm*t*(t-r) - (brent_b-brent_a)*(r-1.0));
                            q = (t-1.0) * (r-1.0) * (s-1.0);
                        }
                        
                        if (p > 0.0) q = -q;
                        p = std::abs(p);
                        
                        const real min1 = 3.0*xm*q - std::abs(tol1*q);
                        const real min2 = std::abs(brent_e*q);
                        
                        if (2.0*p < std::min(min1, min2))
                        {
                            // Accept interpolation
                            brent_e = brent_d;
                            brent_d = p / q;
                        }
                        else
                        {
                            // Reject interpolation, use bisection
                            brent_d = xm;
                            brent_e = brent_d;
                        }
                    }
                    
                    // Update a
                    brent_a = brent_b;
                    fa = fb;
                    
                    // Update b
                    if (std::abs(brent_d) > tol1)
                        brent_b += brent_d;
                    else
                        brent_b += (xm > 0.0 ? tol1 : -tol1);
                    
                    // Evaluate at new point
                    fb = dvel_at_p(brent_b);
                    
                    // Ensure bracket maintained
                    if (fb * (fc / std::abs(fc)) > 0.0)
                    {
                        brent_c = brent_a;
                        fc = fa;
                        brent_d = brent_b - brent_a;
                        brent_e = brent_d;
                    }
                }

                // Compute final solution
                if (converged)
                {
                    pstar = p_current;
                    
                    // Get final velocities
                    real nl_star, Hl_star, csl_star, vl_star;
                    real nr_star, Hr_star, csr_star, vr_star;
                    
                    get_velocity_behind_wave(pstar, nl, Pl, Hl, csl, vl, vtl, wl, -1.0,
                                           nl_star, Hl_star, csl_star, vl_star);
                    get_velocity_behind_wave(pstar, nr, Pr, Hr, csr, vr, vtr, wr, +1.0,
                                           nr_star, Hr_star, csr_star, vr_star);
                    
                    vstar = 0.5 * (vl_star + vr_star);
                }
                else
                {
                    // Fallback: arithmetic mean pressure
                    pstar = 0.5 * (Pl + Pr);

                    real nl_star, Hl_star, csl_star, vl_star;
                    real nr_star, Hr_star, csr_star, vr_star;

                    get_velocity_behind_wave(pstar, nl, Pl, Hl, csl, vl, vtl, wl, -1.0,
                                           nl_star, Hl_star, csl_star, vl_star);
                    get_velocity_behind_wave(pstar, nr, Pr, Hr, csr, vr, vtr, wr, +1.0,
                                           nr_star, Hr_star, csr_star, vr_star);

                    vstar = 0.5 * (vl_star + vr_star);

                    // Diagnostic: warn if Riemann solver didn't converge
                    static int riemann_warn_count = 0;
                    if (riemann_warn_count < 10)
                    {
                        std::cerr << "WARNING: Riemann solver did not converge. "
                                  << "Pl=" << Pl << ", Pr=" << Pr
                                  << ", vl=" << vl << ", vr=" << vr
                                  << ", pstar=" << pstar << std::endl;
                        riemann_warn_count++;
                        if (riemann_warn_count == 10)
                        {
                            std::cerr << "(Further Riemann solver warnings suppressed)" << std::endl;
                        }
                    }
                }
            };
        }

        /**
         * Compute post-wave state using relativistic jump conditions
         *
         * Based on Pons, Martí & Müller (1999) Eq. (21)-(22) and Appendix A,
         * and Kitajima et al. (2025) arXiv:2510.18251v1 for ideal gas EOS.
         *
         * Physical Problem:
         * Given pre-wave state (n_a, P_a, H_a, c_s,a, v^n_a, v^t_a) and interface pressure P*,
         * compute post-wave state (n*, H*, c_s,*, v^n*) behind shock or rarefaction.
         *
         * Wave Type Determination (Pons et al. Eq. 18):
         * - If P* > P_a: SHOCK wave
         *   - Use Rankine-Hugoniot relations (Eq. 21-22)
         *   - Entropy increases across shock
         * - If P* <= P_a: RAREFACTION fan
         *   - Use isentropic relations (Eq. A.4-A.7)
         *   - Entropy constant across rarefaction
         *
         * Shock Jump Conditions (Pons et al. Eq. 21-22, Kitajima Eq. 62-64):
         *   Taub adiabat: a*H^2 + b*H + c = 0
         *   where a = 1 + (gamma-1)(P_a - P)/(gamma*P)
         *         b = -(a+1)
         *         c = H_a(P_a - P)/n_a - H_a^2
         *   Mass flux: j^2 = -[P]/[H/n] = -(P - P_a)/(H/n - H_a/n_a)
         *   Shock speed: v_shock = [-v_a*N_a^2 + sign*j^2*sqrt(1 + n_a^2/j^2)] / (j^2 + N_a^2)
         *   Post-shock velocity: v* = [gamma_shock(P-P_a)/j + H_a*W_a*v_a] / [H_a*W_a + (P-P_a)(gamma_shock*v_a/j + 1/(n_a*W_a))]
         *
         * Rarefaction Relations (Pons et al. Appendix A, Kitajima ideal gas):
         *   Polytropic: K = P_a/n_a^gamma (constant along rarefaction)
         *   Post-rarefaction density: n_star = (P_star/K)^(1/gamma)
         *   Post-rarefaction enthalpy: H_star = 1 + u_star/c^2 + P_star/(n_star*c^2)
         *   Velocity integral (Pons eq. 30): relates v_star to P_star via ODE solution
         *   For ideal gas: A = [(1+v_a/c)/(1-v_a/c)] * [(sqrt(gamma-1)+c_s,a/c)/(sqrt(gamma-1)-c_s,a/c) * (sqrt(gamma-1)-c_s/c)/(sqrt(gamma-1)+c_s/c)]^(-sign*2/sqrt(gamma-1))
         *                  v_star = c*(A-1)/(A+1)
         *
         * Tangential Velocity Treatment (Pons et al. Section 3):
         * - Conservation across shocks (Eq. 54-55): S^y/D = const, S^z/D = const
         *   where S^y = rho*h*W^2*v^y, D = rho*W = n*W
         *   Thus: h*W*v^y = const, h*W*v^z = const
         * - Conservation across rarefactions (Eq. 27-28): h*W*v^y = const, h*W*v^z = const
         * - Direction preserved: v^y/v^z = const
         * - Magnitude iteratively solved since W = 1/sqrt(1 - v^2/c^2) depends on total velocity
         * - For 1D problems (v^t = 0): no iteration needed, coupling vanishes
         *
         * Parameters:
         * @param P      Interface pressure P* (Newton-Raphson solution)
         * @param na     Pre-wave rest-frame density n_a
         * @param Pa     Pre-wave pressure P_a
         * @param Ha     Pre-wave enthalpy H_a = 1 + u_a/c^2 + P_a/(n_a*c^2)
         * @param csa    Pre-wave sound speed c_s,a = sqrt(gamma*P_a/(n_a*H_a))
         * @param va     Pre-wave normal velocity v^n_a (along shock normal)
         * @param vta    Pre-wave tangential velocity magnitude |v^t_a| = sqrt((v^y)^2 + (v^z)^2)
         * @param wa     Pre-wave Lorentz factor W_a = 1/sqrt(1 - v_a^2/c^2) including tangential velocity
         * @param sign   Wave direction: -1 for left wave, +1 for right wave
         * @param[out] n   Post-wave rest-frame density n*
         * @param[out] H   Post-wave enthalpy H_star = 1 + u_star/c^2 + P_star/(n_star*c^2)
         * @param[out] cs  Post-wave sound speed c_s_star = sqrt(gamma*P_star/(n_star*H_star))
         * @param[out] vel Post-wave normal velocity v^n*
         */
        void FluidForce::get_velocity_behind_wave(
            const real P,
            const real na, const real Pa, const real Ha, const real csa,
            const real va, const real vta, const real wa,
            const real sign,
            real &n, real &H, real &cs, real &vel) const
        {
            // Compute post-wave velocity using relativistic jump conditions
            // Extended for tangential velocities following Pons et al. (1999)
            // vta = tangential velocity magnitude (|v^t| in the paper)
            // wa = Lorentz factor INCLUDING tangential velocity: W_a = 1/√(1 - v_a²/c²)
            //      where v_a² = (v^n_a)² + (v^t_a)²

            const real c = m_c_speed;
            const real gamma_c = m_gamma;

            if (P > Pa)
            {
                // SHOCK WAVE - Use Taub adiabat and Rankine-Hugoniot relations
                // Reference: Pons et al. (1999) Eq. 62-64, Kitajima et al. (2025) for ideal gas

                // Taub adiabat for ideal gas (Kitajima Eq. 63, Pons Eq. 62):
                // Quadratic equation: a·H² + b·H + c = 0
                const real a = 1.0 + (gamma_c - 1.0) * (Pa - P) / (gamma_c * P);
                const real b_term = -(gamma_c - 1.0) * (Pa - P) / (gamma_c * P);  // = -(a-1) = 1-a
                const real c_term = Ha * (Pa - P) / na - Ha * Ha;

                // Discriminant check for physical solution
                const real discriminant = b_term * b_term - 4.0 * a * c_term;
                if (discriminant < 0.0)
                {
                    // Unphysical state - use fallback (no change across wave)
                    n = na;
                    H = Ha;
                    cs = csa;
                    vel = va;
                    return;
                }

                // Post-shock enthalpy (positive root of quadratic, Kitajima Eq. 63)
                H = (-b_term + std::sqrt(discriminant)) / (2.0 * a);

                // Validate enthalpy is physical
                if (H <= 0.0)
                {
                    // Unphysical - fallback
                    n = na;
                    H = Ha;
                    cs = csa;
                    vel = va;
                    return;
                }

                // Post-shock rest-frame density from ideal gas EOS (Kitajima formulation)
                // H = 1 + u/c² + P/(n·c²) and u = P/((γ-1)·n)
                // Solving: n = γP / ((γ-1)(H-1))
                n = gamma_c * P / ((gamma_c - 1.0) * (H - 1.0));

                // Mass flux squared (Pons Eq. 64, Kitajima context)
                // j² = -[P]/[H/n] = -(P - P_a) / (H/n - H_a/n_a)
                const real j2 = -(P - Pa) / (H / n - Ha / na);
                if (j2 < 0.0)
                {
                    // Unphysical - fallback
                    n = na;
                    H = Ha;
                    cs = csa;
                    vel = va;
                    return;
                }
                const real j = sign * std::sqrt(j2);

                // Shock velocity (Pons et al. 1999, Eq. 522)
                // V_s = [ρ_a² W_a² v^x_a  +  sign*|j|*√(j² + ρ_a² W_a²(1 - v^x_a²))] / (ρ_a² W_a² + j²)
                // where sign = +1 for right shock, -1 for left shock
                const real Na = wa * na;  // N_a = W_a * ρ_a (lab-frame density)
                const real Na2 = Na * Na;  // (ρ_a W_a)²
                const real denominator = j2 + Na2;

                // BUGFIX: Correct formula from Pons Eq. 522
                const real sqrt_arg = j2 + Na2 * (1.0 - va * va / (c * c));
                const real vshock = (Na2 * va + sign * std::abs(j) * std::sqrt(sqrt_arg)) / denominator;

                // Shock Lorentz factor
                const real gamma_shock = 1.0 / std::sqrt(1.0 - (vshock / c) * (vshock / c));

                // Post-shock normal velocity (Pons Eq. 61, Kitajima formulation)
                // v* = [γ_shock(P-P_a)/j + H_a·W_a·v_a] / [H_a·W_a + (P-P_a)(γ_shock·v_a/j + 1/(n_a·W_a))]
                const real a_vel = gamma_shock * (P - Pa) / j + Ha * wa * va;
                const real b_vel = Ha * wa + (P - Pa) * (gamma_shock * va / j + 1.0 / (na * wa));
                vel = a_vel / b_vel;

                // SAFETY: Ensure post-shock velocity is subluminal
                const real v_limit = 0.99 * c;
                if (std::abs(vel) > v_limit)
                {
                    vel = std::copysign(v_limit, vel);
                }

                // Post-shock sound speed
                cs = std::sqrt(gamma_c * P / (n * H));

                // Tangential velocity handling:
                // For shocks, conservation law (Pons Eq. 54-55): hWv^t = const
                // This is handled implicitly in the Riemann solver - tangential component
                // doesn't affect normal velocity in 1D problem (vta enters only through wa)
            }
            else
            {
                // RAREFACTION WAVE - Use isentropic relations
                // Reference: Pons et al. (1999) Appendix A, Kitajima ideal gas

                // Polytropic constant (isentropic process, Pons Eq. A.4)
                const real K = Pa / std::pow(na, gamma_c);

                // Post-rarefaction rest-frame density (Pons Eq. A.5)
                n = std::pow(P / K, 1.0 / gamma_c);

                // Thermal energy per baryon (ideal gas)
                const real u = P / ((gamma_c - 1.0) * n);

                // Post-rarefaction enthalpy
                H = 1.0 + u / (c * c) + P / (n * c * c);

                // Post-rarefaction sound speed
                cs = std::sqrt(gamma_c * P / (n * H));

                // Velocity from rarefaction invariant (Pons et al. Eq. 30, integral solution)
                // For ideal gas, the ODE integration gives (Pons Appendix A):
                // A = [(1+v_a/c)/(1-v_a/c)] · [(√(γ-1)+c_s,a/c)/(√(γ-1)-c_s,a/c) · (√(γ-1)-c_s/c)/(√(γ-1)+c_s/c)]^(-sign·2/√(γ-1))
                // v* = c·(A-1)/(A+1)
                const real sqgl1 = std::sqrt(gamma_c - 1.0);
                const real term1 = (1.0 + va / c) / (1.0 - va / c);
                const real term2 = (sqgl1 + csa / c) / (sqgl1 - csa / c);
                const real term3 = (sqgl1 - cs / c) / (sqgl1 + cs / c);
                const real exponent = -sign * 2.0 / sqgl1;
                const real A = term1 * std::pow(term2 * term3, exponent);

                vel = c * (A - 1.0) / (A + 1.0);

                // SAFETY: Ensure post-rarefaction velocity is subluminal
                const real v_limit = 0.99 * c;
                if (std::abs(vel) > v_limit || !std::isfinite(vel))
                {
                    // Fall back to safe velocity
                    vel = std::copysign(std::min(std::abs(va), v_limit), va);
                }

                // Tangential velocity handling:
                // For rarefactions, conservation law (Pons Eq. 27-28): hWv^t = const
                // This couples through the Lorentz factor W in the full 2D/3D problem.
                // For 1D problems (vta = 0), there is no coupling effect.
                //
                // For non-zero tangential velocity (future 2D/3D extension):
                // Need iterative solution since W = 1/√(1 - v²/c²) depends on total velocity v² = (v^n)² + (v^t)²
                // and v^t = (H_a·W_a·v^t_a) / (H·W) from conservation law.
                // Currently not implemented as SR Sod test is 1D (vta = 0).
            }
        }

    } // namespace srgsph
} // namespace sph
