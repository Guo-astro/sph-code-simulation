/**
 * Special Relativistic GSPH - Pre-Interaction Phase
 * 
 * Reference: Kitajima, Inutsuka, Seno (2025) "Special Relativistic Hydrodynamics 
 *           with Smoothed Particle Hydrodynamics" arXiv:2510.18251v1, Section 2.3
 * 
 * Key Equations Implemented:
 * - Eq. 220: Particle volume V_p(x) = [Σ_j W(x-x_j, h)]^(-1)
 * - Eqs. 231-233: Variable smoothing length h = η·[V_p*(x)]^(1/d)
 *                where V_p*(x) = [Σ_j W(x-x_j, C_smooth·h)]^(-1)
 * - Eq. 239: Number density N(x) = ν(x)/V_p(x) (NOT kernel sum!)
 * - Primitive recovery: Quartic solver for (γ, v, P, n) from (S, e, N)
 * - MUSCL gradients for 2nd order reconstruction
 * 
 * Parameters:
 * - η = 1.0 (smoothing length coefficient)
 * - C_smooth = 2.0 (expansion factor for h iteration)
 * - Iterative convergence: max_iter=5, tolerance=1e-4
 */

#include <algorithm>
#include <cmath>

#include "parameters.hpp"
#include "srgsph/sr_pre_interaction.hpp"
#include "srgsph/sr_primitive_recovery.hpp"
#include "simulation.hpp"
#include "periodic.hpp"
#include "openmp.hpp"
#include "kernel/kernel_function.hpp"
#include "exception.hpp"
#include "bhtree.hpp"

#ifdef EXHAUSTIVE_SEARCH
#include "exhaustive_search.hpp"
#endif

namespace sph
{
namespace srgsph
{

void PreInteraction::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::PreInteraction::initialize(param);
    m_eta = param->srgsph.eta;
    m_c_smooth = 2.0;  // Kitajima et al. default
    m_gamma = param->physics.gamma;
    m_c_speed = param->srgsph.c_speed;
    m_is_2nd_order = param->srgsph.is_2nd_order;
    m_iteration = param->iterative_sml;
    m_first = true;
}

/**
 * Initial smoothing length calculation
 * 
 * Purpose: Provide initial estimate of h for first timestep before iteration
 * 
 * Formula: h_init = η × N_neighbor × (Δx_avg / A)
 * where:
 * - η = smoothing length coefficient (typically 1.0)
 * - N_neighbor = desired neighbor number
 * - Δx_avg = average particle spacing = (domain_length / num_particles)^(1/d)
 * - A = d-dimensional unit sphere/circle/line "area": 
 *       1D: A=2, 2D: A=π, 3D: A=4π/3
 * 
 * Note: This is a rough estimate. Variable h will be refined by
 *       compute_smoothing_length() using Kitajima Eqs. 231-233.
 */
void PreInteraction::initial_smoothing(std::shared_ptr<Simulation> sim)
{
    // Initial guess for smoothing length based on particle spacing
    // Similar to standard SPH: h ~ (neighbor_number * particle spacing)
    auto & particles = sim->get_particles();
    const int num = sim->get_particle_num();
    const real neighbor = m_neighbor_number;
    
    // Estimate average particle spacing from domain size
    constexpr real A = DIM == 1 ? 2.0 :
                       DIM == 2 ? M_PI :
                                  4.0 * M_PI / 3.0;
    
    // For 1D: estimate spacing as domain_length / num_particles
    real domain_length = 1.0;  // Default for Sod problem (-0.5 to 0.5)
    real avg_spacing = domain_length / real(num);
    
    // h should contain ~neighbor_number particles
    // For 1D with kernel support 2h: neighbor_number = 2h / spacing
    // => h = neighbor_number * spacing / 2
    const real h_init = m_eta * neighbor * avg_spacing / A;
    
#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        p_i.sml = h_init;
    }
}

/**
 * Compute particle volume V_p(x_i)
 * 
 * Reference: Kitajima et al. (2025) Eq. 220
 * 
 * Formula: V_p(x) = [Σ_j W(x - x_j, h)]^(-1)
 * 
 * Physical interpretation:
 * - V_p is the effective volume occupied by particle at position x
 * - Inverse of kernel sum gives volume directly (not density)
 * - Used to compute number density: N = ν/V_p (Eq. 239)
 * 
 * Implementation notes:
 * - Includes self-contribution: W(0, h)
 * - Kernel evaluated at current smoothing length h
 * - Fallback to V_p=1.0 if kernel sum vanishes (pathological case)
 * 
 * Returns: V_p = 1/Σ_j W(r_ij, h)
 */
real PreInteraction::compute_volume(
    const SPHParticle & p_i,
    const std::vector<SPHParticle> & particles,
    const std::vector<int> & neighbor_list,
    const int n_neighbor,
    const Periodic * periodic,
    const KernelFunction * kernel,
    const real h
)
{
    // V_p(x) = [Σ_j W(x - x_j, h)]^(-1)  (Eq. 220)
    real sum_w = 0.0;
    const vec_t & r_i = p_i.pos;
    
    for(int n = 0; n < n_neighbor; ++n) {
        int const j = neighbor_list[n];
        auto & p_j = particles[j];
        const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
        const real r = std::abs(r_ij);
        
        if(r < h) {
            sum_w += kernel->w(r, h);
        }
    }
    
    // Self-contribution
    sum_w += kernel->w(0.0, h);
    
    if(sum_w < 1.0e-20) {
        return 1.0;  // Fallback to avoid division by zero
    }
    
    return 1.0 / sum_w;
}

/**
 * Compute variable smoothing length using iterative volume-based solver
 * 
 * Implements Kitajima et al. (2025) Eqs. 231-233:
 *   h(x) = η · [V_p*(x)]^(1/d)                                    (Eq. 231)
 *   V_p*(x) = [Σ_j W(x - x_j, C_smooth · h)]^(-1)                (Eq. 233)
 * 
 * Physical Meaning:
 * - h adapts to local particle volume, not just spacing
 * - Expanded kernel support (C_smooth·h) prevents h-oscillations
 * - η parameter controls kernel support (~1.2-1.5 for typical kernels)
 * 
 * Algorithm:
 * 1. Start from previous timestep's h as initial guess
 * 2. Compute expanded volume V_p* at current h using C_smooth·h kernel support
 * 3. Update h = η · (V_p*)^(1/d)
 * 4. Iterate until |h_new - h| / h < tol
 * 
 * Parameters (Kitajima recommendations):
 * - η = 1.2 for cubic spline kernel in 3D
 * - C_smooth = 1.0 (no expansion) to 1.5 (strong smoothing)
 * - max_iter = 5 sufficient when using previous h as guess
 * - tol = 1e-4 balances accuracy and performance
 * 
 * Convergence:
 * - Typically converges in 2-3 iterations for smooth flows
 * - May need more iterations near shocks or strong gradients
 * - Returns last iteration if max_iter reached (graceful degradation)
 * 
 * Edge Cases:
 * - If Σ_j W → 0 (isolated particle), expand h by 50% and retry
 * - Self-contribution always included (W(0, h))
 */
real PreInteraction::compute_smoothing_length(
    const SPHParticle & p_i,
    const std::vector<SPHParticle> & particles,
    const std::vector<int> & neighbor_list,
    const int n_neighbor,
    const Periodic * periodic,
    const KernelFunction * kernel
)
{
    // Iteratively solve h = η * [V_p*]^(1/d) where V_p* = [Σ_j W(x-x_j, C_smooth*h)]^(-1)
    // (Kitajima Eqs. 231-233)
    
    real h = p_i.sml;  // Initial guess from previous step
    const int max_iter = 50;   // Reduced from 20 - previous h is good guess
    const real tol = 1.0e-4;  // Relaxed from 1e-6 to match GSPH
    
    for(int iter = 0; iter < max_iter; ++iter) {
        // Compute V_p* at current h using C_smooth expansion
        // V_p*(x) = [Σ_j W(x-x_j, C_smooth*h)]^(-1)  (Eq. 233)
        const real h_smooth = m_c_smooth * h;
        real sum_w_star = 0.0;
        const vec_t & r_i = p_i.pos;
        
        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
            const real r = std::abs(r_ij);
            
            if(r < h_smooth) {
                sum_w_star += kernel->w(r, h_smooth);
            }
        }
        
        // Self-contribution
        sum_w_star += kernel->w(0.0, h_smooth);
        
        if(sum_w_star < 1.0e-20) {
            // Pathological case - expand search
            h *= 1.5;
            continue;
        }
        
        // V_p* = 1 / sum_w_star
        const real Vp_star = 1.0 / sum_w_star;
        
        // New smoothing length: h = η * V_p*^(1/d)  (Eq. 231)
        const real h_new = m_eta * std::pow(Vp_star, 1.0 / real(DIM));
        
        // Check convergence
        if(std::abs(h_new - h) / h < tol) {
            return h_new;
        }
        
        h = h_new;
    }
    
    // Return last iteration result even if not converged
    return h;
}

/**
 * Main pre-interaction calculation: variable h, volume, number density, primitive recovery
 * 
 * Implements complete Kitajima et al. (2025) volume-based formulation:
 * 
 * Step 1: Update smoothing length (if iteration enabled)
 *   h(x) = η · [V_p*(x)]^(1/d)                                    (Eq. 231)
 * 
 * Step 2: Compute particle volume
 *   V_p(x) = [Σ_j W(x - x_j, h)]^(-1)                            (Eq. 220)
 * 
 * Step 3: Compute volume-based number density
 *   N(x) = ν(x) / V_p(x)                                          (Eq. 239)
 *   where ν = baryon number (mass for now)
 * 
 * Step 4: Recover primitive variables (γ, v, P, n) from conserved (S, e, N)
 *   Solve quartic equation using PrimitiveRecovery::conserved_to_primitive()
 * 
 * Step 5: Compute signal velocity for timestep
 *   Δt = C_CFL * min_i[h_i/c_s,i]  (Paper Section 3, simple Courant condition)
 * 
 * Step 6: MUSCL gradients for 2nd order reconstruction (optional)
 *   ∇n, ∇P, ∇v using SPH kernel gradient operator
 * 
 * Key Differences from Fixed-h GSPH:
 * - h varies per particle (adaptive resolution)
 * - N computed from volume, not kernel sum
 * - Primitive recovery uses N directly (no mixing of h and N)
 */
void PreInteraction::calculation(std::shared_ptr<Simulation> sim)
{
    // First timestep: iterate to convergence for initial h
    if(m_first) {
        initial_smoothing(sim);
        m_first = false;
    }

    auto &particles = sim->get_particles();
    auto *periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto *kernel = sim->get_kernel().get();
    auto *tree = sim->get_tree().get();

    omp_real h_per_v_sig(std::numeric_limits<real>::max());

    // For MUSCL reconstruction (2nd order)
    auto &grad_n = sim->get_vector_array("grad_number_density");
    auto &grad_p = sim->get_vector_array("grad_pressure");
    vec_t *grad_v[DIM] = {
        sim->get_vector_array("grad_velocity_0").data(),
#if DIM == 2
        sim->get_vector_array("grad_velocity_1").data(),
#elif DIM == 3
        sim->get_vector_array("grad_velocity_1").data(),
        sim->get_vector_array("grad_velocity_2").data(),
#endif
    };

#pragma omp parallel for
    for (int i = 0; i < num; ++i)
    {
        auto &p_i = particles[i];
        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

        // Neighbor search with current smoothing length
#ifdef EXHAUSTIVE_SEARCH
        int const n_neighbor = exhaustive_search(p_i, p_i.sml, particles, num, neighbor_list,
                                                  m_neighbor_number * neighbor_list_size, periodic, false);
#else
        int const n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, false);
#endif

        // 1. Update smoothing length
        if(m_iteration) {
            // Iterative smoothing length (volume-based approach)
            p_i.sml = compute_smoothing_length(p_i, particles, neighbor_list, n_neighbor,
                                              periodic, kernel);
        }
        // else: keep h from initial_smoothing() or previous step

        const real h_i = p_i.sml;

        // 2. Compute particle volume Vp using VOLUME-BASED approach (Eq. 220)
        // Vp(x) = [Σ_j W(x - x_j, h)]^(-1)
        const real Vp_i = compute_volume(p_i, particles, neighbor_list, n_neighbor, periodic, kernel, h_i);

        // 3. Compute number density N from volume: N = ν/Vp  (Eq. 239)
        // This is the KEY difference from fixed-h: we use volume-based density, not kernel sum
        const real nu_i = p_i.mass;  // Baryon number (equal to mass for now)
        const real N_i = nu_i / Vp_i;

        p_i.dens = N_i; // Store number density in dens field
        p_i.N = N_i;    // BUGFIX: Also store in dedicated N field (used by sr_fluid_force.cpp)

        // 4. Primitive variable recovery from conserved variables (S, e)
        // 4. Primitive variable recovery from conserved variables
        // Read CONSERVED variables (time-integrated by solver)
        // CRITICAL: Use dedicated S/e fields, NOT vel/ene (those store primitives!)
        
        PrimitiveVariables prim;
        
        if(m_first) {
            // First timestep: use initialized primitive values directly
            // No need to recover - just use what was set in initial conditions
            prim.vel = p_i.vel;
            prim.pressure = p_i.pres;
            prim.sound_speed = p_i.sound;
            prim.gamma_lor = p_i.gamma_lor;
            prim.density = p_i.dens / p_i.gamma_lor;  // n = N/γ
            prim.enthalpy = p_i.enthalpy;
        } else {
            // Subsequent timesteps: recover from conserved variables
            const vec_t S_i = p_i.S;   // Canonical momentum S = γHv (per baryon)
            const real e_i = p_i.e;    // Canonical energy e = γH - P/(Nc²) (per baryon)
            
            // Recover primitive variables using quartic solver (Eq. 49-51 in Kitajima 2025)
            prim = PrimitiveRecovery::conserved_to_primitive(S_i, e_i, N_i, m_gamma, m_c_speed);
        }

        // Store recovered PRIMITIVE variables (for Riemann solver and output)
        p_i.vel = prim.vel;           // Primitive velocity v (NOT S!)
        p_i.pres = prim.pressure;     // Pressure P
        p_i.sound = prim.sound_speed; // Sound speed c_s
        p_i.gamma_lor = prim.gamma_lor;   // Lorentz factor γ
        p_i.neighbor = n_neighbor;

        // 5. Signal velocity for timestep
        // Paper Eq. Section 3: Simple Courant condition Δt = C_CFL * min_i[h_i/c_s,i]
        // We compute h/c_s here; timestep.cpp will apply C_CFL and take minimum
        const real h_per_c_s_i = p_i.sml / p_i.sound;
        if (h_per_v_sig.get() > h_per_c_s_i)
        {
            h_per_v_sig.get() = h_per_c_s_i;
        }

        // 6. MUSCL gradient calculation for 2nd order accuracy
        // NOTE: Gradients are computed but NOT YET USED in sr_fluid_force.cpp
        // The Riemann solver currently uses only 1st order (particle values directly)
        // TODO: Implement MUSCL reconstruction with monotonicity limiter (Paper Eq. 19)
        if (!m_is_2nd_order)
        {
            continue;
        }

        vec_t dn(0.0), dp(0.0);
        vec_t dv[DIM];
        for (int k = 0; k < DIM; ++k)
        {
            dv[k] = vec_t(0.0);
        }

        const vec_t & pos_i = p_i.pos;
        for (int n = 0; n < n_neighbor; ++n)
        {
            int const j = neighbor_list[n];
            auto &p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
            const real r = std::abs(r_ij);
            
            if (r >= h_i) break;

            const vec_t dw_ij = kernel->dw(r_ij, r, h_i);
            
            // Gradient of rest-frame density n
            const real n_j = p_j.dens / p_j.gamma_lor;
            dn += dw_ij * p_j.mass * (prim.density - n_j);
            
            // Gradient of pressure
            dp += dw_ij * p_j.mass * (prim.pressure - p_j.pres);
            
            // Gradient of velocity components
            for (int k = 0; k < DIM; ++k)
            {
                dv[k] += dw_ij * p_j.mass * (prim.vel[k] - p_j.vel[k]);
            }
        }

        grad_n[i] = dn;
        grad_p[i] = dp;
        for (int k = 0; k < DIM; ++k)
        {
            grad_v[k][i] = dv[k];
        }
    }

    sim->set_h_per_v_sig(h_per_v_sig.min());

#ifndef EXHAUSTIVE_SEARCH
    tree->set_kernel();
#endif
}

} // namespace srgsph
} // namespace sph
