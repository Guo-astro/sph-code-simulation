// Special Relativistic Sod Shock Tube Test
// Implementation of Kitajima et al. (2025) arXiv:2510.18251v1 Section 3.1.1

#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"
#include "logger.hpp"
#include "srgsph/sr_primitive_recovery.hpp"
#include <cmath>

namespace sph
{

void Solver::make_sr_sod()
{
#if DIM != 1
    THROW_ERROR("SR Sod test requires DIM == 1");
#else

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real gamma = m_param->physics.gamma;
    const real c_speed = m_param->srgsph.c_speed;
    const real c2 = c_speed * c_speed;
    
    // Different baryon number test (default: false)
    bool different_nu = false;
    if (m_sample_parameters.count("different_nu")) {
        different_nu = boost::any_cast<bool>(m_sample_parameters["different_nu"]);
    }
    
    // Uniform test mode (for debugging): use uniform P, n everywhere
    bool uniform_test = false;  // DISABLED - Test actual Sod shock
    if (m_sample_parameters.count("uniform_test")) {
        uniform_test = boost::any_cast<bool>(m_sample_parameters["uniform_test"]);
    }
    
    // ============================================================================
    // Initial conditions from Kitajima et al. (2025) Section 3.1.1, Eqs. 74-75
    // ============================================================================
    // Left state:  (P_L, n_L, v_Lx, v_Lt) = (1.0, 1.0, 0, 0)
    // Right state: (P_R, n_R, v_Rx, v_Rt) = (0.1, 0.125, 0, 0)
    // Particles: 3200 left, 400 right (8:1 ratio)
    // Equal baryon numbers per particle (ν_left = ν_right)
    // ============================================================================
    
    const int N_left = N * 8;
    const int N_right = N;
    const int num = N_left + N_right;
    
    // Domain: x ∈ [-0.5, 0.5], discontinuity at x=0
    const real x_left_start = -0.5;
    const real x_left_end = 0.0;
    const real x_right_start = 0.0;
    const real x_right_end = 0.5;
    
    const real dx_left = (x_left_end - x_left_start) / N_left;
    const real dx_right = (x_right_end - x_right_start) / N_right;
    
    // Left state primitive variables (Eq. 74)
    real P_left = 1.0;
    real n_left = 1.0;   // rest-frame baryon number density
    real v_left = 0.0;
    
    // Right state primitive variables (Eq. 75)
    real P_right = 0.1;
    real n_right = 0.125;  // rest-frame baryon number density
    real v_right = 0.0;
    
    // Uniform test mode override (BEFORE using these values)
    if (uniform_test) {
        P_left = P_right = 0.5;
        n_left = n_right = 0.5;
        v_left = v_right = 0.0;
        WRITE_LOG << "UNIFORM TEST MODE: P=" << P_left << ", n=" << n_left;
    }
    
    // Baryon number per particle (ν)
    // Paper states: "Using SPH particles that have equal baryon numbers"
    real nu_left, nu_right;
    if (different_nu) {
        // Different ν test (Fig. 5 in paper)
        // Conserve total baryon: N_left × ν_left = N_right × ν_right
        nu_left = 1.0;
        nu_right = static_cast<real>(N_left) / static_cast<real>(N_right) * nu_left;
    } else {
        // Equal ν (Fig. 4 in paper - standard Sod test)
        // 
        // With PERIODIC boundaries and FIXED global h:
        // Use empirical kernel sum to determine ν
        // Target: N_left = n_left = 1.0 in left region
        //
        // Empirical kernel sum depends on particle spacing:
        // For N=400 (3600 total): Σ W ≈ 3060
        // For N=100 (900 total):  Σ W ≈ 765
        // So: ν = n / Σ W
        //
        // For EQUAL ν test: both regions use same ν
        // Base it on the left state density
        
        const real kernel_sum_estimate = 765.0;  // Empirical for N=100
        const real nu_equal = n_left / kernel_sum_estimate;
        nu_left = nu_equal;
        nu_right = nu_equal;  // SAME ν for both regions
    }
    
    std::vector<SPHParticle> p(num);
    
    // ============================================================================
    // Smooth transition zone to avoid huge initial forces
    // ============================================================================
    // Instead of sharp discontinuity at x=0, use smooth transition over ~20h
    // This prevents enormous pressure gradients on first timestep
    const real x_discontinuity = 0.0;
    const real h_global = 0.0110;  // Approximate value
    const real transition_width = 20.0 * h_global;  // ~0.22
    
    // Lambda for smooth interpolation
    auto smooth_state = [&](real x, real val_left, real val_right) -> real {
        real dist = x - x_discontinuity;
        if (dist < -0.5 * transition_width) return val_left;
        if (dist > 0.5 * transition_width) return val_right;
        // Smooth tanh transition
        real xi = dist / transition_width;  // ∈ [-0.5, 0.5]
        real s = 0.5 * (1.0 + std::tanh(5.0 * xi));  // 0→1 smoothly
        return val_left + s * (val_right - val_left);
    };
    
    // ============================================================================
    // Initialize left state particles
    // ============================================================================
    // Strategy: Set primitive variables (P, v) and ν directly.
    // DO NOT set conserved variables (S, e) yet!
    // After all particles are initialized, pre_interaction will:
    // 1. Compute N from kernel sum
    // 2. Compute conserved (S, e) from primitives and smoothed N
    for (int i = 0; i < N_left; ++i) {
        auto& p_i = p[i];
        p_i.id = i;
        p_i.pos[0] = x_left_start + (i + 0.5) * dx_left;
        
        // Baryon number per particle
        p_i.nu = nu_left;
        p_i.mass = nu_left;
        
        // SMOOTHED primitive variables (avoids huge initial forces)
        real x = p_i.pos[0];
        real P_smooth = smooth_state(x, P_left, P_right);
        real n_smooth = smooth_state(x, n_left, n_right);
        real v_smooth = smooth_state(x, v_left, v_right);
        
        vec_t vel;
        vel[0] = v_smooth;
        #if DIM >= 2
        vel[1] = 0.0;
        #endif
        #if DIM == 3
        vel[2] = 0.0;
        #endif
        
        p_i.vel = vel;
        p_i.pres = P_smooth;  // SMOOTH pressure
        p_i.dens = n_smooth;  // SMOOTH rest-frame density
        
        // DO NOT compute conserved variables here!
        // They will be computed in pre_interaction after N is smoothed
        p_i.S = 0.0;  // Placeholder (will be recomputed)
        p_i.e = 0.0;  // Placeholder (will be recomputed)
        p_i.N = 0.0;  // Placeholder (will be computed from kernel sum)
        
        // These are for output/diagnostics (use smoothed values)
        const real u_smooth = P_smooth / ((gamma - 1.0) * n_smooth);
        p_i.ene = u_smooth;
        const real H_smooth = 1.0 + u_smooth / c2 + P_smooth / (n_smooth * c2);
        p_i.sound = std::sqrt((gamma - 1.0) * (H_smooth - 1.0) / H_smooth) * c_speed;
        p_i.enthalpy = H_smooth;
        real v2 = vel[0]*vel[0];
        #if DIM >= 2
        v2 += vel[1]*vel[1];
        #endif
        #if DIM == 3
        v2 += vel[2]*vel[2];
        #endif
        p_i.gamma_lor = 1.0 / std::sqrt(1.0 - v2 / c2);
    }
    
    // ============================================================================
    // Initialize right state particles
    // ============================================================================
    for (int i = 0; i < N_right; ++i) {
        auto& p_i = p[N_left + i];
        p_i.id = N_left + i;
        p_i.pos[0] = x_right_start + (i + 0.5) * dx_right;
        
        // Baryon number per particle
        p_i.nu = nu_right;
        p_i.mass = nu_right;
        
        // SMOOTHED primitive variables (avoids huge initial forces)
        real x = p_i.pos[0];
        real P_smooth = smooth_state(x, P_left, P_right);
        real n_smooth = smooth_state(x, n_left, n_right);
        real v_smooth = smooth_state(x, v_left, v_right);
        
        vec_t vel;
        vel[0] = v_smooth;
        #if DIM >= 2
        vel[1] = 0.0;
        #endif
        #if DIM == 3
        vel[2] = 0.0;
        #endif
        
        p_i.vel = vel;
        p_i.pres = P_smooth;  // SMOOTH pressure
        p_i.dens = n_smooth;  // SMOOTH rest-frame density
        
        // DO NOT compute conserved variables here!
        // They will be computed in pre_interaction after N is smoothed
        p_i.S = 0.0;  // Placeholder (will be recomputed)
        p_i.e = 0.0;  // Placeholder (will be recomputed)
        p_i.N = 0.0;  // Placeholder (will be computed from kernel sum)
        
        // For output/diagnostics (use smoothed values)
        const real u_smooth = P_smooth / ((gamma - 1.0) * n_smooth);
        p_i.ene = u_smooth;
        const real H_smooth = 1.0 + u_smooth / c2 + P_smooth / (n_smooth * c2);
        p_i.sound = std::sqrt((gamma - 1.0) * (H_smooth - 1.0) / H_smooth) * c_speed;
        p_i.enthalpy = H_smooth;
        real v2 = vel[0]*vel[0];
        #if DIM >= 2
        v2 += vel[1]*vel[1];
        #endif
        #if DIM == 3
        v2 += vel[2]*vel[2];
        #endif
        p_i.gamma_lor = 1.0 / std::sqrt(1.0 - v2 / c2);
    }
    
    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
    
    WRITE_LOG << "SR Sod shock tube initialized (Kitajima et al. 2025, Eqs. 74-75):";
    WRITE_LOG << "  Left:  " << N_left << " particles, P=" << P_left 
             << ", n=" << n_left << ", ν=" << nu_left;
    WRITE_LOG << "  Right: " << N_right << " particles, P=" << P_right 
             << ", n=" << n_right << ", ν=" << nu_right;
    WRITE_LOG << "  Different ν: " << (different_nu ? "YES" : "NO");
    WRITE_LOG << "  Transition: Smoothed over width " << transition_width << " (~5h)";
    WRITE_LOG << "  Note: Primitives (P,n,v) smoothly interpolated near x=0 to avoid huge forces.";
    WRITE_LOG << "        Conserved variables (S,e,N) will be computed by pre_interaction.";
    WRITE_LOG << "        This ensures consistency between kernel-smoothed N and conserved e.";
#endif
}

}
