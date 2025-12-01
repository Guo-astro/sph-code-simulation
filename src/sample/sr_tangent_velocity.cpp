// Special Relativistic Shock Tube with Tangential Velocity
// Implementation of Kitajima et al. (2025) arXiv:2510.18251v1 Section 3.1.5
//
// Riemann Problem 5: One-Dimensional Shock Tube Problem With Tangential Velocity
//
// In relativistic hydrodynamics, the velocity component along the partition v^t
// affects the solution, unlike non-relativistic cases (Pons et al. 2000).
// When the initial tangential velocity is large, errors in numerical solution
// become large and extremely high resolution is required.
//
// IMPORTANT: This is a 1D problem with DIM=1. The tangent velocity v^t is tracked
// as a separate scalar field (particle.vel_t), NOT as a second velocity component.
// The SPH kernel and density calculations remain 1D.
//
// Initial conditions (strong blast with tangent velocity):
//   Left:  (P, n, v^x, v^t) = (1000, 1.0, 0, {0, 0.9, 0.99})
//   Right: (P, n, v^x, v^t) = (0.01, 1.0, 0, {0, 0.9, 0.99})
//
// The simulation is stopped when the shock front reaches x = 0.2
// (not a fixed end time, as shock velocity depends on tangent velocity)

#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"
#include "logger.hpp"
#include "srgsph/sr_primitive_recovery.hpp"
#include <cmath>
#include <string>

namespace sph
{

void Solver::make_sr_tangent_velocity()
{
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real gamma = m_param->physics.gamma;
    const real c_speed = m_param->srgsph.c_speed;
    const real c2 = c_speed * c_speed;

    // ============================================================================
    // Tangent velocity configuration
    // Following Kitajima et al. (2025) Section 3.1.5 (Riemann Problem 5)
    // ============================================================================
    
    // Tangent velocities v^t_L and v^t_R (default: 0.9)
    real vt_left = 0.9;
    real vt_right = 0.9;
    
    if (m_sample_parameters.count("vt_left")) {
        vt_left = boost::any_cast<real>(m_sample_parameters["vt_left"]);
    }
    if (m_sample_parameters.count("vt_right")) {
        vt_right = boost::any_cast<real>(m_sample_parameters["vt_right"]);
    }

    // ============================================================================
    // Initial conditions from Kitajima Section 3.1.5 (Eqs. 94-95)
    // Same as strong blast, but with tangential velocity
    // ============================================================================
    
    // Normal velocity v^x is zero for this test
    const real vx_left = 0.0;
    const real vx_right = 0.0;
    
    // Strong blast initial conditions
    const real P_left = 1000.0;
    const real n_left = 1.0;
    const real P_right = 0.01;
    const real n_right = 1.0;

    // Equal density on both sides with equal particle count
    // Following paper: 1600 particles each side for standard test
    const int N_left = N;
    const int N_right = N;

    WRITE_LOG << "Test type: TANGENT VELOCITY (Section 3.1.5)";
    WRITE_LOG << "  v^t_L = " << vt_left << ", v^t_R = " << vt_right;

    // ============================================================================
    // Ghost particle configuration
    // ============================================================================
    bool use_ghost_particles = true;
    if (m_sample_parameters.count("useGhostParticles")) {
        use_ghost_particles = boost::any_cast<bool>(m_sample_parameters["useGhostParticles"]);
    }
    
    int ghost_layers = 6;
    if (m_sample_parameters.count("ghostLayers")) {
        ghost_layers = boost::any_cast<int>(m_sample_parameters["ghostLayers"]);
    }

    const int N_ghost_left = use_ghost_particles ? ghost_layers : 0;
    const int N_ghost_right = use_ghost_particles ? ghost_layers : 0;
    const int num = N_left + N_right + N_ghost_left + N_ghost_right;

    // Domain boundaries from config (rangeMin/rangeMax), discontinuity at x=0
    const real x_left_start = m_param->periodic.range_min[0];
    const real x_left_end = 0.0;
    const real x_right_start = 0.0;
    const real x_right_end = m_param->periodic.range_max[0];

    const real dx_left = (x_left_end - x_left_start) / N_left;
    const real dx_right = (x_right_end - x_right_start) / N_right;

    // ============================================================================
    // Lorentz factors including tangential velocity
    // γ = 1/√(1 - (v_x² + v_t²)/c²)
    // ============================================================================
    const real v2_left = vx_left * vx_left + vt_left * vt_left;
    const real v2_right = vx_right * vx_right + vt_right * vt_right;
    
    // Check for superluminal velocities
    if (v2_left >= c2) {
        THROW_ERROR("Left state velocity exceeds speed of light: v² = " + std::to_string(v2_left));
    }
    if (v2_right >= c2) {
        THROW_ERROR("Right state velocity exceeds speed of light: v² = " + std::to_string(v2_right));
    }
    
    const real gamma_left = 1.0 / std::sqrt(1.0 - v2_left / c2);
    const real gamma_right = 1.0 / std::sqrt(1.0 - v2_right / c2);

    // Baryon number per particle (ν = γ * n * dx)
    const real nu_left = gamma_left * n_left * dx_left;
    const real nu_right = gamma_right * n_right * dx_right;

    std::vector<SPHParticle> p(num);

    // ============================================================================
    // Initialize left state particles
    // ============================================================================
    for (int i = 0; i < N_left; ++i) {
        auto& p_i = p[i];
        p_i.id = i;
        
        // Position: 1D along x-axis
        p_i.pos[0] = x_left_start + (i + 0.5) * dx_left;
#if DIM >= 2
        p_i.pos[1] = 0.0;
#endif
#if DIM == 3
        p_i.pos[2] = 0.0;
#endif

        // Baryon number per particle
        p_i.nu = nu_left;
        p_i.mass = nu_left;

        // Velocity: v^x = vx_left (in vel[0]), v^t = vt_left (in vel_t scalar)
        p_i.vel[0] = vx_left;
#if DIM >= 2
        p_i.vel[1] = 0.0;  // Don't put vt here - that changes the physics
#endif
#if DIM == 3
        p_i.vel[2] = 0.0;
#endif
        p_i.vel_t = vt_left;  // Store tangent velocity in separate field

        // Lorentz factor from total velocity (including v_t)
        p_i.gamma_lor = gamma_left;

        // Thermodynamic quantities
        const real u_init = P_left / ((gamma - 1.0) * n_left);
        const real H_init = 1.0 + u_init / c2 + P_left / (n_left * c2);
        p_i.sound = std::sqrt((gamma - 1.0) * (H_init - 1.0) / H_init) * c_speed;
        p_i.enthalpy = H_init;

        // Conserved variables (including tangent momentum)
        const real N_conserved = p_i.gamma_lor * n_left;
        p_i.N = N_conserved;
        
        // S = γ H v (normal component only in vel)
        p_i.S[0] = vx_left * (p_i.gamma_lor * H_init);
#if DIM >= 2
        p_i.S[1] = 0.0;
#endif
#if DIM == 3
        p_i.S[2] = 0.0;
#endif
        // Tangent momentum as separate field
        p_i.S_t = vt_left * (p_i.gamma_lor * H_init);
        
        const real X = gamma / (gamma - 1.0);
        p_i.e = (H_init * (X * p_i.gamma_lor * p_i.gamma_lor - 1.0) + 1.0) / (X * p_i.gamma_lor);
        
        p_i.ene = u_init;
        p_i.pres = P_left;
        p_i.dens = n_left;  // REST-FRAME density (not conserved N = γ*n)

        // Initialize derivatives
        p_i.dS = vec_t(0.0);
        p_i.de = 0.0;
        p_i.dS_old = vec_t(0.0);
        p_i.de_old = 0.0;
        p_i.dS_t = 0.0;
        p_i.dS_t_old = 0.0;

        // Smoothing length initial guess
        p_i.sml = dx_left;
    }

    // ============================================================================
    // Initialize right state particles
    // ============================================================================
    for (int i = 0; i < N_right; ++i) {
        auto& p_i = p[N_left + i];
        p_i.id = N_left + i;
        
        // Position: 1D along x-axis
        p_i.pos[0] = x_right_start + (i + 0.5) * dx_right;
#if DIM >= 2
        p_i.pos[1] = 0.0;
#endif
#if DIM == 3
        p_i.pos[2] = 0.0;
#endif

        // Baryon number per particle
        p_i.nu = nu_right;
        p_i.mass = nu_right;

        // Velocity
        p_i.vel[0] = vx_right;
#if DIM >= 2
        p_i.vel[1] = 0.0;
#endif
#if DIM == 3
        p_i.vel[2] = 0.0;
#endif
        p_i.vel_t = vt_right;

        // Lorentz factor from total velocity
        p_i.gamma_lor = gamma_right;

        // Thermodynamic quantities
        const real u_init = P_right / ((gamma - 1.0) * n_right);
        const real H_init = 1.0 + u_init / c2 + P_right / (n_right * c2);
        p_i.sound = std::sqrt((gamma - 1.0) * (H_init - 1.0) / H_init) * c_speed;
        p_i.enthalpy = H_init;

        // Conserved variables
        const real N_conserved = p_i.gamma_lor * n_right;
        p_i.N = N_conserved;
        
        p_i.S[0] = vx_right * (p_i.gamma_lor * H_init);
#if DIM >= 2
        p_i.S[1] = 0.0;
#endif
#if DIM == 3
        p_i.S[2] = 0.0;
#endif
        p_i.S_t = vt_right * (p_i.gamma_lor * H_init);
        
        const real X = gamma / (gamma - 1.0);
        p_i.e = (H_init * (X * p_i.gamma_lor * p_i.gamma_lor - 1.0) + 1.0) / (X * p_i.gamma_lor);
        
        p_i.ene = u_init;
        p_i.pres = P_right;
        p_i.dens = n_right;  // REST-FRAME density (not conserved N = γ*n)

        // Initialize derivatives
        p_i.dS = vec_t(0.0);
        p_i.de = 0.0;
        p_i.dS_old = vec_t(0.0);
        p_i.de_old = 0.0;
        p_i.dS_t = 0.0;
        p_i.dS_t_old = 0.0;

        // Smoothing length initial guess
        p_i.sml = dx_right;
    }

    // ============================================================================
    // Initialize ghost particles
    // ============================================================================
    if (use_ghost_particles) {
        int ghost_id = N_left + N_right;
        
        // Left ghost particles
        for (int i = 0; i < N_ghost_left; ++i) {
            auto& p_i = p[ghost_id];
            p_i.id = ghost_id;
            p_i.pos[0] = x_left_start - (i + 0.5) * dx_left;
#if DIM >= 2
            p_i.pos[1] = 0.0;
#endif
#if DIM == 3
            p_i.pos[2] = 0.0;
#endif
            p_i.is_ghost = true;

            p_i.nu = nu_left;
            p_i.mass = nu_left;

            // Reflecting boundary: negate normal velocity, keep tangent
            p_i.vel[0] = -vx_left;
#if DIM >= 2
            p_i.vel[1] = 0.0;
#endif
#if DIM == 3
            p_i.vel[2] = 0.0;
#endif
            p_i.vel_t = vt_left;

            p_i.gamma_lor = gamma_left;

            const real u_init = P_left / ((gamma - 1.0) * n_left);
            const real H_init = 1.0 + u_init / c2 + P_left / (n_left * c2);
            p_i.sound = std::sqrt((gamma - 1.0) * (H_init - 1.0) / H_init) * c_speed;
            p_i.enthalpy = H_init;

            const real N_conserved = p_i.gamma_lor * n_left;
            p_i.N = N_conserved;
            p_i.S[0] = (-vx_left) * (p_i.gamma_lor * H_init);
#if DIM >= 2
            p_i.S[1] = 0.0;
#endif
#if DIM == 3
            p_i.S[2] = 0.0;
#endif
            p_i.S_t = vt_left * (p_i.gamma_lor * H_init);
            
            const real X = gamma / (gamma - 1.0);
            p_i.e = (H_init * (X * p_i.gamma_lor * p_i.gamma_lor - 1.0) + 1.0) / (X * p_i.gamma_lor);
            
            p_i.ene = u_init;
            p_i.pres = P_left;
            p_i.dens = n_left;  // REST-FRAME density

            p_i.dS = vec_t(0.0);
            p_i.de = 0.0;
            p_i.dS_old = vec_t(0.0);
            p_i.de_old = 0.0;
            p_i.dS_t = 0.0;
            p_i.dS_t_old = 0.0;
            p_i.sml = dx_left;

            ghost_id++;
        }

        // Right ghost particles
        for (int i = 0; i < N_ghost_right; ++i) {
            auto& p_i = p[ghost_id];
            p_i.id = ghost_id;
            p_i.pos[0] = x_right_end + (i + 0.5) * dx_right;
#if DIM >= 2
            p_i.pos[1] = 0.0;
#endif
#if DIM == 3
            p_i.pos[2] = 0.0;
#endif
            p_i.is_ghost = true;

            p_i.nu = nu_right;
            p_i.mass = nu_right;

            // Reflecting boundary
            p_i.vel[0] = -vx_right;
#if DIM >= 2
            p_i.vel[1] = 0.0;
#endif
#if DIM == 3
            p_i.vel[2] = 0.0;
#endif
            p_i.vel_t = vt_right;

            p_i.gamma_lor = gamma_right;

            const real u_init = P_right / ((gamma - 1.0) * n_right);
            const real H_init = 1.0 + u_init / c2 + P_right / (n_right * c2);
            p_i.sound = std::sqrt((gamma - 1.0) * (H_init - 1.0) / H_init) * c_speed;
            p_i.enthalpy = H_init;

            const real N_conserved = p_i.gamma_lor * n_right;
            p_i.N = N_conserved;
            p_i.S[0] = (-vx_right) * (p_i.gamma_lor * H_init);
#if DIM >= 2
            p_i.S[1] = 0.0;
#endif
#if DIM == 3
            p_i.S[2] = 0.0;
#endif
            p_i.S_t = vt_right * (p_i.gamma_lor * H_init);
            
            const real X = gamma / (gamma - 1.0);
            p_i.e = (H_init * (X * p_i.gamma_lor * p_i.gamma_lor - 1.0) + 1.0) / (X * p_i.gamma_lor);
            
            p_i.ene = u_init;
            p_i.pres = P_right;
            p_i.dens = n_right;  // REST-FRAME density

            p_i.dS = vec_t(0.0);
            p_i.de = 0.0;
            p_i.dS_old = vec_t(0.0);
            p_i.de_old = 0.0;
            p_i.dS_t = 0.0;
            p_i.dS_t_old = 0.0;
            p_i.sml = dx_right;

            ghost_id++;
        }
    }

    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());

    WRITE_LOG << "SR Tangent Velocity shock tube initialized (Kitajima et al. 2025 Section 3.1.5):";
    WRITE_LOG << "  Left:  " << N_left << " particles, P=" << P_left << ", n=" << n_left;
    WRITE_LOG << "         v^x=" << vx_left << ", v^t=" << vt_left << ", γ=" << gamma_left;
    WRITE_LOG << "  Right: " << N_right << " particles, P=" << P_right << ", n=" << n_right;
    WRITE_LOG << "         v^x=" << vx_right << ", v^t=" << vt_right << ", γ=" << gamma_right;
    if (use_ghost_particles) {
        WRITE_LOG << "  Ghost: " << (N_ghost_left + N_ghost_right) << " particles";
    }
    WRITE_LOG << "  Baryon numbers: nu_L=" << nu_left << ", nu_R=" << nu_right;
    WRITE_LOG << "  dx_left=" << dx_left << ", dx_right=" << dx_right;
}

}
