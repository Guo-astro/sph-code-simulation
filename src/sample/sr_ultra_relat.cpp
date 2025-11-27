// Ultra-Relativistic Shock Tests
// Based on Kitajima et al. (2025) arXiv:2510.18251v1 Section 3.2
//
// Tests extreme Lorentz factors by colliding a relativistic flow
// with a stationary fluid.
//
// Paper Initial Conditions (Section 3.2):
//   Left:  (P, n, v^x) = (1.0, 1.0, v)  where v = 0.9 to 0.999999999
//   Right: (P, n, v^x) = (1.0, 1.0, 0)
//   gamma_c = 5/3
//   1000 particles on each side
//   t = 0.3

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

void Solver::make_sr_ultra_relat()
{
#if DIM != 1
    THROW_ERROR("Ultra-relativistic shock test requires DIM == 1");
#else

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real gamma = m_param->physics.gamma;
    const real c_speed = m_param->srgsph.c_speed;
    const real c2 = c_speed * c_speed;

    // Velocity of left state (default: 0.9c)
    // Paper tests: v = 0.9, 0.99, 0.999, 0.9999, 0.999999999
    real v_left_input = 0.9;
    if (m_sample_parameters.count("v_left")) {
        v_left_input = boost::any_cast<real>(m_sample_parameters["v_left"]);
    }

    // ============================================================================
    // Ghost particle configuration
    // Ghost particles extend beyond boundaries to provide proper neighbor support
    // ============================================================================
    bool use_ghost_particles = true;  // Enable ghost particles for better conservation
    if (m_sample_parameters.count("useGhostParticles")) {
        use_ghost_particles = boost::any_cast<bool>(m_sample_parameters["useGhostParticles"]);
    }
    
    // Number of ghost particle layers
    // Need enough layers so boundary particles have same neighbor count as interior
    // Search radius is 6*h, and h ~ dx, so we need ~6 layers on each side
    int ghost_layers = 6;  // Default: 6 layers on each side for full neighbor support
    if (m_sample_parameters.count("ghostLayers")) {
        ghost_layers = boost::any_cast<int>(m_sample_parameters["ghostLayers"]);
    }

    // ============================================================================
    // Initial conditions from paper Section 3.2 (lines 592-595)
    // Left:  (P, n, v^x) = (1.0, 1.0, v)  - moving at relativistic speed
    // Right: (P, n, v^x) = (1.0, 1.0, 0)  - at rest
    // ============================================================================
    const real P_left = 1.0;   // Paper uses P=1.0 (not 10.0!)
    const real n_left = 1.0;
    const real v_left = v_left_input * c_speed;

    const real P_right = 1.0;  // Paper uses P=1.0 (not 10.0!)
    const real n_right = 1.0;
    const real v_right = 0.0;

    // Paper: 1000 particles on each side (line 605)
    const int N_left = N;
    const int N_right = N;
    const int N_ghost_left = use_ghost_particles ? ghost_layers : 0;
    const int N_ghost_right = use_ghost_particles ? ghost_layers : 0;
    const int num = N_left + N_right + N_ghost_left + N_ghost_right;

    // Domain: x ∈ [-0.5, 0.5], discontinuity at x=0
    const real x_left_start = -0.5;
    const real x_left_end = 0.0;
    const real x_right_start = 0.0;
    const real x_right_end = 0.5;

    const real dx_left = (x_left_end - x_left_start) / N_left;
    const real dx_right = (x_right_end - x_right_start) / N_right;

    // Baryon number per particle (volume-based approach)
    const real nu_left = n_left * dx_left;
    const real nu_right = n_right * dx_right;

    // Compute Lorentz factor for left state
    const real gamma_lor_left = 1.0 / std::sqrt(1.0 - v_left * v_left / c2);

    std::vector<SPHParticle> p(num);

    // ============================================================================
    // Initialize left state (relativistic flow)
    // ============================================================================
    for (int i = 0; i < N_left; ++i) {
        auto& p_i = p[i];
        p_i.id = i;
        p_i.pos[0] = x_left_start + (i + 0.5) * dx_left;

        p_i.nu = nu_left;
        p_i.mass = nu_left;

        vec_t vel;
        vel[0] = v_left;

        p_i.gamma_lor = gamma_lor_left;

        // Thermodynamic quantities
        const real u_init = P_left / ((gamma - 1.0) * n_left);
        const real H_init = 1.0 + u_init / c2 + P_left / (n_left * c2);
        p_i.sound = std::sqrt((gamma - 1.0) * (H_init - 1.0) / H_init) * c_speed;
        p_i.enthalpy = H_init;

        // Conserved variables
        const real N_conserved = p_i.gamma_lor * n_left;
        const vec_t S_conserved = vel * (p_i.gamma_lor * H_init);
        const real X = gamma / (gamma - 1.0);
        const real e_conserved = (H_init * (X * p_i.gamma_lor * p_i.gamma_lor - 1.0) + 1.0) / (X * p_i.gamma_lor);

        p_i.N = N_conserved;
        p_i.S = S_conserved;
        p_i.e = e_conserved;
        p_i.vel = vel;
        p_i.ene = u_init;
        p_i.pres = P_left;
        p_i.dens = N_conserved;

        // Initialize derivatives
        p_i.dS = vec_t(0.0);
        p_i.de = 0.0;
        p_i.dS_old = vec_t(0.0);
        p_i.de_old = 0.0;

        p_i.sml = dx_left;
    }

    // ============================================================================
    // Initialize right state (at rest)
    // ============================================================================
    for (int i = 0; i < N_right; ++i) {
        auto& p_i = p[N_left + i];
        p_i.id = N_left + i;
        p_i.pos[0] = x_right_start + (i + 0.5) * dx_right;

        p_i.nu = nu_right;
        p_i.mass = nu_right;

        vec_t vel;
        vel[0] = v_right;

        p_i.gamma_lor = 1.0;

        // Thermodynamic quantities
        const real u_init = P_right / ((gamma - 1.0) * n_right);
        const real H_init = 1.0 + u_init / c2 + P_right / (n_right * c2);
        p_i.sound = std::sqrt((gamma - 1.0) * (H_init - 1.0) / H_init) * c_speed;
        p_i.enthalpy = H_init;

        // Conserved variables
        const real N_conserved = p_i.gamma_lor * n_right;
        const vec_t S_conserved = vel * (p_i.gamma_lor * H_init);
        const real X = gamma / (gamma - 1.0);
        const real e_conserved = (H_init * (X * p_i.gamma_lor * p_i.gamma_lor - 1.0) + 1.0) / (X * p_i.gamma_lor);

        p_i.N = N_conserved;
        p_i.S = S_conserved;
        p_i.e = e_conserved;
        p_i.vel = vel;
        p_i.ene = u_init;
        p_i.pres = P_right;
        p_i.dens = N_conserved;

        // Initialize derivatives
        p_i.dS = vec_t(0.0);
        p_i.de = 0.0;
        p_i.dS_old = vec_t(0.0);
        p_i.de_old = 0.0;

        p_i.sml = dx_right;
    }

    // ============================================================================
    // Initialize ghost particles (left boundary - extend left state)
    // Ghost particles mirror the boundary conditions for proper neighbor support
    // ============================================================================
    if (use_ghost_particles) {
        int ghost_id = N_left + N_right;
        
        // Left ghost particles (extend left state beyond x = -0.5)
        for (int i = 0; i < N_ghost_left; ++i) {
            auto& p_i = p[ghost_id];
            p_i.id = ghost_id;
            // Place ghost particles beyond left boundary, mirroring spacing
            p_i.pos[0] = x_left_start - (i + 0.5) * dx_left;
            p_i.is_ghost = true;  // Mark as ghost particle

            // Copy left state properties
            p_i.nu = nu_left;
            p_i.mass = nu_left;

            vec_t vel;
            vel[0] = -v_left;  // Reflect velocity for wall condition

            const real v2 = vel[0] * vel[0];
            p_i.gamma_lor = 1.0 / std::sqrt(1.0 - v2 / c2);

            const real u_init = P_left / ((gamma - 1.0) * n_left);
            const real H_init = 1.0 + u_init / c2 + P_left / (n_left * c2);
            p_i.sound = std::sqrt((gamma - 1.0) * (H_init - 1.0) / H_init) * c_speed;
            p_i.enthalpy = H_init;

            const real N_conserved = p_i.gamma_lor * n_left;
            const vec_t S_conserved = vel * (p_i.gamma_lor * H_init);
            const real X = gamma / (gamma - 1.0);
            const real e_conserved = (H_init * (X * p_i.gamma_lor * p_i.gamma_lor - 1.0) + 1.0) / (X * p_i.gamma_lor);

            p_i.N = N_conserved;
            p_i.S = S_conserved;
            p_i.e = e_conserved;
            p_i.vel = vel;
            p_i.ene = u_init;
            p_i.pres = P_left;
            p_i.dens = N_conserved;

            p_i.dS = vec_t(0.0);
            p_i.de = 0.0;
            p_i.dS_old = vec_t(0.0);
            p_i.de_old = 0.0;
            p_i.sml = dx_left;

            ghost_id++;
        }

        // Right ghost particles (extend right state beyond x = 0.5)
        for (int i = 0; i < N_ghost_right; ++i) {
            auto& p_i = p[ghost_id];
            p_i.id = ghost_id;
            // Place ghost particles beyond right boundary
            p_i.pos[0] = x_right_end + (i + 0.5) * dx_right;
            p_i.is_ghost = true;  // Mark as ghost particle

            // Copy right state properties
            p_i.nu = nu_right;
            p_i.mass = nu_right;

            vec_t vel;
            vel[0] = -v_right;  // Reflect velocity for wall condition (v_right is 0)

            const real v2 = vel[0] * vel[0];
            p_i.gamma_lor = 1.0 / std::sqrt(1.0 - v2 / c2);

            const real u_init = P_right / ((gamma - 1.0) * n_right);
            const real H_init = 1.0 + u_init / c2 + P_right / (n_right * c2);
            p_i.sound = std::sqrt((gamma - 1.0) * (H_init - 1.0) / H_init) * c_speed;
            p_i.enthalpy = H_init;

            const real N_conserved = p_i.gamma_lor * n_right;
            const vec_t S_conserved = vel * (p_i.gamma_lor * H_init);
            const real X = gamma / (gamma - 1.0);
            const real e_conserved = (H_init * (X * p_i.gamma_lor * p_i.gamma_lor - 1.0) + 1.0) / (X * p_i.gamma_lor);

            p_i.N = N_conserved;
            p_i.S = S_conserved;
            p_i.e = e_conserved;
            p_i.vel = vel;
            p_i.ene = u_init;
            p_i.pres = P_right;
            p_i.dens = N_conserved;

            p_i.dS = vec_t(0.0);
            p_i.de = 0.0;
            p_i.dS_old = vec_t(0.0);
            p_i.de_old = 0.0;
            p_i.sml = dx_right;

            ghost_id++;
        }
    }

    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());

    WRITE_LOG << "Ultra-relativistic shock initialized (Kitajima Section 3.2):";
    WRITE_LOG << "  Left:  v=" << v_left / c_speed << "c, gamma=" << gamma_lor_left
              << ", P=" << P_left << ", n=" << n_left;
    WRITE_LOG << "  Right: v=0, gamma=1, P=" << P_right << ", n=" << n_right;
    WRITE_LOG << "  Particles: " << N_left << " left + " << N_right << " right";
    if (use_ghost_particles) {
        WRITE_LOG << "  Ghost: " << (N_ghost_left + N_ghost_right) << " particles (" 
                  << N_ghost_left << " left, " << N_ghost_right << " right)";
    }
#endif
}

}
