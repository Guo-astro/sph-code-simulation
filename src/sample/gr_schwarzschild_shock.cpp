/**
 * Schwarzschild Radial Shock Tube
 *
 * 1D shock tube aligned with the radial direction in Schwarzschild spacetime.
 * Tests the gravitational source terms in GR-GSPH.
 *
 * Based on Liptai & Price (2019) GRSPH benchmarks.
 *
 * Key physics:
 * - Gravitational redshift affects wave speeds
 * - Outgoing waves are slowed (redshifted)
 * - Ingoing waves are accelerated (blueshifted)
 * - Comparison with Minkowski case shows GR effects
 *
 * Domain: r ∈ [r_min, r_max] where r_min > 2M (outside horizon)
 * Discontinuity at r = r_disc
 *
 * Initial conditions (Rosswog Test 1 in curved spacetime):
 *   Left (r < r_disc):  (n, v, P) = (10, 0, 40/3)
 *   Right (r > r_disc): (n, v, P) = (1, 0, 1e-6)
 */

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

void Solver::make_gr_schwarzschild_shock()
{
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real gamma = m_param->physics.gamma;
    const real c_speed = m_param->srgsph.c_speed;
    const real c2 = c_speed * c_speed;

    // Black hole mass from config (or default)
    real M = 1.0;
    if (m_sample_parameters.count("blackHoleMass")) {
        M = boost::any_cast<real>(m_sample_parameters["blackHoleMass"]);
    }

    // Domain in units of M (Schwarzschild radius = 2M)
    // Must be outside horizon: r > 2M
    real r_min = 3.0 * M;  // Inner edge (safe margin above horizon)
    real r_max = 15.0 * M; // Outer edge
    real r_disc = 6.0 * M; // Discontinuity location

    if (m_sample_parameters.count("rMin")) {
        r_min = boost::any_cast<real>(m_sample_parameters["rMin"]);
    }
    if (m_sample_parameters.count("rMax")) {
        r_max = boost::any_cast<real>(m_sample_parameters["rMax"]);
    }
    if (m_sample_parameters.count("rDiscontinuity")) {
        r_disc = boost::any_cast<real>(m_sample_parameters["rDiscontinuity"]);
    }

    // Validate domain
    const real r_horizon = 2.0 * M;
    if (r_min <= r_horizon) {
        THROW_ERROR("Schwarzschild shock tube: r_min must be > 2M (horizon)");
    }
    if (r_disc <= r_min || r_disc >= r_max) {
        THROW_ERROR("Schwarzschild shock tube: r_disc must be in (r_min, r_max)");
    }

    WRITE_LOG << "GR-GSPH Schwarzschild Radial Shock Tube:";
    WRITE_LOG << "  Black hole mass M = " << M;
    WRITE_LOG << "  Schwarzschild radius r_s = " << r_horizon;
    WRITE_LOG << "  Domain: r ∈ [" << r_min << ", " << r_max << "]";
    WRITE_LOG << "  Discontinuity at r = " << r_disc;

    // ============================================================================
    // Initial conditions: Rosswog Test 1
    // ============================================================================
    const real n_left = 10.0;
    const real v_left = 0.0;
    const real P_left = 40.0 / 3.0;

    const real n_right = 1.0;
    const real v_right = 0.0;
    const real P_right = 1.0e-6;

    WRITE_LOG << "  Left state (r < " << r_disc << "): n=" << n_left << ", v=" << v_left << ", P=" << P_left;
    WRITE_LOG << "  Right state (r > " << r_disc << "): n=" << n_right << ", v=" << v_right << ", P=" << P_right;

    // ============================================================================
    // Particle distribution
    // ============================================================================
    // Split particles between left and right based on density ratio
    const real L_left = r_disc - r_min;
    const real L_right = r_max - r_disc;
    const real L_total = L_left + L_right;

    // Equal mass per particle means more particles in denser region
    const real mass_left = n_left * L_left;
    const real mass_right = n_right * L_right;
    const real mass_total = mass_left + mass_right;

    const int N_left = static_cast<int>(N * mass_left / mass_total);
    const int N_right = N - N_left;

    WRITE_LOG << "  Particles: " << N_left << " (left) + " << N_right << " (right) = " << N;

    // Particle spacing
    const real dx_left = L_left / N_left;
    const real dx_right = L_right / N_right;

    // Total baryon number (mass in SR/GR units)
    const real nu = (mass_left + mass_right) / N;

    WRITE_LOG << "  Baryon number per particle nu = " << nu;

    // ============================================================================
    // Create particles
    // ============================================================================
    auto& particles = m_sim->get_particles();
    particles.clear();
    particles.reserve(N);

    int particle_id = 0;

    // Left particles (high density region, r < r_disc)
    for (int i = 0; i < N_left; ++i) {
        SPHParticle p;
        p.id = particle_id++;

        // Position (1D along x-axis = radial direction)
        const real r = r_min + (i + 0.5) * dx_left;
        p.pos[0] = r;
#if DIM >= 2
        p.pos[1] = 0.0;
#endif
#if DIM == 3
        p.pos[2] = 0.0;
#endif

        // Primitive variables
        p.dens = n_left;       // Rest-frame number density
        p.pres = P_left;
        p.vel[0] = v_left;     // Radial velocity (Eulerian V^r)
#if DIM >= 2
        p.vel[1] = 0.0;
#endif
#if DIM == 3
        p.vel[2] = 0.0;
#endif
        p.vel_t = 0.0;  // No tangent velocity

        // Baryon number (constant per particle)
        p.nu = nu;

        // Specific internal energy: u = P / ((gamma-1) * n)
        p.ene = P_left / ((gamma - 1.0) * n_left);

        // Sound speed (relativistic): c_s² = γ P / (ρ h)
        // where h = 1 + ε + P/ρ = 1 + γ/(γ-1) * P/n for ideal gas
        const real h_left = 1.0 + gamma / (gamma - 1.0) * P_left / n_left;
        p.sound = std::sqrt(gamma * P_left / (n_left * h_left));

        // Lorentz factor (initially 1.0 since v=0)
        p.gamma_lor = 1.0;

        // Lab-frame density N = n * Γ
        p.N = n_left * p.gamma_lor;

        // Initial smoothing length estimate
        p.sml = 1.2 * dx_left;

        p.is_ghost = false;

        particles.push_back(p);
    }

    // Right particles (low density region, r > r_disc)
    for (int i = 0; i < N_right; ++i) {
        SPHParticle p;
        p.id = particle_id++;

        // Position
        const real r = r_disc + (i + 0.5) * dx_right;
        p.pos[0] = r;
#if DIM >= 2
        p.pos[1] = 0.0;
#endif
#if DIM == 3
        p.pos[2] = 0.0;
#endif

        // Primitive variables
        p.dens = n_right;
        p.pres = P_right;
        p.vel[0] = v_right;
#if DIM >= 2
        p.vel[1] = 0.0;
#endif
#if DIM == 3
        p.vel[2] = 0.0;
#endif
        p.vel_t = 0.0;

        p.nu = nu;
        p.ene = P_right / ((gamma - 1.0) * n_right);

        const real h_right = 1.0 + gamma / (gamma - 1.0) * P_right / n_right;
        p.sound = std::sqrt(gamma * P_right / (n_right * h_right));

        p.gamma_lor = 1.0;
        p.N = n_right * p.gamma_lor;
        p.sml = 1.2 * dx_right;
        p.is_ghost = false;

        particles.push_back(p);
    }

    m_sim->set_particle_num(static_cast<int>(particles.size()));

    // ============================================================================
    // Add ghost particles at boundaries
    // ============================================================================
    if (m_param->periodic.is_valid == false) {
        // Get ghost layers from config
        int ghost_layers = 6;
        if (m_sample_parameters.count("ghostLayers")) {
            ghost_layers = boost::any_cast<int>(m_sample_parameters["ghostLayers"]);
        }

        const int num_real = m_sim->get_particle_num();

        // Inner boundary ghosts (near r_min)
        for (int i = 0; i < ghost_layers; ++i) {
            SPHParticle ghost;
            ghost.id = particle_id++;

            // Mirror position inside r_min
            const real r = r_min - (i + 0.5) * dx_left;
            ghost.pos[0] = r;
#if DIM >= 2
            ghost.pos[1] = 0.0;
#endif
#if DIM == 3
            ghost.pos[2] = 0.0;
#endif

            // Copy left state
            ghost.dens = n_left;
            ghost.pres = P_left;
            ghost.vel[0] = -v_left;  // Reflect velocity (inflow boundary)
#if DIM >= 2
            ghost.vel[1] = 0.0;
#endif
#if DIM == 3
            ghost.vel[2] = 0.0;
#endif
            ghost.vel_t = 0.0;
            ghost.nu = nu;
            ghost.ene = P_left / ((gamma - 1.0) * n_left);
            ghost.sound = particles[0].sound;
            ghost.gamma_lor = 1.0;
            ghost.N = n_left;
            ghost.sml = dx_left;
            ghost.is_ghost = true;

            particles.push_back(ghost);
        }

        // Outer boundary ghosts (near r_max)
        for (int i = 0; i < ghost_layers; ++i) {
            SPHParticle ghost;
            ghost.id = particle_id++;

            const real r = r_max + (i + 0.5) * dx_right;
            ghost.pos[0] = r;
#if DIM >= 2
            ghost.pos[1] = 0.0;
#endif
#if DIM == 3
            ghost.pos[2] = 0.0;
#endif

            // Copy right state (outflow boundary)
            ghost.dens = n_right;
            ghost.pres = P_right;
            ghost.vel[0] = v_right;
#if DIM >= 2
            ghost.vel[1] = 0.0;
#endif
#if DIM == 3
            ghost.vel[2] = 0.0;
#endif
            ghost.vel_t = 0.0;
            ghost.nu = nu;
            ghost.ene = P_right / ((gamma - 1.0) * n_right);
            ghost.sound = particles[num_real - 1].sound;
            ghost.gamma_lor = 1.0;
            ghost.N = n_right;
            ghost.sml = dx_right;
            ghost.is_ghost = true;

            particles.push_back(ghost);
        }

        WRITE_LOG << "  Added " << (2 * ghost_layers) << " ghost particles at boundaries";
    }

    // Store sample parameters for output
    m_sample_parameters["r_min"] = r_min;
    m_sample_parameters["r_max"] = r_max;
    m_sample_parameters["r_disc"] = r_disc;
    m_sample_parameters["n_left"] = n_left;
    m_sample_parameters["n_right"] = n_right;
    m_sample_parameters["P_left"] = P_left;
    m_sample_parameters["P_right"] = P_right;

    WRITE_LOG << "Schwarzschild shock tube initialized with " << particles.size() << " total particles";
}

} // namespace sph
