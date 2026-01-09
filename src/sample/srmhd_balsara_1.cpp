#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"

#include <cmath>

namespace sph
{

/**
 * SRMHD Balsara Test 1: Relativistic MHD Shock Tube
 *
 * Based on Balsara (2001), ApJS, 132, 83 - Test Problem 1
 * A standard benchmark for relativistic MHD codes.
 *
 * Initial conditions (Balsara Test 1):
 * Left (x < 0):
 *   rho = 1.0, P = 1.0, vx = 0, vy = 0, vz = 0
 *   Bx = 0.5, By = 1.0, Bz = 0
 *
 * Right (x > 0):
 *   rho = 0.125, P = 0.1, vx = 0, vy = 0, vz = 0
 *   Bx = 0.5, By = -1.0, Bz = 0
 *
 * Parameters:
 *   gamma = 2.0 (relativistic adiabatic index)
 *   c = 1.0 (speed of light)
 *
 * Domain: [-0.5, 0.5]
 * End time: t = 0.4
 *
 * The solution consists of:
 * - Left-going fast rarefaction
 * - Left-going slow shock
 * - Contact discontinuity
 * - Right-going slow shock
 * - Right-going fast rarefaction
 *
 * This test checks:
 * - Relativistic fast and slow MHD waves
 * - Tangential B-field discontinuities
 * - Contact discontinuity tracking
 */
void Solver::make_srmhd_balsara_1()
{
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real gamma = m_param->physics.gamma;
    const real c_speed = m_param->srgsph.c_speed;

    // Left state (Balsara Test 1)
    const real rho_L = 1.0;
    const real P_L = 1.0;
    const real vx_L = 0.0;
    const real vy_L = 0.0;
    const real vz_L = 0.0;

    // Right state
    const real rho_R = 0.125;
    const real P_R = 0.1;
    const real vx_R = 0.0;
    const real vy_R = 0.0;
    const real vz_R = 0.0;

    // Magnetic field - read from config, default to Balsara values
    real Bx_L = 0.5, By_L = 1.0, Bz_L = 0.0;
    real Bx_R = 0.5, By_R = -1.0, Bz_R = 0.0;
    if (m_sample_parameters.count("useMHD")) {
        bool use_mhd = boost::any_cast<bool>(m_sample_parameters["useMHD"]);
        if (!use_mhd) {
            // Pure hydro test (B=0)
            Bx_L = Bx_R = By_L = By_R = Bz_L = Bz_R = 0.0;
            std::cout << "[SRMHD] Pure hydro test (B=0)" << std::endl;
        }
    }

    // Domain
    const real x_min = -0.5;
    const real x_max = 0.5;
    const real x_interface = 0.0;

    // Particle spacing
    const real domain_length = x_max - x_min;
    const real dx_R = domain_length / N;
    const real dx_L = dx_R * (rho_R / rho_L);  // Adjust for density ratio

    // Smoothing length: use dx directly like SR-GSPH (sr_sod.cpp)
    // The pre-interaction module will iterate to find proper h if enabled
    const real h_L = dx_L;
    const real h_R = dx_R;

    // Ghost particle configuration
    bool use_ghost_particles = true;
    if (m_sample_parameters.count("useGhostParticles")) {
        use_ghost_particles = boost::any_cast<bool>(m_sample_parameters["useGhostParticles"]);
    }

    int ghost_layers = 6;
    if (m_sample_parameters.count("ghostLayers")) {
        ghost_layers = boost::any_cast<int>(m_sample_parameters["ghostLayers"]);
    }

    const int n_ghost_layers = use_ghost_particles ? ghost_layers : 0;

    // Estimate particle count
    const real L_left = x_interface - x_min;
    const real L_right = x_max - x_interface;
    const int N_left = static_cast<int>(L_left / dx_L) + 1;
    const int N_right = static_cast<int>(L_right / dx_R) + 1;
    const int N_ghost_left = n_ghost_layers;
    const int N_ghost_right = n_ghost_layers;
    const int num_total = N_left + N_right + N_ghost_left + N_ghost_right;

    std::vector<SPHParticle> particles(num_total);

    // Mass per particle (1D mass = rho * dx)
    // For SR, mass is baryon number nu
    const real nu_L = rho_L * dx_L;  // Baryon number per particle
    const real nu_R = rho_R * dx_R;

    int idx = 0;

    // Compute Lorentz factor (all particles initially at rest)
    auto compute_lorentz = [&c_speed](real vx, real vy, real vz) {
        const real v2 = vx * vx + vy * vy + vz * vz;
        return 1.0 / std::sqrt(1.0 - v2 / (c_speed * c_speed));
    };

    // Compute enthalpy
    auto compute_enthalpy = [&gamma, &c_speed](real rho, real P) {
        const real u = P / ((gamma - 1.0) * rho);
        return 1.0 + u / (c_speed * c_speed) + P / (rho * c_speed * c_speed);
    };

    // Left region (x < 0)
    real x = x_min + dx_L * 0.5;
    while (x < x_interface && idx < num_total) {
        auto& p = particles[idx];

        const real gamma_lor = compute_lorentz(vx_L, vy_L, vz_L);
        const real H = compute_enthalpy(rho_L, P_L);
        const real N_lab = gamma_lor * rho_L;  // Lab-frame density

#if DIM == 1
        p.pos = vec_t(x);
        p.vel = vec_t(vx_L);
        p.S = p.vel * (gamma_lor * H);
#elif DIM == 2
        p.pos = vec_t(x, 0.0);
        p.vel = vec_t(vx_L, vy_L);
        p.S = p.vel * (gamma_lor * H);
#else
        p.pos = vec_t(x, 0.0, 0.0);
        p.vel = vec_t(vx_L, vy_L, vz_L);
        p.S = p.vel * (gamma_lor * H);
#endif
        p.vy_mhd = vy_L;
        p.vz_mhd = vz_L;
        p.dens = rho_L;
        p.pres = P_L;
        p.nu = nu_L;
        p.mass = nu_L;
        p.N = N_lab;
        p.ene = P_L / ((gamma - 1.0) * rho_L);
        p.gamma_lor = gamma_lor;
        p.enthalpy = H;
        p.B = vec3_t(Bx_L, By_L, Bz_L);
        p.sml = h_L;
        p.id = idx;

        // Canonical energy (from SR-GSPH formulation)
        const real X_L = gamma / (gamma - 1.0);
        p.e = (H * (X_L * gamma_lor * gamma_lor - 1.0) + 1.0) / (X_L * gamma_lor);

        // Initialize derivatives to zero (critical for first timestep)
        p.dS = vec_t(0.0);
        p.de = 0.0;
        p.dS_t = 0.0;
        p.dS_old = vec_t(0.0);
        p.de_old = 0.0;
        p.dS_t_old = 0.0;
        p.dB = vec3_t(0.0, 0.0, 0.0);
        p.dB_old = vec3_t(0.0, 0.0, 0.0);

        ++idx;
        x += dx_L;
    }

    // Right region (x >= 0)
    x = x_interface + dx_R * 0.5;
    while (x < x_max && idx < num_total) {
        auto& p = particles[idx];

        const real gamma_lor = compute_lorentz(vx_R, vy_R, vz_R);
        const real H = compute_enthalpy(rho_R, P_R);
        const real N_lab = gamma_lor * rho_R;

#if DIM == 1
        p.pos = vec_t(x);
        p.vel = vec_t(vx_R);
        p.S = p.vel * (gamma_lor * H);
#elif DIM == 2
        p.pos = vec_t(x, 0.0);
        p.vel = vec_t(vx_R, vy_R);
        p.S = p.vel * (gamma_lor * H);
#else
        p.pos = vec_t(x, 0.0, 0.0);
        p.vel = vec_t(vx_R, vy_R, vz_R);
        p.S = p.vel * (gamma_lor * H);
#endif
        p.vy_mhd = vy_R;
        p.vz_mhd = vz_R;
        p.dens = rho_R;
        p.pres = P_R;
        p.nu = nu_R;
        p.mass = nu_R;
        p.N = N_lab;
        p.ene = P_R / ((gamma - 1.0) * rho_R);
        p.gamma_lor = gamma_lor;
        p.enthalpy = H;
        p.B = vec3_t(Bx_R, By_R, Bz_R);
        p.sml = h_R;
        p.id = idx;

        // Canonical energy (from SR-GSPH formulation)
        const real X_R = gamma / (gamma - 1.0);
        p.e = (H * (X_R * gamma_lor * gamma_lor - 1.0) + 1.0) / (X_R * gamma_lor);

        // Initialize derivatives to zero (critical for first timestep)
        p.dS = vec_t(0.0);
        p.de = 0.0;
        p.dS_t = 0.0;
        p.dS_old = vec_t(0.0);
        p.de_old = 0.0;
        p.dS_t_old = 0.0;
        p.dB = vec3_t(0.0, 0.0, 0.0);
        p.dB_old = vec3_t(0.0, 0.0, 0.0);

        ++idx;
        x += dx_R;
    }

    // Add left ghost particles
    for (int g = 0; g < n_ghost_layers && idx < num_total; ++g) {
        auto& p = particles[idx];
        real ghost_x = x_min - dx_L * (0.5 + g);

        const real gamma_lor = compute_lorentz(vx_L, vy_L, vz_L);
        const real H = compute_enthalpy(rho_L, P_L);
        const real N_lab = gamma_lor * rho_L;

#if DIM == 1
        p.pos = vec_t(ghost_x);
        p.vel = vec_t(vx_L);
        p.S = p.vel * (gamma_lor * H);
#elif DIM == 2
        p.pos = vec_t(ghost_x, 0.0);
        p.vel = vec_t(vx_L, vy_L);
        p.S = p.vel * (gamma_lor * H);
#else
        p.pos = vec_t(ghost_x, 0.0, 0.0);
        p.vel = vec_t(vx_L, vy_L, vz_L);
        p.S = p.vel * (gamma_lor * H);
#endif
        p.vy_mhd = vy_L;
        p.vz_mhd = vz_L;
        p.dens = rho_L;
        p.pres = P_L;
        p.nu = nu_L;
        p.mass = nu_L;
        p.N = N_lab;
        p.ene = P_L / ((gamma - 1.0) * rho_L);
        p.gamma_lor = gamma_lor;
        p.enthalpy = H;
        p.B = vec3_t(Bx_L, By_L, Bz_L);
        p.sml = h_L;
        p.is_ghost = true;
        p.id = idx;

        // Canonical energy for ghost
        const real X_gL = gamma / (gamma - 1.0);
        p.e = (H * (X_gL * gamma_lor * gamma_lor - 1.0) + 1.0) / (X_gL * gamma_lor);

        // Initialize derivatives to zero
        p.dS = vec_t(0.0);
        p.de = 0.0;
        p.dS_t = 0.0;
        p.dS_old = vec_t(0.0);
        p.de_old = 0.0;
        p.dS_t_old = 0.0;
        p.dB = vec3_t(0.0, 0.0, 0.0);
        p.dB_old = vec3_t(0.0, 0.0, 0.0);

        ++idx;
    }

    // Add right ghost particles
    for (int g = 0; g < n_ghost_layers && idx < num_total; ++g) {
        auto& p = particles[idx];
        real ghost_x = x_max + dx_R * (0.5 + g);

        const real gamma_lor = compute_lorentz(vx_R, vy_R, vz_R);
        const real H = compute_enthalpy(rho_R, P_R);
        const real N_lab = gamma_lor * rho_R;

#if DIM == 1
        p.pos = vec_t(ghost_x);
        p.vel = vec_t(vx_R);
        p.S = p.vel * (gamma_lor * H);
#elif DIM == 2
        p.pos = vec_t(ghost_x, 0.0);
        p.vel = vec_t(vx_R, vy_R);
        p.S = p.vel * (gamma_lor * H);
#else
        p.pos = vec_t(ghost_x, 0.0, 0.0);
        p.vel = vec_t(vx_R, vy_R, vz_R);
        p.S = p.vel * (gamma_lor * H);
#endif
        p.vy_mhd = vy_R;
        p.vz_mhd = vz_R;
        p.dens = rho_R;
        p.pres = P_R;
        p.nu = nu_R;
        p.mass = nu_R;
        p.N = N_lab;
        p.ene = P_R / ((gamma - 1.0) * rho_R);
        p.gamma_lor = gamma_lor;
        p.enthalpy = H;
        p.B = vec3_t(Bx_R, By_R, Bz_R);
        p.sml = h_R;
        p.is_ghost = true;
        p.id = idx;

        // Canonical energy for ghost
        const real X_gR = gamma / (gamma - 1.0);
        p.e = (H * (X_gR * gamma_lor * gamma_lor - 1.0) + 1.0) / (X_gR * gamma_lor);

        // Initialize derivatives to zero
        p.dS = vec_t(0.0);
        p.de = 0.0;
        p.dS_t = 0.0;
        p.dS_old = vec_t(0.0);
        p.de_old = 0.0;
        p.dS_t_old = 0.0;
        p.dB = vec3_t(0.0, 0.0, 0.0);
        p.dB_old = vec3_t(0.0, 0.0, 0.0);

        ++idx;
    }

    // Resize to actual count
    particles.resize(idx);

    int n_real = 0, n_ghost = 0;
    for (const auto& p : particles) {
        if (p.is_ghost) ++n_ghost;
        else ++n_real;
    }

    std::cout << "[SRMHD Balsara Test 1]" << std::endl;
    std::cout << "  Real particles: " << n_real << std::endl;
    std::cout << "  Ghost particles: " << n_ghost << std::endl;
    std::cout << "  Total particles: " << idx << std::endl;
    std::cout << "  Left:  N=" << N_left << ", dx=" << dx_L << std::endl;
    std::cout << "  Right: N=" << N_right << ", dx=" << dx_R << std::endl;
    std::cout << "  Bx = " << Bx_L << " (constant)" << std::endl;
    std::cout << "  By_L = " << By_L << ", By_R = " << By_R << std::endl;
    std::cout << "  c = " << c_speed << ", gamma = " << gamma << std::endl;

    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());
}

}
