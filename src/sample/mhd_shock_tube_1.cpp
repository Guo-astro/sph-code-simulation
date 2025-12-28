#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"

#include <cmath>

namespace sph
{

/**
 * MHD Shock Tube 1: Dai-Woodward test (Iwasaki & Inutsuka 2011, Fig. 5)
 *
 * This is the standard MHD shock tube test from Dai & Woodward (1994) and Ryu & Jones (1995).
 * The exact solution consists of:
 * - Two fast shocks
 * - Two rotational discontinuities (Alfven waves)
 * - Two slow shocks
 * - One contact discontinuity
 *
 * Initial conditions:
 * Left (x < 0):
 *   rho = 1.08, P = 0.95, vx = 1.2, vy = 0.01, vz = 0.5
 *   By = 3.6/sqrt(4*pi), Bz = 2/sqrt(4*pi)
 *
 * Right (x > 0):
 *   rho = 1.0, P = 1.0, vx = 0, vy = 0, vz = 0
 *   By = 4/sqrt(4*pi), Bz = 2/sqrt(4*pi)
 *
 * Common: Bx = 2/sqrt(4*pi), gamma = 5/3
 *
 * Domain: [-0.74, 0.5]
 * End time: t = 0.2
 *
 * Note: Uses vec3_t for B field and separate vy_mhd, vz_mhd for velocity
 * components to support DIM=1 builds with full 3D MHD physics.
 */
void Solver::make_mhd_shock_tube_1()
{
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real gamma = m_param->physics.gamma;

    // Normalization: sqrt(4*pi) for Gaussian units -> 1 for our code
    const real sqrt_4pi = std::sqrt(4.0 * M_PI);

    // Left state (Dai-Woodward)
    const real rho_L = 1.08;
    const real P_L = 0.95;
    const real vx_L = 1.2;
    const real vy_L = 0.01;
    const real vz_L = 0.5;
    const real Bx_L = 2.0 / sqrt_4pi;
    const real By_L = 3.6 / sqrt_4pi;
    const real Bz_L = 2.0 / sqrt_4pi;

    // Right state
    const real rho_R = 1.0;
    const real P_R = 1.0;
    const real vx_R = 0.0;
    const real vy_R = 0.0;
    const real vz_R = 0.0;
    const real Bx_R = 2.0 / sqrt_4pi;
    const real By_R = 4.0 / sqrt_4pi;
    const real Bz_R = 2.0 / sqrt_4pi;

    // Domain (1D line along x-axis)
    const real x_min = -0.74;
    const real x_max = 0.5;
    const real x_interface = 0.0;

    // Particle spacing
    const real domain_length = x_max - x_min;
    const real dx_R = domain_length / N;
    const real dx_L = dx_R * (rho_R / rho_L);  // Adjust for density ratio

    // Ghost particle configuration from config
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
    const real mass_L = rho_L * dx_L;
    const real mass_R = rho_R * dx_R;

    int idx = 0;

    // Left region (x < 0)
    real x = x_min + dx_L * 0.5;
    while (x < x_interface && idx < num_total) {
        auto& p = particles[idx];
#if DIM == 1
        p.pos = vec_t(x);
        p.vel = vec_t(vx_L);
#elif DIM == 2
        p.pos = vec_t(x, 0.0);
        p.vel = vec_t(vx_L, vy_L);
#else
        p.pos = vec_t(x, 0.0, 0.0);
        p.vel = vec_t(vx_L, vy_L, vz_L);
#endif
        p.vy_mhd = vy_L;
        p.vz_mhd = vz_L;
        p.dens = rho_L;
        p.pres = P_L;
        p.mass = mass_L;
        p.ene = P_L / ((gamma - 1.0) * rho_L);
        p.B = vec3_t(Bx_L, By_L, Bz_L);
        p.id = idx;
        ++idx;
        x += dx_L;
    }

    // Right region (x >= 0)
    x = x_interface + dx_R * 0.5;
    while (x < x_max && idx < num_total) {
        auto& p = particles[idx];
#if DIM == 1
        p.pos = vec_t(x);
        p.vel = vec_t(vx_R);
#elif DIM == 2
        p.pos = vec_t(x, 0.0);
        p.vel = vec_t(vx_R, vy_R);
#else
        p.pos = vec_t(x, 0.0, 0.0);
        p.vel = vec_t(vx_R, vy_R, vz_R);
#endif
        p.vy_mhd = vy_R;
        p.vz_mhd = vz_R;
        p.dens = rho_R;
        p.pres = P_R;
        p.mass = mass_R;
        p.ene = P_R / ((gamma - 1.0) * rho_R);
        p.B = vec3_t(Bx_R, By_R, Bz_R);
        p.id = idx;
        ++idx;
        x += dx_R;
    }

    // Add left ghost particles (extend leftward from x_min)
    for (int g = 0; g < n_ghost_layers && idx < num_total; ++g) {
        auto& p = particles[idx];
        real ghost_x = x_min - dx_L * (0.5 + g);
#if DIM == 1
        p.pos = vec_t(ghost_x);
        p.vel = vec_t(vx_L);
#elif DIM == 2
        p.pos = vec_t(ghost_x, 0.0);
        p.vel = vec_t(vx_L, vy_L);
#else
        p.pos = vec_t(ghost_x, 0.0, 0.0);
        p.vel = vec_t(vx_L, vy_L, vz_L);
#endif
        p.vy_mhd = vy_L;
        p.vz_mhd = vz_L;
        p.dens = rho_L;
        p.pres = P_L;
        p.mass = mass_L;
        p.ene = P_L / ((gamma - 1.0) * rho_L);
        p.B = vec3_t(Bx_L, By_L, Bz_L);
        p.sml = dx_L;  // Set initial smoothing length
        p.is_ghost = true;  // Mark as ghost particle
        p.id = idx;
        ++idx;
    }

    // Add right ghost particles (extend rightward from x_max)
    for (int g = 0; g < n_ghost_layers && idx < num_total; ++g) {
        auto& p = particles[idx];
        real ghost_x = x_max + dx_R * (0.5 + g);
#if DIM == 1
        p.pos = vec_t(ghost_x);
        p.vel = vec_t(vx_R);
#elif DIM == 2
        p.pos = vec_t(ghost_x, 0.0);
        p.vel = vec_t(vx_R, vy_R);
#else
        p.pos = vec_t(ghost_x, 0.0, 0.0);
        p.vel = vec_t(vx_R, vy_R, vz_R);
#endif
        p.vy_mhd = vy_R;
        p.vz_mhd = vz_R;
        p.dens = rho_R;
        p.pres = P_R;
        p.mass = mass_R;
        p.ene = P_R / ((gamma - 1.0) * rho_R);
        p.B = vec3_t(Bx_R, By_R, Bz_R);
        p.sml = dx_R;  // Set initial smoothing length
        p.is_ghost = true;  // Mark as ghost particle
        p.id = idx;
        ++idx;
    }

    // Resize to actual count
    particles.resize(idx);

    int n_real = 0, n_ghost = 0;
    for (const auto& p : particles) {
        if (p.is_ghost) ++n_ghost;
        else ++n_real;
    }

    std::cout << "[MHD Shock Tube 1 - Dai-Woodward]" << std::endl;
    std::cout << "  Real particles: " << n_real << std::endl;
    std::cout << "  Ghost particles: " << n_ghost << std::endl;
    std::cout << "  Total particles: " << idx << std::endl;
    std::cout << "  Left:  N=" << N_left << ", dx=" << dx_L << std::endl;
    std::cout << "  Right: N=" << N_right << ", dx=" << dx_R << std::endl;
    std::cout << "  Bx = " << Bx_L << " (constant)" << std::endl;
    std::cout << "  By_L = " << By_L << ", By_R = " << By_R << std::endl;

    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());
}

}
