#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"

#include <cmath>

namespace sph
{

/**
 * MHD Shock Tube 2: Strong shock test (Iwasaki & Inutsuka 2011, Fig. 6)
 *
 * This is a more challenging MHD shock tube test from Dai & Woodward (1994),
 * Ryu & Jones (1995), and Toth (2000).
 * The exact solution consists of:
 * - Two fast shocks
 * - A left-propagating slow rarefaction wave
 * - A right-propagating slow shock
 * - One contact discontinuity
 *
 * Initial conditions:
 * Left (x < 0):
 *   rho = 1, P = 20, vx = 10, vy = 0, vz = 0
 *   By = 5/sqrt(4*pi), Bz = 0
 *
 * Right (x > 0):
 *   rho = 1, P = 1, vx = -10, vy = 0, vz = 0
 *   By = 5/sqrt(4*pi), Bz = 0
 *
 * Common: Bx = 5/sqrt(4*pi), gamma = 5/3
 *
 * Domain: [-1, 1]
 * End time: t = 0.06
 */
void Solver::make_mhd_shock_tube_2()
{
    // MHD shock tubes require 3D vectors (Bx, By, Bz and vx, vy, vz)
    // even though spatial variation is only in x direction
#if DIM != 3
    THROW_ERROR("MHD Shock Tube 2 requires DIM == 3 (for full B and v vectors)");
#else
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real gamma = m_param->physics.gamma;

    // Normalization
    const real sqrt_4pi = std::sqrt(4.0 * M_PI);

    // Left state (Strong shock)
    const real rho_L = 1.0;
    const real P_L = 20.0;
    const real vx_L = 10.0;
    const real vy_L = 0.0;
    const real vz_L = 0.0;
    const real Bx_L = 5.0 / sqrt_4pi;
    const real By_L = 5.0 / sqrt_4pi;
    const real Bz_L = 0.0;

    // Right state
    const real rho_R = 1.0;
    const real P_R = 1.0;
    const real vx_R = -10.0;
    const real vy_R = 0.0;
    const real vz_R = 0.0;
    const real Bx_R = 5.0 / sqrt_4pi;
    const real By_R = 5.0 / sqrt_4pi;
    const real Bz_R = 0.0;

    // Domain (1D line along x-axis)
    const real x_min = -1.0;
    const real x_max = 1.0;
    const real x_interface = 0.0;

    // Equal particle spacing (same density on both sides)
    const real domain_length = x_max - x_min;
    const real dx = domain_length / N;

    // Count particles on each side
    const real L_left = x_interface - x_min;
    const real L_right = x_max - x_interface;
    const int N_left = static_cast<int>(L_left / dx) + 1;
    const int N_right = static_cast<int>(L_right / dx) + 1;
    const int num_total = N_left + N_right;

    std::vector<SPHParticle> particles(num_total);

    // Mass per particle (1D mass = rho * dx)
    const real mass = rho_L * dx;  // Same density everywhere

    int idx = 0;

    // Left region (x < 0)
    real x = x_min + dx * 0.5;
    while (x < x_interface && idx < num_total) {
        auto& p = particles[idx];
        p.pos = vec_t{x, 0.0, 0.0};
        p.vel = vec_t{vx_L, vy_L, vz_L};
        p.dens = rho_L;
        p.pres = P_L;
        p.mass = mass;
        p.ene = P_L / ((gamma - 1.0) * rho_L);
        p.B = vec3_t{Bx_L, By_L, Bz_L};
        p.id = idx;
        ++idx;
        x += dx;
    }

    // Right region (x >= 0)
    x = x_interface + dx * 0.5;
    while (x < x_max && idx < num_total) {
        auto& p = particles[idx];
        p.pos = vec_t{x, 0.0, 0.0};
        p.vel = vec_t{vx_R, vy_R, vz_R};
        p.dens = rho_R;
        p.pres = P_R;
        p.mass = mass;
        p.ene = P_R / ((gamma - 1.0) * rho_R);
        p.B = vec3_t{Bx_R, By_R, Bz_R};
        p.id = idx;
        ++idx;
        x += dx;
    }

    // Resize to actual count
    particles.resize(idx);

    std::cout << "[MHD Shock Tube 2 - Strong Shock]" << std::endl;
    std::cout << "  Total particles: " << idx << std::endl;
    std::cout << "  dx = " << dx << std::endl;
    std::cout << "  Bx = By = " << Bx_L << " (constant)" << std::endl;
    std::cout << "  Colliding velocities: vx = +/- " << std::abs(vx_L) << std::endl;

    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());
#endif
}

}
