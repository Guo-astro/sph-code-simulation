/**
 * @file shock_tube_ghost.cpp
 * @brief Implementation of ghost particle management for shock tube boundaries
 */

#include "shock_tube_ghost.hpp"
#include <cmath>
#include <algorithm>

namespace sph {

// Default configuration for Sod shock tube
static ShockTubeGhostConfig default_config;

void create_shock_tube_ghost_particles(
    std::vector<SPHParticle>& particles,
    const ShockTubeGhostConfig& config,
    real dx, real dy, real dz,
    real Ly, real Lz)
{
    // Calculate number of particles in y and z directions
    int Ny = static_cast<int>(std::round(Ly / dy));
    int Nz = static_cast<int>(std::round(Lz / dz));

    if (Ny < 1) Ny = 1;
    if (Nz < 1) Nz = 1;

    // Number of ghost particles at each boundary
    int ghosts_per_boundary = config.n_ghost_layers * Ny * Nz;

    // Reserve space for new ghost particles
    size_t original_size = particles.size();
    particles.reserve(original_size + 2 * ghosts_per_boundary);

    // Find the maximum existing particle ID
    int max_id = 0;
    for (const auto& p : particles) {
        if (p.id > max_id) max_id = p.id;
    }
    int next_id = max_id + 1;

    // Create LEFT boundary ghosts (x < x_min)
    // These have LEFT state properties (high density/pressure)
    for (int layer = 0; layer < config.n_ghost_layers; ++layer) {
        real x = config.x_min - (layer + 0.5) * dx;

        for (int iz = 0; iz < Nz; ++iz) {
            for (int iy = 0; iy < Ny; ++iy) {
                SPHParticle ghost;

                ghost.pos[0] = x;
#if DIM >= 2
                ghost.pos[1] = (iy + 0.5) * dy;
#endif
#if DIM >= 3
                ghost.pos[2] = (iz + 0.5) * dz;
#endif

                ghost.vel = vec_t(0.0);
                ghost.vel[0] = config.vel_left;

                ghost.dens = config.rho_left;
                ghost.pres = config.pres_left;
                ghost.mass = config.rho_left * dx * dy * dz;
                ghost.ene = config.pres_left / ((config.gamma - 1.0) * config.rho_left);
                ghost.sound = std::sqrt(config.gamma * config.pres_left / config.rho_left);

                ghost.id = next_id++;
                ghost.is_ghost = true;
                // Set smoothing length based on particle spacing and expected neighbor count
                // sml = (N_neighbor * mass / (dens * A))^(1/DIM) where A is volume factor
                constexpr real A = DIM == 1 ? 2.0 : DIM == 2 ? M_PI : 4.0 * M_PI / 3.0;
                constexpr int N_neighbor = 50;  // Default neighbor number
                ghost.sml = std::pow(N_neighbor * ghost.mass / (ghost.dens * A), 1.0 / DIM);
                ghost.acc = vec_t(0.0);

                particles.push_back(ghost);
            }
        }
    }

    // Create RIGHT boundary ghosts (x > x_max)
    // These have RIGHT state properties (low density/pressure)
    for (int layer = 0; layer < config.n_ghost_layers; ++layer) {
        real x = config.x_max + (layer + 0.5) * dx;

        for (int iz = 0; iz < Nz; ++iz) {
            for (int iy = 0; iy < Ny; ++iy) {
                SPHParticle ghost;

                ghost.pos[0] = x;
#if DIM >= 2
                ghost.pos[1] = (iy + 0.5) * dy;
#endif
#if DIM >= 3
                ghost.pos[2] = (iz + 0.5) * dz;
#endif

                ghost.vel = vec_t(0.0);
                ghost.vel[0] = config.vel_right;

                ghost.dens = config.rho_right;
                ghost.pres = config.pres_right;
                ghost.mass = config.rho_right * dx * dy * dz;
                ghost.ene = config.pres_right / ((config.gamma - 1.0) * config.rho_right);
                ghost.sound = std::sqrt(config.gamma * config.pres_right / config.rho_right);

                ghost.id = next_id++;
                ghost.is_ghost = true;
                // Set smoothing length based on particle spacing and expected neighbor count
                constexpr real A = DIM == 1 ? 2.0 : DIM == 2 ? M_PI : 4.0 * M_PI / 3.0;
                constexpr int N_neighbor = 50;  // Default neighbor number
                ghost.sml = std::pow(N_neighbor * ghost.mass / (ghost.dens * A), 1.0 / DIM);
                ghost.acc = vec_t(0.0);

                particles.push_back(ghost);
            }
        }
    }
}

void create_shock_tube_ghost_particles(
    std::vector<SPHParticle>& particles,
    int n_ghost_layers,
    real dx, real dy, real dz,
    real Ly, real Lz,
    real gamma)
{
    ShockTubeGhostConfig config;
    config.n_ghost_layers = n_ghost_layers;
    config.gamma = gamma;
    // Use default Sod shock tube values for other parameters

    create_shock_tube_ghost_particles(particles, config, dx, dy, dz, Ly, Lz);
}

void update_shock_tube_ghost_particles(
    std::vector<SPHParticle>& particles,
    const ShockTubeGhostConfig& config)
{
    for (auto& p : particles) {
        if (!p.is_ghost) continue;

        // Determine if this is a left or right ghost based on position
        if (p.pos[0] < config.x_min) {
            // Left ghost - restore LEFT state
            p.dens = config.rho_left;
            p.pres = config.pres_left;
            p.vel = vec_t(0.0);
            p.vel[0] = config.vel_left;
            p.ene = config.pres_left / ((config.gamma - 1.0) * config.rho_left);
            p.sound = std::sqrt(config.gamma * config.pres_left / config.rho_left);
            p.acc = vec_t(0.0);
        } else if (p.pos[0] > config.x_max) {
            // Right ghost - restore RIGHT state
            p.dens = config.rho_right;
            p.pres = config.pres_right;
            p.vel = vec_t(0.0);
            p.vel[0] = config.vel_right;
            p.ene = config.pres_right / ((config.gamma - 1.0) * config.rho_right);
            p.sound = std::sqrt(config.gamma * config.pres_right / config.rho_right);
            p.acc = vec_t(0.0);
        }
    }
}

void update_shock_tube_ghost_particles(
    std::vector<SPHParticle>& particles,
    real gamma)
{
    ShockTubeGhostConfig config;
    config.gamma = gamma;

    update_shock_tube_ghost_particles(particles, config);
}

bool needs_shock_tube_ghost_boundaries(const std::string& sample_name)
{
    return (sample_name == "shock_tube" ||
            sample_name == "shock_tube_2d" ||
            sample_name == "shock_tube_3d");
}

} // namespace sph
