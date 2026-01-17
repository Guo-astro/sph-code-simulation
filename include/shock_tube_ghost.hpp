#pragma once

/**
 * @file shock_tube_ghost.hpp
 * @brief Ghost particle management for shock tube boundary conditions
 *
 * Provides fixed boundary conditions for shock tube tests by creating
 * ghost particles at domain boundaries that maintain constant properties.
 *
 * For the Sod shock tube:
 * - Left boundary (x < x_min): Ghost particles with LEFT state (high density/pressure)
 * - Right boundary (x > x_max): Ghost particles with RIGHT state (low density/pressure)
 */

#include "particle.hpp"
#include "vector_type.hpp"
#include <vector>

namespace sph {

/**
 * @brief Configuration for shock tube ghost boundaries
 */
struct ShockTubeGhostConfig {
    // Domain boundaries
    real x_min = 0.0;
    real x_max = 1.0;

    // Left state (x < discontinuity)
    real rho_left = 1.0;
    real pres_left = 1.0;
    real vel_left = 0.0;

    // Right state (x > discontinuity)
    real rho_right = 0.125;
    real pres_right = 0.1;
    real vel_right = 0.0;

    // Discontinuity position (for reference)
    real x_discontinuity = 0.5;

    // Thermodynamic parameter
    real gamma = 1.4;

    // Ghost particle parameters
    int n_ghost_layers = 3;  // Number of ghost particle layers at each boundary
};

/**
 * @brief Create ghost particles for shock tube boundaries
 *
 * Adds ghost particle layers at x < x_min and x > x_max with fixed
 * properties corresponding to the left and right states respectively.
 *
 * @param particles Vector of existing real particles (modified to add ghosts)
 * @param config Ghost boundary configuration
 * @param dx Particle spacing in x direction
 * @param dy Particle spacing in y direction (for 2D/3D)
 * @param dz Particle spacing in z direction (for 3D)
 * @param Ly Domain size in y
 * @param Lz Domain size in z
 */
void create_shock_tube_ghost_particles(
    std::vector<SPHParticle>& particles,
    const ShockTubeGhostConfig& config,
    real dx, real dy, real dz,
    real Ly, real Lz);

/**
 * @brief Simplified interface for creating ghost particles
 */
void create_shock_tube_ghost_particles(
    std::vector<SPHParticle>& particles,
    int n_ghost_layers,
    real dx, real dy, real dz,
    real Ly, real Lz,
    real gamma);

/**
 * @brief Update/restore ghost particle properties
 *
 * Called during simulation to ensure ghost particles maintain their
 * fixed boundary condition properties. Ghost positions are not modified
 * (position fixing is done by skipping ghosts in the integrator).
 *
 * @param particles Vector of particles
 * @param config Ghost boundary configuration
 */
void update_shock_tube_ghost_particles(
    std::vector<SPHParticle>& particles,
    const ShockTubeGhostConfig& config);

/**
 * @brief Simplified interface for updating ghost particles
 */
void update_shock_tube_ghost_particles(
    std::vector<SPHParticle>& particles,
    real gamma);

/**
 * @brief Check if boundary ghost particles are needed for this sample type
 *
 * @param sample_name Name of the sample (e.g., "shock_tube", "shock_tube_2d", "shock_tube_3d")
 * @return true if ghost boundaries should be used
 */
bool needs_shock_tube_ghost_boundaries(const std::string& sample_name);

} // namespace sph
