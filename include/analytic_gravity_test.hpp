#pragma once

#include <vector>
#include "defines.hpp"
#include "vector_type.hpp"
#include "particle.hpp"

/**
 * @file analytic_gravity_test.hpp
 * @brief Test utilities for validating gravity implementations against analytic solutions
 *
 * Provides functions to:
 * 1. Create uniform sphere particle distributions
 * 2. Compute direct summation gravity with various softening kernels
 * 3. Compare numerical results against analytic uniform sphere solutions
 *
 * TDD STUB: These declarations exist only to allow tests to compile.
 * Implementation is pending (TDD red phase).
 */

namespace sph {

/**
 * @enum GravitySofteningType
 * @brief Types of softening kernels for gravity calculation
 */
enum class GravitySofteningType {
    HERNQUIST_KATZ,  // Cubic spline (Hernquist & Katz 1989)
    WENDLAND_C4      // Wendland C4 kernel
};

/**
 * @struct DirectGravityResult
 * @brief Result of direct summation gravity calculation
 */
struct DirectGravityResult {
    vec_t acceleration;  // Gravitational acceleration vector
    real potential;      // Gravitational potential (scalar)
};

/**
 * @brief Create uniform density sphere particle distribution
 * @param N Number of particles
 * @param R Sphere radius
 * @param rho Density (total mass M = 4/3 * pi * rho * R^3)
 * @param seed Random seed for reproducibility
 * @return Vector of SPH particles with random positions inside sphere
 *
 * Uses rejection sampling or radial inversion to ensure uniform density.
 * Each particle has mass = M/N.
 *
 * TDD STUB: Returns empty vector to cause test failures
 */
std::vector<SPHParticle> create_uniform_sphere_distribution(
    int N, double R, double rho, unsigned seed = 42);

/**
 * @brief Compute direct summation gravity at a test position
 * @param pos Test position
 * @param particles Source particles
 * @param softening Softening length (epsilon or h depending on kernel)
 * @param type Softening kernel type
 * @return DirectGravityResult containing acceleration and potential
 *
 * This is O(N) for a single test position.
 * Uses the specified softening kernel for close encounters.
 *
 * TDD STUB: Returns zero acceleration to cause test failures
 */
DirectGravityResult compute_direct_gravity(
    const vec_t& pos,
    const std::vector<SPHParticle>& particles,
    real softening,
    GravitySofteningType type);

/**
 * @brief Compute direct summation gravity for all particles
 * @param particles Particles (modified in place with grav_acc and phi)
 * @param softening Softening length
 * @param type Softening kernel type
 * @param G Gravitational constant
 *
 * This is O(N^2) for all particles.
 *
 * TDD STUB: Does nothing to cause test failures
 */
void compute_direct_gravity_all(
    std::vector<SPHParticle>& particles,
    real softening,
    GravitySofteningType type,
    real G = 1.0);

} // namespace sph
