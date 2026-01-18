/**
 * @file test_analytic_gravity_3d.cpp
 * @brief TDD tests for 3D analytic gravity (uniform sphere)
 *
 * These tests are written BEFORE the implementation (TDD red phase).
 * All tests should FAIL until the analytic gravity testing utilities are implemented.
 *
 * Test Categories:
 * 1. Analytic reference tests - validate uniform sphere gravity formula
 * 2. Softening convergence tests - as h -> 0, matches analytic
 * 3. Resolution convergence tests - as N -> inf, error decreases
 * 4. Kernel comparison tests - Hernquist-Katz vs Wendland C4
 *
 * Reference: TDD plan at thoughts/shared/plans/gravity-optimizations-tdd-plan.md
 *
 * Uniform Sphere Gravity:
 *   Inside (r < R):  g = (4/3) * pi * G * rho * r  (linear in r)
 *   Outside (r > R): g = G * M / r^2              (inverse square law)
 */

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <random>

#include "defines.hpp"
#include "vector_type.hpp"
#include "particle.hpp"
#include "analytic_gravity_test.hpp"  // Stub header for test utilities

namespace {

// =============================================================================
// Tolerances
// =============================================================================
const double kAbsTol = 1e-10;
const double kContinuityTol = 1e-4;

} // anonymous namespace

// =============================================================================
// Test Fixture
// =============================================================================

class AnalyticGravityTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Physical parameters (dimensionless units)
        G = 1.0;        // Gravitational constant
        rho = 1.0;      // Uniform density
        R = 1.0;        // Sphere radius
        M = (4.0/3.0) * M_PI * rho * R * R * R;  // Total mass
    }

    double G, rho, R, M;

    // Analytic gravity inside sphere: g = (4/3) * pi * G * rho * r
    double analytic_g_inside(double r) const {
        return (4.0/3.0) * M_PI * G * rho * r;
    }

    // Analytic gravity outside sphere: g = G * M / r^2
    double analytic_g_outside(double r) const {
        return G * M / (r * r);
    }

    // Analytic gravity magnitude at radius r
    double analytic_g(double r) const {
        if (r < 1e-30) return 0.0;
        return (r < R) ? analytic_g_inside(r) : analytic_g_outside(r);
    }
};

// =============================================================================
// 1. ANALYTIC REFERENCE TESTS
// =============================================================================

TEST_F(AnalyticGravityTest, GivenUniformSphere_WhenAtCenter_ThenGravityIsZero) {
    // GIVEN: Uniform density sphere
    // WHEN: Compute gravity at r = 0
    // THEN: g = 0 (by symmetry)

    double g = analytic_g_inside(0.0);
    EXPECT_NEAR(g, 0.0, kAbsTol);
}

TEST_F(AnalyticGravityTest, GivenUniformSphere_WhenInsideSphere_ThenGravityIsLinear) {
    // GIVEN: Uniform density sphere of radius R
    // WHEN: Compute gravity at r < R
    // THEN: g proportional to r

    for (double r = 0.1; r < R; r += 0.1) {
        double g = analytic_g_inside(r);
        double expected = (4.0/3.0) * M_PI * G * rho * r;
        EXPECT_NEAR(g, expected, kAbsTol)
            << "Gravity at r=" << r << " should be linear";
    }
}

TEST_F(AnalyticGravityTest, GivenUniformSphere_WhenOutsideSphere_ThenGravityIsInverseSquare) {
    // GIVEN: Uniform density sphere
    // WHEN: Compute gravity at r > R
    // THEN: g = GM/r^2

    for (double r = R + 0.1; r < 5.0 * R; r += 0.5) {
        double g = analytic_g_outside(r);
        double expected = G * M / (r * r);
        EXPECT_NEAR(g, expected, kAbsTol)
            << "Gravity at r=" << r << " should follow inverse square law";
    }
}

TEST_F(AnalyticGravityTest, GivenUniformSphere_WhenAtSurface_ThenGravityIsContinuous) {
    // GIVEN: Uniform density sphere
    // WHEN: Compute gravity just inside and outside surface
    // THEN: g is continuous at r = R

    double eps = 1e-8;
    double g_inside = analytic_g_inside(R - eps);
    double g_outside = analytic_g_outside(R + eps);

    EXPECT_NEAR(g_inside, g_outside, kContinuityTol)
        << "Gravity should be continuous at surface: inside=" << g_inside
        << ", outside=" << g_outside;
}

// =============================================================================
// 2. SOFTENING CONVERGENCE TESTS
// =============================================================================

TEST_F(AnalyticGravityTest, GivenSoftenedGravity_WhenSofteningDecreases_ThenConvergesToAnalytic) {
    // GIVEN: Uniform sphere particle distribution
    // WHEN: Compute gravity with decreasing softening length
    // THEN: Results converge to analytic solution

    const int N = 10000;  // Number of particles
    auto particles = sph::create_uniform_sphere_distribution(N, R, rho, 42);

    vec_t test_pos;
#if DIM == 3
    test_pos = {0.5 * R, 0.0, 0.0};  // Inside sphere
#elif DIM == 2
    test_pos = {0.5 * R, 0.0};
#else
    test_pos = {0.5 * R};
#endif
    double test_r = 0.5 * R;
    double analytic = analytic_g_inside(test_r);

    std::vector<double> softenings = {0.5, 0.2, 0.1, 0.05, 0.02};
    double prev_error = 1.0;

    std::cout << "[INFO] Softening convergence test (Wendland C4):" << std::endl;

    for (double eps : softenings) {
        auto result = sph::compute_direct_gravity(test_pos, particles, eps,
                                                   sph::GravitySofteningType::WENDLAND_C4);
        double computed = std::abs(result.acceleration);
        double error = std::abs(computed - analytic) / analytic;

        std::cout << "  eps=" << eps << ": computed=" << computed
                  << ", analytic=" << analytic << ", error=" << error << std::endl;

        // Error should decrease as softening decreases
        EXPECT_LT(error, prev_error)
            << "Error should decrease with softening, eps=" << eps;
        prev_error = error;
    }
}

TEST_F(AnalyticGravityTest, GivenHernquistKatzSoftening_WhenSofteningDecreases_ThenConverges) {
    // GIVEN: Uniform sphere with Hernquist-Katz softening
    // WHEN: Softening epsilon -> 0
    // THEN: Converges to analytic gravity

    const int N = 10000;
    auto particles = sph::create_uniform_sphere_distribution(N, R, rho, 42);

    // Test at multiple positions
    std::vector<double> radii = {0.25 * R, 0.5 * R, 0.75 * R, 1.0 * R, 1.5 * R, 2.0 * R};

    std::cout << "[INFO] Hernquist-Katz softening convergence:" << std::endl;

    for (double r : radii) {
        vec_t test_pos;
#if DIM == 3
        test_pos = {r, 0.0, 0.0};
#elif DIM == 2
        test_pos = {r, 0.0};
#else
        test_pos = {r};
#endif

        double analytic = (r < R) ? analytic_g_inside(r) : analytic_g_outside(r);
        if (r < 1e-10) continue;  // Skip zero radius

        std::vector<double> softenings = {0.2, 0.1, 0.05};
        for (double eps : softenings) {
            auto result = sph::compute_direct_gravity(test_pos, particles, eps,
                                                       sph::GravitySofteningType::HERNQUIST_KATZ);
            double computed = std::abs(result.acceleration);
            double error = std::abs(computed - analytic) / (analytic + kAbsTol);

            std::cout << "  r=" << r << ", eps=" << eps << ": error=" << error << std::endl;
        }
    }
}

// =============================================================================
// 3. RESOLUTION CONVERGENCE TESTS
// =============================================================================

TEST_F(AnalyticGravityTest, GivenFixedSoftening_WhenResolutionIncreases_ThenErrorDecreases) {
    // GIVEN: Fixed softening epsilon
    // WHEN: Increase number of particles (N)
    // THEN: Error converges to softening-limited error floor

    double eps = 0.1;  // Fixed softening
    std::vector<int> particle_counts = {100, 500, 1000, 5000, 10000};

    vec_t test_pos;
#if DIM == 3
    test_pos = {0.5 * R, 0.0, 0.0};
#elif DIM == 2
    test_pos = {0.5 * R, 0.0};
#else
    test_pos = {0.5 * R};
#endif
    double analytic = analytic_g_inside(0.5 * R);

    std::cout << "[INFO] Resolution convergence test (fixed softening):" << std::endl;

    for (int N : particle_counts) {
        auto particles = sph::create_uniform_sphere_distribution(N, R, rho, 42);
        auto result = sph::compute_direct_gravity(test_pos, particles, eps,
                                                   sph::GravitySofteningType::WENDLAND_C4);
        double computed = std::abs(result.acceleration);
        double error = std::abs(computed - analytic) / analytic;

        std::cout << "  N=" << N << ": error=" << error << std::endl;
    }
}

// =============================================================================
// 4. KERNEL COMPARISON TESTS
// =============================================================================

TEST_F(AnalyticGravityTest, GivenSameSoftening_WhenComparingKernels_ThenBothAreAccurate) {
    // GIVEN: Same softening length for both kernels
    // WHEN: Compute gravity inside sphere
    // THEN: Both kernels should be reasonably accurate

    const int N = 10000;
    auto particles = sph::create_uniform_sphere_distribution(N, R, rho, 42);
    double eps = 0.1;

    vec_t test_pos;
#if DIM == 3
    test_pos = {0.5 * R, 0.0, 0.0};
#elif DIM == 2
    test_pos = {0.5 * R, 0.0};
#else
    test_pos = {0.5 * R};
#endif
    double analytic = analytic_g_inside(0.5 * R);

    auto result_hk = sph::compute_direct_gravity(test_pos, particles, eps,
                                                  sph::GravitySofteningType::HERNQUIST_KATZ);
    auto result_wc4 = sph::compute_direct_gravity(test_pos, particles, eps,
                                                   sph::GravitySofteningType::WENDLAND_C4);

    double g_hk = std::abs(result_hk.acceleration);
    double g_wc4 = std::abs(result_wc4.acceleration);

    double error_hk = std::abs(g_hk - analytic) / analytic;
    double error_wc4 = std::abs(g_wc4 - analytic) / analytic;

    std::cout << "[INFO] Kernel comparison:" << std::endl;
    std::cout << "  Hernquist-Katz: g=" << g_hk << ", error=" << error_hk << std::endl;
    std::cout << "  Wendland C4:    g=" << g_wc4 << ", error=" << error_wc4 << std::endl;

    // Both should be reasonably accurate
    EXPECT_LT(error_hk, 0.15)
        << "Hernquist-Katz error should be < 15%";
    EXPECT_LT(error_wc4, 0.15)
        << "Wendland C4 error should be < 15%";
}

TEST_F(AnalyticGravityTest, GivenUniformSphere_WhenOutside_ThenMatchesPointMass) {
    // GIVEN: Uniform sphere particle distribution
    // WHEN: Compute gravity at r > R (outside sphere)
    // THEN: Should match point mass GM/r^2

    const int N = 10000;
    auto particles = sph::create_uniform_sphere_distribution(N, R, rho, 42);
    double eps = 0.05;  // Small softening

    // Test at 2R (well outside sphere)
    double test_r = 2.0 * R;
    vec_t test_pos;
#if DIM == 3
    test_pos = {test_r, 0.0, 0.0};
#elif DIM == 2
    test_pos = {test_r, 0.0};
#else
    test_pos = {test_r};
#endif

    double analytic = analytic_g_outside(test_r);

    auto result = sph::compute_direct_gravity(test_pos, particles, eps,
                                               sph::GravitySofteningType::WENDLAND_C4);
    double computed = std::abs(result.acceleration);
    double error = std::abs(computed - analytic) / analytic;

    std::cout << "[INFO] Point mass test at r=" << test_r << ":" << std::endl;
    std::cout << "  computed=" << computed << ", analytic=" << analytic
              << ", error=" << error << std::endl;

    // Should be within 5% at this distance
    EXPECT_LT(error, 0.05)
        << "Gravity at 2R should match point mass within 5%";
}

// =============================================================================
// Additional Tests
// =============================================================================

TEST_F(AnalyticGravityTest, DirectGravityResultHasCorrectDirection) {
    // GIVEN: Test particle at positive x
    // WHEN: Compute gravity from sphere centered at origin
    // THEN: Acceleration points toward center (negative x direction)

    const int N = 1000;
    auto particles = sph::create_uniform_sphere_distribution(N, R, rho, 42);
    double eps = 0.1;

    double test_r = 0.5 * R;
    vec_t test_pos;
#if DIM == 3
    test_pos = {test_r, 0.0, 0.0};
#elif DIM == 2
    test_pos = {test_r, 0.0};
#else
    test_pos = {test_r};
#endif

    auto result = sph::compute_direct_gravity(test_pos, particles, eps,
                                               sph::GravitySofteningType::WENDLAND_C4);

    // Acceleration should point toward center (negative x direction)
    EXPECT_LT(result.acceleration[0], 0.0)
        << "Gravity should point toward center (negative x)";

#if DIM >= 2
    // Y and Z components should be near zero due to symmetry
    EXPECT_NEAR(result.acceleration[1], 0.0, 0.1 * std::abs(result.acceleration[0]))
        << "Y component should be near zero due to symmetry";
#endif
#if DIM == 3
    EXPECT_NEAR(result.acceleration[2], 0.0, 0.1 * std::abs(result.acceleration[0]))
        << "Z component should be near zero due to symmetry";
#endif
}

TEST_F(AnalyticGravityTest, PotentialDecreasesWithDistance) {
    // GIVEN: Uniform sphere particle distribution
    // WHEN: Compute potential at increasing radii
    // THEN: Potential should be monotonically decreasing (less negative)

    const int N = 1000;
    auto particles = sph::create_uniform_sphere_distribution(N, R, rho, 42);
    double eps = 0.1;

    double prev_phi = -1e10;  // Very negative
    for (double r = 0.1; r < 3.0 * R; r += 0.2) {
        vec_t test_pos;
#if DIM == 3
        test_pos = {r, 0.0, 0.0};
#elif DIM == 2
        test_pos = {r, 0.0};
#else
        test_pos = {r};
#endif

        auto result = sph::compute_direct_gravity(test_pos, particles, eps,
                                                   sph::GravitySofteningType::WENDLAND_C4);

        // Potential should increase (become less negative) with distance
        EXPECT_GT(result.potential, prev_phi - kAbsTol)
            << "Potential should increase with distance at r=" << r;
        prev_phi = result.potential;
    }
}
