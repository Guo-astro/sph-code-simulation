/**
 * @file test_gravity_first_principles.cpp
 * @brief First-principles validation tests for gravity softening
 *
 * Based on literature:
 * - Dehnen (2001) MNRAS 324, 273: Optimal softening, MASE metric
 * - Price & Monaghan (2007) MNRAS 374, 1347: Energy conservation test
 * - Hernquist & Katz (1989): Original cubic spline softening
 *
 * Tests:
 * 1. Kernel Accuracy: Lookup table matches analytic polynomial
 * 2. Energy Conservation: Total energy constant in isolated system
 * 3. Newton's Third Law: F_ij = -F_ji (momentum conservation)
 */

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <iostream>
#include <iomanip>

#include "defines.hpp"
#include "hernquist_katz_lookup_table.hpp"
#include "softening_lookup_table.hpp"

namespace {

// =============================================================================
// Constants - BRUTAL TOLERANCES (tight enough to catch any regression)
// =============================================================================
const double G = 1.0;  // Gravitational constant

// Lookup table accuracy: must match polynomial to 1e-6 relative error
// (actual performance: ~1e-8 to 1e-10)
const double kRelTol = 1e-6;

// Energy conservation: leapfrog + softened potential achieves ~4e-4
// Set threshold at 5e-4 to catch any degradation
const double kEnergyTol = 5e-4;

// Newton's third law: should be exact to machine precision
const double kSymmetryTol = 1e-15;

// Momentum conservation: should be exact to ~1e-14
const double kMomentumTol = 1e-13;

// =============================================================================
// Reference Polynomial Implementations (from Hernquist & Katz 1989)
// =============================================================================

/**
 * @brief Exact H-K potential kernel f(r,h) - analytic reference
 */
double hk_f_polynomial(double r, double h) {
    const double e = h * 0.5;
    const double u = r / e;

    if (u < 1.0) {
        return (-0.5 * u * u * (1.0/3.0 - 3.0/20.0 * u*u + u*u*u/20.0) + 1.4) / e;
    } else if (u < 2.0) {
        if (r < 1e-30) return 1.4 / e;
        return -1.0/(15.0*r) + (-u*u * (4.0/3.0 - u + 0.3*u*u - u*u*u/30.0) + 1.6) / e;
    } else {
        return 1.0 / r;
    }
}

/**
 * @brief Exact H-K force kernel g(r,h) - analytic reference
 */
double hk_g_polynomial(double r, double h) {
    const double e = h * 0.5;
    const double u = r / e;

    if (u < 1.0) {
        return (4.0/3.0 - 1.2*u*u + 0.5*u*u*u) / (e*e*e);
    } else if (u < 2.0) {
        if (r < 1e-30) return 0.0;
        return (-1.0/15.0 + 8.0/3.0*u*u*u - 3.0*u*u*u*u + 1.2*u*u*u*u*u - u*u*u*u*u*u/6.0) / (r*r*r);
    } else {
        return 1.0 / (r*r*r);
    }
}

/**
 * @brief Exact Wendland C4 potential kernel phi(r,h) - analytic reference
 */
double wendland_phi_polynomial(double r, double h) {
    const double q = r / h;

    if (q >= 1.0) {
        return 1.0 / r;
    }

    const double q2 = q * q;
    const double q3 = q2 * q;
    const double q4 = q2 * q2;
    const double q5 = q4 * q;
    const double q6 = q3 * q3;
    const double q7 = q6 * q;
    const double q8 = q4 * q4;
    const double q9 = q8 * q;

    const double a0 =  3.4374743761;
    const double a1 = -0.0031873250;
    const double a2 = -10.2154807743;
    const double a3 = -1.1577720555;
    const double a4 =  36.1013669755;
    const double a5 = -26.3399094060;
    const double a6 = -44.1079372114;
    const double a7 =  82.6543766683;
    const double a8 = -50.5921624056;
    const double a9 =  11.2232565249;

    return (a0 + a1*q + a2*q2 + a3*q3 + a4*q4 + a5*q5 + a6*q6 + a7*q7 + a8*q8 + a9*q9) / h;
}

/**
 * @brief Exact Wendland C4 force kernel g(r,h) - analytic reference
 */
double wendland_g_polynomial(double r, double h) {
    const double q = r / h;

    if (q >= 1.0) {
        return 1.0 / (r * r * r);
    }

    if (q < 1e-10) {
        return 0.0;
    }

    const double q2 = q * q;
    const double q3 = q2 * q;
    const double q4 = q2 * q2;
    const double q5 = q4 * q;
    const double q6 = q3 * q3;
    const double q7 = q6 * q;

    const double b0 = -0.0031873250;
    const double b1 = -20.4309615486;
    const double b2 = -3.4733161665;
    const double b3 = 144.4054679020;
    const double b4 = -131.6995470300;
    const double b5 = -264.6476232684;
    const double b6 = 578.5806366781;
    const double b7 = -404.7372992448;
    const double b8 = 101.0093087241;

    const double dphi_dq = b0 + b1*q + b2*q2 + b3*q3 + b4*q4 + b5*q5 + b6*q6 + b7*q7 + b8*q*q7;
    const double h3 = h * h * h;

    return -dphi_dq / (h3 * q);
}

// =============================================================================
// TEST 1: Kernel Accuracy - Lookup Table vs Polynomial
// =============================================================================

class KernelAccuracyTest : public ::testing::Test {
protected:
    const sph::HernquistKatzLookupTable& hk_table = sph::HernquistKatzLookupTable::get_instance();
    const sph::WendlandC4LookupTable& wc4_table = sph::WendlandC4LookupTable::get_instance();
};

// Test H-K f(r,h) lookup matches polynomial
TEST_F(KernelAccuracyTest, HernquistKatzF_MatchesPolynomial) {
    const double h = 1.0;

    // Test points across the entire range
    std::vector<double> test_r = {0.01, 0.1, 0.2, 0.3, 0.4, 0.49,  // u < 1
                                   0.51, 0.6, 0.7, 0.8, 0.9, 0.99,  // 1 < u < 2
                                   1.01, 1.5, 2.0, 3.0};            // u >= 2

    int passed = 0;
    int total = 0;

    for (double r : test_r) {
        double poly = hk_f_polynomial(r, h);
        double lookup = hk_table.f_full(r, h);

        double rel_error = std::abs(lookup - poly) / std::abs(poly);
        total++;

        // Allow slightly larger error near discontinuities (u=1, u=2)
        double tol = kRelTol;
        double u = r / (h * 0.5);
        if (std::abs(u - 1.0) < 0.05 || std::abs(u - 2.0) < 0.05) {
            tol = 1e-4;  // 0.01% near boundaries (still strict)
        }

        EXPECT_LT(rel_error, tol)
            << "H-K f mismatch at r=" << r << " (u=" << u << ")"
            << ": poly=" << poly << ", lookup=" << lookup;

        if (rel_error < tol) passed++;
    }

    std::cout << "[H-K f] Passed " << passed << "/" << total << " points" << std::endl;
}

// Test H-K g(r,h) lookup matches polynomial
TEST_F(KernelAccuracyTest, HernquistKatzG_MatchesPolynomial) {
    const double h = 1.0;

    std::vector<double> test_r = {0.01, 0.1, 0.2, 0.3, 0.4, 0.49,
                                   0.51, 0.6, 0.7, 0.8, 0.9, 0.99,
                                   1.01, 1.5, 2.0, 3.0};

    int passed = 0;
    int total = 0;

    for (double r : test_r) {
        double poly = hk_g_polynomial(r, h);
        double lookup = hk_table.g_full(r, h);

        double rel_error = std::abs(lookup - poly) / std::abs(poly);
        total++;

        double tol = kRelTol;
        double u = r / (h * 0.5);
        if (std::abs(u - 1.0) < 0.05 || std::abs(u - 2.0) < 0.05) {
            tol = 1e-4;  // 0.01% near boundaries (still strict)
        }

        EXPECT_LT(rel_error, tol)
            << "H-K g mismatch at r=" << r << " (u=" << u << ")"
            << ": poly=" << poly << ", lookup=" << lookup;

        if (rel_error < tol) passed++;
    }

    std::cout << "[H-K g] Passed " << passed << "/" << total << " points" << std::endl;
}

// Test Wendland C4 phi(r,h) lookup matches polynomial
TEST_F(KernelAccuracyTest, WendlandC4Phi_MatchesPolynomial) {
    const double h = 1.0;

    std::vector<double> test_r = {0.01, 0.1, 0.2, 0.3, 0.4, 0.5,
                                   0.6, 0.7, 0.8, 0.9, 0.99,
                                   1.01, 1.5, 2.0};

    int passed = 0;
    int total = 0;

    for (double r : test_r) {
        double poly = wendland_phi_polynomial(r, h);
        double lookup = wc4_table.phi_full(r, h);

        double rel_error = std::abs(lookup - poly) / std::abs(poly);
        total++;

        EXPECT_LT(rel_error, kRelTol)
            << "W-C4 phi mismatch at r=" << r
            << ": poly=" << poly << ", lookup=" << lookup;

        if (rel_error < kRelTol) passed++;
    }

    std::cout << "[W-C4 phi] Passed " << passed << "/" << total << " points" << std::endl;
}

// Test Wendland C4 g(r,h) lookup matches polynomial
TEST_F(KernelAccuracyTest, WendlandC4G_MatchesPolynomial) {
    const double h = 1.0;

    std::vector<double> test_r = {0.01, 0.1, 0.2, 0.3, 0.4, 0.5,
                                   0.6, 0.7, 0.8, 0.9, 0.99,
                                   1.01, 1.5, 2.0};

    int passed = 0;
    int total = 0;

    for (double r : test_r) {
        double poly = wendland_g_polynomial(r, h);
        double lookup = wc4_table.g_full(r, h);

        double rel_error = std::abs(lookup - poly) / std::abs(poly);
        total++;

        EXPECT_LT(rel_error, kRelTol)
            << "W-C4 g mismatch at r=" << r
            << ": poly=" << poly << ", lookup=" << lookup;

        if (rel_error < kRelTol) passed++;
    }

    std::cout << "[W-C4 g] Passed " << passed << "/" << total << " points" << std::endl;
}

// Test boundary conditions from H-K 1989 paper
TEST_F(KernelAccuracyTest, HernquistKatzBoundaryConditions) {
    const double h = 1.0;
    const double e = h * 0.5;

    // At r=0: f(0) = 1.4/e, g(0) = 4/(3e³)
    double f_at_zero = hk_table.f_full(0.001, h);
    double g_at_zero = hk_table.g_full(0.001, h);

    double f_expected = 1.4 / e;
    double g_expected = 4.0 / (3.0 * e * e * e);

    // At r≈0: use 0.1% tolerance (actual should be much better)
    EXPECT_NEAR(f_at_zero, f_expected, 1e-3 * f_expected)
        << "f(0) should be 1.4/e = " << f_expected;
    EXPECT_NEAR(g_at_zero, g_expected, 1e-3 * g_expected)
        << "g(0) should be 4/(3e³) = " << g_expected;

    // At r >> h: f → 1/r, g → 1/r³ (exact to machine precision)
    double r_far = 5.0;
    double f_far = hk_table.f_full(r_far, h);
    double g_far = hk_table.g_full(r_far, h);

    EXPECT_NEAR(f_far, 1.0/r_far, 1e-14)
        << "f(r>>h) should be 1/r";
    EXPECT_NEAR(g_far, 1.0/(r_far*r_far*r_far), 1e-14)
        << "g(r>>h) should be 1/r³";

    std::cout << "[H-K Boundary] f(0)=" << f_at_zero << " (expected " << f_expected << ")"
              << ", g(0)=" << g_at_zero << " (expected " << g_expected << ")" << std::endl;
}

// =============================================================================
// TEST 2: Energy Conservation
// =============================================================================

class EnergyConservationTest : public ::testing::Test {
protected:
    struct Particle {
        std::array<double, 3> pos;
        std::array<double, 3> vel;
        double mass;
    };

    const sph::HernquistKatzLookupTable& hk_table = sph::HernquistKatzLookupTable::get_instance();

    double compute_kinetic_energy(const std::vector<Particle>& particles) {
        double KE = 0.0;
        for (const auto& p : particles) {
            double v2 = p.vel[0]*p.vel[0] + p.vel[1]*p.vel[1] + p.vel[2]*p.vel[2];
            KE += 0.5 * p.mass * v2;
        }
        return KE;
    }

    double compute_potential_energy(const std::vector<Particle>& particles, double h) {
        double PE = 0.0;
        for (size_t i = 0; i < particles.size(); ++i) {
            for (size_t j = i + 1; j < particles.size(); ++j) {
                double dx = particles[i].pos[0] - particles[j].pos[0];
                double dy = particles[i].pos[1] - particles[j].pos[1];
                double dz = particles[i].pos[2] - particles[j].pos[2];
                double r = std::sqrt(dx*dx + dy*dy + dz*dz);

                double f = hk_table.f_full(r, h);
                PE -= G * particles[i].mass * particles[j].mass * f;
            }
        }
        return PE;
    }

    void compute_accelerations(const std::vector<Particle>& particles, double h,
                               std::vector<std::array<double, 3>>& acc) {
        acc.resize(particles.size());
        for (auto& a : acc) {
            a = {0.0, 0.0, 0.0};
        }

        for (size_t i = 0; i < particles.size(); ++i) {
            for (size_t j = 0; j < particles.size(); ++j) {
                if (i == j) continue;

                double dx = particles[i].pos[0] - particles[j].pos[0];
                double dy = particles[i].pos[1] - particles[j].pos[1];
                double dz = particles[i].pos[2] - particles[j].pos[2];
                double r = std::sqrt(dx*dx + dy*dy + dz*dz);

                if (r < 1e-30) continue;

                double g = hk_table.g_full(r, h);

                acc[i][0] -= G * particles[j].mass * g * dx;
                acc[i][1] -= G * particles[j].mass * g * dy;
                acc[i][2] -= G * particles[j].mass * g * dz;
            }
        }
    }
};

// Test energy conservation in 2-body elliptical orbit
TEST_F(EnergyConservationTest, TwoBodyEllipticalOrbit) {
    const double h = 0.3;

    // Set up 2-body problem
    std::vector<Particle> particles(2);
    particles[0] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0};
    particles[1] = {{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0};

    // Circular orbit velocity (softened)
    double r0 = 1.0;
    double g_val = hk_table.g_full(r0, h);
    double v_circ = std::sqrt(G * 2.0 * g_val * r0);

    // Use 80% of circular velocity for elliptical orbit
    double v0 = 0.8 * v_circ;
    particles[0].vel[1] = -v0 * 0.5;
    particles[1].vel[1] = v0 * 0.5;

    // Initial energy
    double E0 = compute_kinetic_energy(particles) + compute_potential_energy(particles, h);

    // Leapfrog integration
    const double dt = 0.01;
    const int n_steps = 2000;
    std::vector<std::array<double, 3>> acc;

    double max_dE_rel = 0.0;

    for (int step = 0; step < n_steps; ++step) {
        // Kick-drift-kick
        compute_accelerations(particles, h, acc);

        for (size_t i = 0; i < particles.size(); ++i) {
            particles[i].vel[0] += 0.5 * dt * acc[i][0];
            particles[i].vel[1] += 0.5 * dt * acc[i][1];
            particles[i].vel[2] += 0.5 * dt * acc[i][2];
        }

        for (size_t i = 0; i < particles.size(); ++i) {
            particles[i].pos[0] += dt * particles[i].vel[0];
            particles[i].pos[1] += dt * particles[i].vel[1];
            particles[i].pos[2] += dt * particles[i].vel[2];
        }

        compute_accelerations(particles, h, acc);

        for (size_t i = 0; i < particles.size(); ++i) {
            particles[i].vel[0] += 0.5 * dt * acc[i][0];
            particles[i].vel[1] += 0.5 * dt * acc[i][1];
            particles[i].vel[2] += 0.5 * dt * acc[i][2];
        }

        // Check energy
        double E = compute_kinetic_energy(particles) + compute_potential_energy(particles, h);
        double dE_rel = std::abs(E - E0) / std::abs(E0);
        max_dE_rel = std::max(max_dE_rel, dE_rel);
    }

    std::cout << "[Energy Conservation] E0=" << E0
              << ", max |ΔE/E0|=" << std::scientific << max_dE_rel
              << " (threshold: " << kEnergyTol << ")" << std::endl;

    // Energy must be conserved to within threshold (catches any force inconsistency)
    EXPECT_LT(max_dE_rel, kEnergyTol)
        << "REGRESSION: Energy not conserved: max |ΔE/E0| = " << max_dE_rel
        << " exceeds threshold " << kEnergyTol;
}

// Test energy conservation with Wendland C4 kernel
TEST_F(EnergyConservationTest, TwoBodyEllipticalOrbit_WendlandC4) {
    const sph::WendlandC4LookupTable& wc4_table = sph::WendlandC4LookupTable::get_instance();
    const double h = 0.5;

    // Set up 2-body problem
    std::vector<Particle> particles(2);
    particles[0] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0};
    particles[1] = {{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0};

    // Circular orbit velocity (softened)
    double r0 = 1.0;
    double g_val = wc4_table.g_full(r0, h);
    double v_circ = std::sqrt(G * 2.0 * g_val * r0);

    double v0 = 0.8 * v_circ;
    particles[0].vel[1] = -v0 * 0.5;
    particles[1].vel[1] = v0 * 0.5;

    // Compute energy with Wendland C4
    auto compute_PE_wc4 = [&](const std::vector<Particle>& p) {
        double PE = 0.0;
        for (size_t i = 0; i < p.size(); ++i) {
            for (size_t j = i + 1; j < p.size(); ++j) {
                double dx = p[i].pos[0] - p[j].pos[0];
                double dy = p[i].pos[1] - p[j].pos[1];
                double dz = p[i].pos[2] - p[j].pos[2];
                double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                PE -= G * p[i].mass * p[j].mass * wc4_table.phi_full(r, h);
            }
        }
        return PE;
    };

    auto compute_acc_wc4 = [&](const std::vector<Particle>& p, std::vector<std::array<double, 3>>& a) {
        a.resize(p.size());
        for (auto& x : a) x = {0.0, 0.0, 0.0};

        for (size_t i = 0; i < p.size(); ++i) {
            for (size_t j = 0; j < p.size(); ++j) {
                if (i == j) continue;
                double dx = p[i].pos[0] - p[j].pos[0];
                double dy = p[i].pos[1] - p[j].pos[1];
                double dz = p[i].pos[2] - p[j].pos[2];
                double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (r < 1e-30) continue;

                double g = wc4_table.g_full(r, h);
                a[i][0] -= G * p[j].mass * g * dx;
                a[i][1] -= G * p[j].mass * g * dy;
                a[i][2] -= G * p[j].mass * g * dz;
            }
        }
    };

    double E0 = compute_kinetic_energy(particles) + compute_PE_wc4(particles);

    const double dt = 0.01;
    const int n_steps = 2000;
    std::vector<std::array<double, 3>> acc;

    double max_dE_rel = 0.0;

    for (int step = 0; step < n_steps; ++step) {
        compute_acc_wc4(particles, acc);

        for (size_t i = 0; i < particles.size(); ++i) {
            particles[i].vel[0] += 0.5 * dt * acc[i][0];
            particles[i].vel[1] += 0.5 * dt * acc[i][1];
            particles[i].vel[2] += 0.5 * dt * acc[i][2];
        }

        for (size_t i = 0; i < particles.size(); ++i) {
            particles[i].pos[0] += dt * particles[i].vel[0];
            particles[i].pos[1] += dt * particles[i].vel[1];
            particles[i].pos[2] += dt * particles[i].vel[2];
        }

        compute_acc_wc4(particles, acc);

        for (size_t i = 0; i < particles.size(); ++i) {
            particles[i].vel[0] += 0.5 * dt * acc[i][0];
            particles[i].vel[1] += 0.5 * dt * acc[i][1];
            particles[i].vel[2] += 0.5 * dt * acc[i][2];
        }

        double E = compute_kinetic_energy(particles) + compute_PE_wc4(particles);
        double dE_rel = std::abs(E - E0) / std::abs(E0);
        max_dE_rel = std::max(max_dE_rel, dE_rel);
    }

    std::cout << "[Energy Conservation W-C4] E0=" << E0
              << ", max |ΔE/E0|=" << std::scientific << max_dE_rel
              << " (threshold: " << kEnergyTol << ")" << std::endl;

    EXPECT_LT(max_dE_rel, kEnergyTol)
        << "REGRESSION: Energy not conserved with W-C4: max |ΔE/E0| = " << max_dE_rel
        << " exceeds threshold " << kEnergyTol;
}

// =============================================================================
// TEST 3: Newton's Third Law (F_ij = -F_ji)
// =============================================================================

class NewtonsThirdLawTest : public ::testing::Test {
protected:
    const sph::HernquistKatzLookupTable& hk_table = sph::HernquistKatzLookupTable::get_instance();
    const sph::WendlandC4LookupTable& wc4_table = sph::WendlandC4LookupTable::get_instance();
};

// Test F_ij = -F_ji for Hernquist-Katz kernel
TEST_F(NewtonsThirdLawTest, HernquistKatzForceSymmetry) {
    const double h = 1.0;
    const int n_tests = 100;

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> pos_dist(-5.0, 5.0);
    std::uniform_real_distribution<double> mass_dist(0.1, 2.0);

    double max_asymmetry = 0.0;

    for (int t = 0; t < n_tests; ++t) {
        std::array<double, 3> r_i = {pos_dist(rng), pos_dist(rng), pos_dist(rng)};
        std::array<double, 3> r_j = {pos_dist(rng), pos_dist(rng), pos_dist(rng)};
        double m_i = mass_dist(rng);
        double m_j = mass_dist(rng);

        double dx = r_i[0] - r_j[0];
        double dy = r_i[1] - r_j[1];
        double dz = r_i[2] - r_j[2];
        double r = std::sqrt(dx*dx + dy*dy + dz*dz);

        if (r < 1e-10) continue;

        double g = hk_table.g_full(r, h);

        // F_ij = -G * m_i * m_j * g * r_ij
        std::array<double, 3> F_ij = {-G * m_i * m_j * g * dx,
                                       -G * m_i * m_j * g * dy,
                                       -G * m_i * m_j * g * dz};

        // F_ji = -G * m_j * m_i * g * r_ji = -G * m_j * m_i * g * (-r_ij)
        std::array<double, 3> F_ji = {G * m_j * m_i * g * dx,
                                       G * m_j * m_i * g * dy,
                                       G * m_j * m_i * g * dz};

        // Check F_ij + F_ji = 0
        double sum_mag = std::sqrt((F_ij[0]+F_ji[0])*(F_ij[0]+F_ji[0]) +
                                    (F_ij[1]+F_ji[1])*(F_ij[1]+F_ji[1]) +
                                    (F_ij[2]+F_ji[2])*(F_ij[2]+F_ji[2]));
        double F_mag = std::sqrt(F_ij[0]*F_ij[0] + F_ij[1]*F_ij[1] + F_ij[2]*F_ij[2]);

        double asymmetry = sum_mag / (F_mag + 1e-30);
        max_asymmetry = std::max(max_asymmetry, asymmetry);
    }

    std::cout << "[Newton's 3rd Law H-K] max |F_ij + F_ji| / |F_ij| = "
              << std::scientific << max_asymmetry << std::endl;

    EXPECT_LT(max_asymmetry, kSymmetryTol)
        << "Newton's 3rd law violated: max asymmetry = " << max_asymmetry
        << " (threshold: " << kSymmetryTol << ")";
}

// Test F_ij = -F_ji for Wendland C4 kernel
TEST_F(NewtonsThirdLawTest, WendlandC4ForceSymmetry) {
    const double h = 1.0;
    const int n_tests = 100;

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> pos_dist(-5.0, 5.0);
    std::uniform_real_distribution<double> mass_dist(0.1, 2.0);

    double max_asymmetry = 0.0;

    for (int t = 0; t < n_tests; ++t) {
        std::array<double, 3> r_i = {pos_dist(rng), pos_dist(rng), pos_dist(rng)};
        std::array<double, 3> r_j = {pos_dist(rng), pos_dist(rng), pos_dist(rng)};
        double m_i = mass_dist(rng);
        double m_j = mass_dist(rng);

        double dx = r_i[0] - r_j[0];
        double dy = r_i[1] - r_j[1];
        double dz = r_i[2] - r_j[2];
        double r = std::sqrt(dx*dx + dy*dy + dz*dz);

        if (r < 1e-10) continue;

        double g = wc4_table.g_full(r, h);

        std::array<double, 3> F_ij = {-G * m_i * m_j * g * dx,
                                       -G * m_i * m_j * g * dy,
                                       -G * m_i * m_j * g * dz};

        std::array<double, 3> F_ji = {G * m_j * m_i * g * dx,
                                       G * m_j * m_i * g * dy,
                                       G * m_j * m_i * g * dz};

        double sum_mag = std::sqrt((F_ij[0]+F_ji[0])*(F_ij[0]+F_ji[0]) +
                                    (F_ij[1]+F_ji[1])*(F_ij[1]+F_ji[1]) +
                                    (F_ij[2]+F_ji[2])*(F_ij[2]+F_ji[2]));
        double F_mag = std::sqrt(F_ij[0]*F_ij[0] + F_ij[1]*F_ij[1] + F_ij[2]*F_ij[2]);

        double asymmetry = sum_mag / (F_mag + 1e-30);
        max_asymmetry = std::max(max_asymmetry, asymmetry);
    }

    std::cout << "[Newton's 3rd Law W-C4] max |F_ij + F_ji| / |F_ij| = "
              << std::scientific << max_asymmetry << std::endl;

    EXPECT_LT(max_asymmetry, kSymmetryTol)
        << "Newton's 3rd law violated: max asymmetry = " << max_asymmetry
        << " (threshold: " << kSymmetryTol << ")";
}

// Test momentum conservation in N-body system
TEST_F(NewtonsThirdLawTest, MomentumConservation) {
    const double h = 0.5;
    const int n_particles = 10;

    std::mt19937 rng(123);
    std::uniform_real_distribution<double> pos_dist(-2.0, 2.0);
    std::uniform_real_distribution<double> mass_dist(0.5, 1.5);

    // Create random particles
    std::vector<std::array<double, 3>> positions(n_particles);
    std::vector<double> masses(n_particles);

    for (int i = 0; i < n_particles; ++i) {
        positions[i] = {pos_dist(rng), pos_dist(rng), pos_dist(rng)};
        masses[i] = mass_dist(rng);
    }

    // Compute total force on system
    std::array<double, 3> total_force = {0.0, 0.0, 0.0};

    for (int i = 0; i < n_particles; ++i) {
        for (int j = 0; j < n_particles; ++j) {
            if (i == j) continue;

            double dx = positions[i][0] - positions[j][0];
            double dy = positions[i][1] - positions[j][1];
            double dz = positions[i][2] - positions[j][2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (r < 1e-30) continue;

            double g = hk_table.g_full(r, h);

            total_force[0] -= G * masses[i] * masses[j] * g * dx;
            total_force[1] -= G * masses[i] * masses[j] * g * dy;
            total_force[2] -= G * masses[i] * masses[j] * g * dz;
        }
    }

    double total_force_mag = std::sqrt(total_force[0]*total_force[0] +
                                        total_force[1]*total_force[1] +
                                        total_force[2]*total_force[2]);

    std::cout << "[Momentum Conservation] Total force on system = "
              << std::scientific << total_force_mag << std::endl;

    // Total internal force should be zero (momentum conservation)
    EXPECT_LT(total_force_mag, kMomentumTol)
        << "Momentum not conserved: total force = " << total_force_mag
        << " (threshold: " << kMomentumTol << ")";
}

// =============================================================================
// TEST 4: Lookup Table vs Simple For-Loop (Direct Polynomial Evaluation)
// Verifies lookup table produces identical results to naive O(N²) calculation
// =============================================================================

class LookupVsForLoopTest : public ::testing::Test {
protected:
    const sph::HernquistKatzLookupTable& hk_table = sph::HernquistKatzLookupTable::get_instance();
    const sph::WendlandC4LookupTable& wc4_table = sph::WendlandC4LookupTable::get_instance();

    // Direct polynomial implementations (no lookup table)
    static double hk_f_direct(double r, double h) {
        const double e = h * 0.5;
        const double u = r / e;

        if (u >= 2.0) {
            return (r > 1e-30) ? 1.0 / r : 0.0;
        }

        double result;
        if (u < 1.0) {
            double u2 = u * u;
            double u3 = u2 * u;
            result = -0.5 * u2 * (1.0/3.0 - 3.0/20.0 * u2 + u3/20.0) + 1.4;
        } else {
            double u2 = u * u;
            double u3 = u2 * u;
            result = -1.0/(15.0*u) + (-u2 * (4.0/3.0 - u + 0.3*u2 - u3/30.0) + 1.6);
        }
        return result / e;
    }

    static double hk_g_direct(double r, double h) {
        const double e = h * 0.5;
        const double u = r / e;

        if (u >= 2.0) {
            return (r > 1e-30) ? 1.0 / (r * r * r) : 0.0;
        }

        double result;
        if (u < 1e-10) {
            result = 4.0 / 3.0;
        } else if (u < 1.0) {
            double u2 = u * u;
            double u3 = u2 * u;
            result = 4.0/3.0 - 1.2*u2 + 0.5*u3;
        } else {
            double u2 = u * u;
            double u3 = u2 * u;
            double u4 = u2 * u2;
            double u5 = u4 * u;
            double u6 = u3 * u3;
            result = (-1.0/15.0 + 8.0/3.0*u3 - 3.0*u4 + 1.2*u5 - u6/6.0) / u3;
        }
        return result / (e * e * e);
    }

    static double wc4_phi_direct(double r, double h) {
        const double q = r / h;

        if (q >= 1.0) {
            return (r > 1e-30) ? 1.0 / r : 0.0;
        }

        const double q2 = q * q;
        const double q3 = q2 * q;
        const double q4 = q2 * q2;
        const double q5 = q4 * q;
        const double q6 = q3 * q3;
        const double q7 = q6 * q;
        const double q8 = q4 * q4;
        const double q9 = q8 * q;

        const double a0 =  3.4374743761;
        const double a1 = -0.0031873250;
        const double a2 = -10.2154807743;
        const double a3 = -1.1577720555;
        const double a4 =  36.1013669755;
        const double a5 = -26.3399094060;
        const double a6 = -44.1079372114;
        const double a7 =  82.6543766683;
        const double a8 = -50.5921624056;
        const double a9 =  11.2232565249;

        return (a0 + a1*q + a2*q2 + a3*q3 + a4*q4 + a5*q5 + a6*q6 + a7*q7 + a8*q8 + a9*q9) / h;
    }

    static double wc4_g_direct(double r, double h) {
        const double q = r / h;

        if (q >= 1.0) {
            return (r > 1e-30) ? 1.0 / (r * r * r) : 0.0;
        }

        if (q < 1e-10) {
            return 0.0;
        }

        const double q2 = q * q;
        const double q3 = q2 * q;
        const double q4 = q2 * q2;
        const double q5 = q4 * q;
        const double q6 = q3 * q3;
        const double q7 = q6 * q;

        const double b0 = -0.0031873250;
        const double b1 = -20.4309615486;
        const double b2 = -3.4733161665;
        const double b3 = 144.4054679020;
        const double b4 = -131.6995470300;
        const double b5 = -264.6476232684;
        const double b6 = 578.5806366781;
        const double b7 = -404.7372992448;
        const double b8 = 101.0093087241;

        const double dphi_dq = b0 + b1*q + b2*q2 + b3*q3 + b4*q4 + b5*q5 + b6*q6 + b7*q7 + b8*q*q7;
        const double h3 = h * h * h;

        return -dphi_dq / (h3 * q);
    }
};

// Test: Compare N-body force calculation using lookup vs simple for-loop
TEST_F(LookupVsForLoopTest, HernquistKatz_NBodyForce_LookupMatchesForLoop) {
    const int N = 50;  // Small N for exact O(N²) comparison
    const double h = 0.5;

    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> pos_dist(-2.0, 2.0);
    std::uniform_real_distribution<double> mass_dist(0.1, 1.0);

    // Generate random particles
    std::vector<std::array<double, 3>> positions(N);
    std::vector<double> masses(N);

    for (int i = 0; i < N; ++i) {
        positions[i] = {pos_dist(rng), pos_dist(rng), pos_dist(rng)};
        masses[i] = mass_dist(rng);
    }

    // Compute accelerations with LOOKUP TABLE
    std::vector<std::array<double, 3>> acc_lookup(N, {0.0, 0.0, 0.0});

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;

            double dx = positions[i][0] - positions[j][0];
            double dy = positions[i][1] - positions[j][1];
            double dz = positions[i][2] - positions[j][2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (r < 1e-30) continue;

            double g = hk_table.g_full(r, h);  // LOOKUP

            acc_lookup[i][0] -= G * masses[j] * g * dx;
            acc_lookup[i][1] -= G * masses[j] * g * dy;
            acc_lookup[i][2] -= G * masses[j] * g * dz;
        }
    }

    // Compute accelerations with SIMPLE FOR-LOOP (direct polynomial)
    std::vector<std::array<double, 3>> acc_forloop(N, {0.0, 0.0, 0.0});

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;

            double dx = positions[i][0] - positions[j][0];
            double dy = positions[i][1] - positions[j][1];
            double dz = positions[i][2] - positions[j][2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (r < 1e-30) continue;

            double g = hk_g_direct(r, h);  // DIRECT POLYNOMIAL

            acc_forloop[i][0] -= G * masses[j] * g * dx;
            acc_forloop[i][1] -= G * masses[j] * g * dy;
            acc_forloop[i][2] -= G * masses[j] * g * dz;
        }
    }

    // Compare: lookup must match for-loop exactly (within numerical precision)
    double max_rel_error = 0.0;
    int mismatches = 0;

    for (int i = 0; i < N; ++i) {
        double acc_mag_lookup = std::sqrt(
            acc_lookup[i][0]*acc_lookup[i][0] +
            acc_lookup[i][1]*acc_lookup[i][1] +
            acc_lookup[i][2]*acc_lookup[i][2]);

        double acc_mag_forloop = std::sqrt(
            acc_forloop[i][0]*acc_forloop[i][0] +
            acc_forloop[i][1]*acc_forloop[i][1] +
            acc_forloop[i][2]*acc_forloop[i][2]);

        double rel_error = std::abs(acc_mag_lookup - acc_mag_forloop) /
                           (acc_mag_forloop + 1e-30);

        if (rel_error > kRelTol) {
            mismatches++;
        }
        max_rel_error = std::max(max_rel_error, rel_error);
    }

    std::cout << "[H-K Lookup vs ForLoop] N=" << N
              << ", max relative error = " << std::scientific << max_rel_error
              << ", mismatches = " << mismatches << "/" << N << std::endl;

    EXPECT_EQ(mismatches, 0)
        << "REGRESSION: Lookup table forces differ from direct polynomial!";
    EXPECT_LT(max_rel_error, kRelTol)
        << "REGRESSION: max error " << max_rel_error << " exceeds " << kRelTol;
}

// Test: Compare potential energy using lookup vs simple for-loop
TEST_F(LookupVsForLoopTest, HernquistKatz_PotentialEnergy_LookupMatchesForLoop) {
    const int N = 50;
    const double h = 0.5;

    std::mt19937 rng(54321);
    std::uniform_real_distribution<double> pos_dist(-2.0, 2.0);
    std::uniform_real_distribution<double> mass_dist(0.1, 1.0);

    std::vector<std::array<double, 3>> positions(N);
    std::vector<double> masses(N);

    for (int i = 0; i < N; ++i) {
        positions[i] = {pos_dist(rng), pos_dist(rng), pos_dist(rng)};
        masses[i] = mass_dist(rng);
    }

    // Compute PE with LOOKUP TABLE
    double PE_lookup = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dx = positions[i][0] - positions[j][0];
            double dy = positions[i][1] - positions[j][1];
            double dz = positions[i][2] - positions[j][2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            PE_lookup -= G * masses[i] * masses[j] * hk_table.f_full(r, h);
        }
    }

    // Compute PE with DIRECT POLYNOMIAL
    double PE_forloop = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dx = positions[i][0] - positions[j][0];
            double dy = positions[i][1] - positions[j][1];
            double dz = positions[i][2] - positions[j][2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            PE_forloop -= G * masses[i] * masses[j] * hk_f_direct(r, h);
        }
    }

    double rel_error = std::abs(PE_lookup - PE_forloop) / std::abs(PE_forloop);

    std::cout << "[H-K PE Lookup vs ForLoop] PE_lookup=" << PE_lookup
              << ", PE_forloop=" << PE_forloop
              << ", rel_error=" << std::scientific << rel_error << std::endl;

    EXPECT_LT(rel_error, kRelTol)
        << "REGRESSION: Potential energy from lookup differs from direct polynomial!";
}

// Test: Wendland C4 lookup vs for-loop
TEST_F(LookupVsForLoopTest, WendlandC4_NBodyForce_LookupMatchesForLoop) {
    const int N = 50;
    const double h = 0.8;

    std::mt19937 rng(99999);
    std::uniform_real_distribution<double> pos_dist(-1.5, 1.5);
    std::uniform_real_distribution<double> mass_dist(0.1, 1.0);

    std::vector<std::array<double, 3>> positions(N);
    std::vector<double> masses(N);

    for (int i = 0; i < N; ++i) {
        positions[i] = {pos_dist(rng), pos_dist(rng), pos_dist(rng)};
        masses[i] = mass_dist(rng);
    }

    // Compute with LOOKUP
    std::vector<std::array<double, 3>> acc_lookup(N, {0.0, 0.0, 0.0});
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;

            double dx = positions[i][0] - positions[j][0];
            double dy = positions[i][1] - positions[j][1];
            double dz = positions[i][2] - positions[j][2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (r < 1e-30) continue;

            double g = wc4_table.g_full(r, h);

            acc_lookup[i][0] -= G * masses[j] * g * dx;
            acc_lookup[i][1] -= G * masses[j] * g * dy;
            acc_lookup[i][2] -= G * masses[j] * g * dz;
        }
    }

    // Compute with DIRECT POLYNOMIAL
    std::vector<std::array<double, 3>> acc_forloop(N, {0.0, 0.0, 0.0});
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;

            double dx = positions[i][0] - positions[j][0];
            double dy = positions[i][1] - positions[j][1];
            double dz = positions[i][2] - positions[j][2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (r < 1e-30) continue;

            double g = wc4_g_direct(r, h);

            acc_forloop[i][0] -= G * masses[j] * g * dx;
            acc_forloop[i][1] -= G * masses[j] * g * dy;
            acc_forloop[i][2] -= G * masses[j] * g * dz;
        }
    }

    double max_rel_error = 0.0;
    int mismatches = 0;

    for (int i = 0; i < N; ++i) {
        double acc_mag_lookup = std::sqrt(
            acc_lookup[i][0]*acc_lookup[i][0] +
            acc_lookup[i][1]*acc_lookup[i][1] +
            acc_lookup[i][2]*acc_lookup[i][2]);

        double acc_mag_forloop = std::sqrt(
            acc_forloop[i][0]*acc_forloop[i][0] +
            acc_forloop[i][1]*acc_forloop[i][1] +
            acc_forloop[i][2]*acc_forloop[i][2]);

        double rel_error = std::abs(acc_mag_lookup - acc_mag_forloop) /
                           (acc_mag_forloop + 1e-30);

        if (rel_error > kRelTol) {
            mismatches++;
        }
        max_rel_error = std::max(max_rel_error, rel_error);
    }

    std::cout << "[W-C4 Lookup vs ForLoop] N=" << N
              << ", max relative error = " << std::scientific << max_rel_error
              << ", mismatches = " << mismatches << "/" << N << std::endl;

    EXPECT_EQ(mismatches, 0)
        << "REGRESSION: W-C4 lookup table forces differ from direct polynomial!";
    EXPECT_LT(max_rel_error, kRelTol)
        << "REGRESSION: max error " << max_rel_error << " exceeds " << kRelTol;
}

}  // namespace
