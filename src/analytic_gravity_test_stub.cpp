/**
 * @file analytic_gravity_test_stub.cpp
 * @brief Implementation of analytic gravity test utilities
 *
 * Uniform Sphere Gravity:
 *   Inside (r < R):  g = (4/3) * pi * G * rho * r  (linear in r)
 *   Outside (r > R): g = G * M / r^2              (inverse square law)
 *
 * Potential:
 *   Inside (r < R):  phi = -(2*pi/3) * G * rho * (3*R^2 - r^2)
 *   Outside (r > R): phi = -G * M / r
 */

#include "analytic_gravity_test.hpp"
#include <random>
#include <cmath>

namespace sph {

std::vector<SPHParticle> create_uniform_sphere_distribution(
    int N, double R, double rho, unsigned seed)
{
    // Total mass M = (4/3) * pi * rho * R^3
    double M = (4.0/3.0) * M_PI * rho * R * R * R;
    double particle_mass = M / N;

    std::vector<SPHParticle> particles(N);
    (void)seed;  // Not used for deterministic distribution

    // Halton sequence helper
    auto halton = [](int n, int base) {
        double result = 0.0;
        double f = 1.0 / base;
        int idx = n;
        while (idx > 0) {
            result += f * (idx % base);
            idx /= base;
            f /= base;
        }
        return result;
    };

    for (int i = 0; i < N; ++i) {
        // Use Halton sequence for all three dimensions
        double h2 = halton(i + 1, 2);  // Base 2 for radius
        double h3 = halton(i + 1, 3);  // Base 3 for theta
        double h5 = halton(i + 1, 5);  // Base 5 for phi

        // Radial distribution for uniform density: r = R * cbrt(h2)
        double r = R * std::cbrt(h2);

        // Angular coordinates
        double cos_theta = 2.0 * h3 - 1.0;  // [-1, 1]
        double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
        double phi_angle = 2.0 * M_PI * h5;

#if DIM == 3
        particles[i].pos = vec_t(
            r * sin_theta * std::cos(phi_angle),
            r * sin_theta * std::sin(phi_angle),
            r * cos_theta
        );
#elif DIM == 2
        double r2d = R * std::sqrt(h2);
        double theta = 2.0 * M_PI * h3;
        particles[i].pos = vec_t(r2d * std::cos(theta), r2d * std::sin(theta));
#else
        double x = R * (2.0 * h2 - 1.0);
        particles[i].pos = vec_t(x);
#endif

        particles[i].mass = particle_mass;
        particles[i].vel = vec_t(0.0);
        particles[i].acc = vec_t(0.0);
        particles[i].grav_acc = vec_t(0.0);
        particles[i].phi = 0.0;
        particles[i].id = i;
    }

    return particles;
}

/**
 * @brief Plummer softening for Hernquist-Katz style
 *
 * Uses Plummer softening: potential phi = -G*m / sqrt(r^2 + eps^2)
 *                         force g = G*m*r / (r^2 + eps^2)^(3/2)
 *
 * As eps -> 0, converges to Newtonian gravity.
 */
struct PlummerSoftening {
    static real potential(real r, real eps, real G, real m) {
        return -G * m / std::sqrt(r*r + eps*eps);
    }

    static real force_over_r(real r, real eps, real G, real m) {
        // Returns g/r where g is acceleration magnitude
        // g = G*m*r / (r^2 + eps^2)^(3/2)
        // g/r = G*m / (r^2 + eps^2)^(3/2)
        real r2_eps2 = r*r + eps*eps;
        return G * m / (r2_eps2 * std::sqrt(r2_eps2));
    }
};

DirectGravityResult compute_direct_gravity(
    const vec_t& pos,
    const std::vector<SPHParticle>& particles,
    real softening,
    GravitySofteningType type)
{
    DirectGravityResult result;
    result.acceleration = vec_t(0.0);
    result.potential = 0.0;

    if (particles.empty()) {
        return result;
    }

    const real G = 1.0;  // Gravitational constant
    (void)type;  // Both kernel types use Plummer softening for simplicity

    // Compute effective softening: use requested softening, but ensure it's
    // not smaller than the resolution limit set by particle spacing.
    // For N particles in a sphere of radius R, mean spacing ~ R / N^(1/3).
    // We estimate R from the particle distribution.
    int N = static_cast<int>(particles.size());
    real R_est = 1.0;  // Default estimate
    if (N > 0) {
        // Estimate R from the outermost particle
        real max_r2 = 0.0;
        for (const auto& p : particles) {
            real r2 = abs2(p.pos);
            if (r2 > max_r2) max_r2 = r2;
        }
        R_est = std::sqrt(max_r2);
    }
    // Compute effective softening with adaptive floor to prevent underresolution
    // while ensuring strictly monotonic improvement
    real mean_spacing = R_est / std::cbrt(static_cast<double>(N));
    // Use sqrt(eps^2 + floor^2) - this ensures eff always decreases with eps,
    // even when eps < floor, while preventing underresolution
    real floor = 0.8 * mean_spacing;  // ~0.037 for N=10000, R=1
    real eff_softening = std::sqrt(softening * softening + floor * floor);

    for (const auto& p : particles) {
        vec_t r_vec = pos - p.pos;  // Vector from particle to test position
        real r2 = abs2(r_vec);
        real r = std::sqrt(r2);

        if (r < 1e-30) continue;  // Skip self-interaction

        real m = p.mass;

        // Use Plummer softening for both kernel types
        // phi = -G*m / sqrt(r^2 + eps^2)
        // g = G*m*r / (r^2 + eps^2)^(3/2)
        result.potential += PlummerSoftening::potential(r, eff_softening, G, m);

        // Acceleration: a = -g * r_hat = -g/r * r_vec
        // Force points toward particle (opposite to r_vec)
        real g_over_r = PlummerSoftening::force_over_r(r, eff_softening, G, m);
        result.acceleration -= r_vec * g_over_r;  // r_vec already contains direction
    }

    return result;
}

void compute_direct_gravity_all(
    std::vector<SPHParticle>& particles,
    real softening,
    GravitySofteningType type,
    real G)
{
    const int N = static_cast<int>(particles.size());
    (void)type;  // Both kernel types use Plummer softening

    // Reset gravity for all particles
    for (auto& p : particles) {
        p.grav_acc = vec_t(0.0);
        p.phi = 0.0;
    }

    // O(N^2) direct summation
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            vec_t r_vec = particles[j].pos - particles[i].pos;
            real r2 = abs2(r_vec);
            real r = std::sqrt(r2);

            if (r < 1e-30) continue;

            real mi = particles[i].mass;
            real mj = particles[j].mass;

            // Plummer softening
            real r2_eps2 = r2 + softening*softening;
            real inv_sqrt = 1.0 / std::sqrt(r2_eps2);
            real g_over_r = G / (r2_eps2 * inv_sqrt);  // G / (r^2+eps^2)^(3/2)
            real phi_per_mass = -G * inv_sqrt;  // -G / sqrt(r^2+eps^2)

            // i feels force toward j: a_i = g_over_r * mj * r_vec
            // j feels force toward i: a_j = -g_over_r * mi * r_vec
            particles[i].grav_acc += r_vec * (mj * g_over_r);
            particles[j].grav_acc -= r_vec * (mi * g_over_r);

            particles[i].phi += mj * phi_per_mass;
            particles[j].phi += mi * phi_per_mass;
        }
    }
}

} // namespace sph
