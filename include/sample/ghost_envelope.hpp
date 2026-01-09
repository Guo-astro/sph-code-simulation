/**
 * @file ghost_envelope.hpp
 * @brief Ghost boundary envelope for pressure confinement (SSOT)
 *
 * Single Source of Truth for ghost envelope particle generation.
 * Used for pressure confinement during relaxation of self-gravitating clouds.
 *
 * Physics:
 *   - Ghost particles are excluded from time integration (is_ghost = true)
 *   - They provide neighbors for edge particles to compute correct density
 *   - Pressure gradient at cloud edge prevents expansion
 *
 * First Principles (Wendland C4 kernel):
 *   - Kernel gradient peaks at q = r/h ≈ 0.5
 *   - Optimal first layer distance: 0.5h from cloud edge
 *   - Layer spacing: 0.5h for dense packing
 *   - Number of layers: 4 (to cover 0.5h to 2.0h, full kernel support)
 *
 * Reference: Derived from SPH kernel analysis for Wendland C4
 */

#pragma once

#include "particle.hpp"
#include "vector_type.hpp"
#include <vector>
#include <cmath>
#include <iostream>

namespace sph {

/**
 * @brief Configuration for ghost envelope generation
 */
struct GhostEnvelopeConfig {
    real R_cloud;           // Cloud radius [code units]
    real rho_edge;          // Density at cloud edge [code units]
    real u_envelope;        // Internal energy for envelope [code units]
    real particle_mass;     // Mass per particle [code units]
    int N_neighbor;         // Target neighbor number (typically 50)
    int num_layers = 4;     // Number of envelope layers (default: 4)
    real envelope_factor = 2.0;  // Envelope extends to R_cloud * (1 + factor * h/R)
};

/**
 * @brief Ghost envelope generator (SSOT)
 *
 * Creates ghost boundary particles for pressure confinement.
 * All formulas derived from first principles of SPH kernel theory.
 */
class GhostEnvelopeGenerator {
public:
    /**
     * @brief Compute smoothing length at edge from first principles
     *
     * From SPH neighbor condition:
     *   N_neighbor = (4/3)π h³ ρ / m
     *
     * Solving for h:
     *   h = (3 N_neighbor m / (4π ρ))^(1/3)
     */
    static real compute_h_edge(real particle_mass, real rho_edge, int N_neighbor) {
        constexpr real A_vol = 4.0 * M_PI / 3.0;
        return std::pow(N_neighbor * particle_mass / (rho_edge * A_vol), 1.0 / 3.0);
    }

    /**
     * @brief Compute optimal layer spacing from kernel gradient analysis
     *
     * Wendland C4 kernel gradient peaks at q = r/h ≈ 0.47
     * Optimal spacing: 0.5h for dense packing with maximum pressure gradient
     */
    static real compute_layer_spacing(real h_edge) {
        return 0.5 * h_edge;
    }

    /**
     * @brief Compute first layer distance from cloud edge
     *
     * First layer at R_cloud + 0.5h where Wendland C4 gradient peaks.
     * This provides maximum pressure force contribution.
     */
    static real compute_first_layer_distance(real h_edge) {
        return 0.5 * h_edge;
    }

    /**
     * @brief Generate ghost envelope particles
     *
     * Creates spherical shells of ghost particles outside cloud edge.
     * Uses Fibonacci sphere (golden spiral) for uniform distribution.
     *
     * @param config Envelope configuration
     * @return Vector of ghost particles (with is_ghost = true)
     */
    static std::vector<SPHParticle> generate(const GhostEnvelopeConfig& config) {
        std::vector<SPHParticle> envelope_particles;

        // Compute h_edge from first principles
        real h_edge = compute_h_edge(config.particle_mass, config.rho_edge, config.N_neighbor);

        // Layer spacing: 0.5h for dense packing
        real dr_layer = compute_layer_spacing(h_edge);

        // Envelope outer radius (extend to 2h from cloud edge for full kernel support)
        real R_envelope = config.R_cloud + 2.0 * h_edge;

        // Golden ratio for Fibonacci sphere
        const real golden_ratio = (1.0 + std::sqrt(5.0)) / 2.0;
        const real golden_angle = 2.0 * M_PI / (golden_ratio * golden_ratio);

        int particle_id = 0;

        for (int layer = 0; layer < config.num_layers; ++layer) {
            // Radius for this shell:
            // First layer at R_cloud + 0.5h (optimal for pressure gradient)
            // Subsequent layers at R_cloud + 1.0h, 1.5h, 2.0h, ...
            real r_shell = config.R_cloud + (layer + 1.0) * dr_layer;

            if (r_shell > R_envelope) break;

            // Particles per shell: based on target volume density ρ_edge
            // Shell volume: dV = 4πr² × dr = 4πr² × 0.5h = 2πr²h
            // Required mass: dM = ρ_edge × dV
            // Number of particles: N = dM / m_p = 2πρ_edge × r²h / m_p
            int N_shell = static_cast<int>(2.0 * M_PI * config.rho_edge * r_shell * r_shell * h_edge / config.particle_mass);
            N_shell = std::max(N_shell, 20);  // At least 20 particles per shell

            // Create particles on shell using golden spiral
            for (int i = 0; i < N_shell; ++i) {
                // Fibonacci sphere distribution
                real y = 1.0 - (2.0 * i + 1.0) / N_shell;  // -1 to 1
                real radius_at_y = std::sqrt(std::max(0.0, 1.0 - y * y));
                real theta = golden_angle * i;

                real x = radius_at_y * std::cos(theta);
                real z = radius_at_y * std::sin(theta);

                SPHParticle p;
#if DIM == 3
                p.pos = vec_t{r_shell * x, r_shell * y, r_shell * z};
#elif DIM == 2
                // For 2D, create ring instead of sphere
                real angle = 2.0 * M_PI * i / N_shell;
                p.pos = vec_t{r_shell * std::cos(angle), r_shell * std::sin(angle)};
#else
                p.pos = vec_t{r_shell};  // 1D: just the distance
#endif
                p.vel = 0.0;
                p.mass = config.particle_mass;
                p.dens = config.rho_edge;  // Assign edge density
                p.ene = config.u_envelope;
                p.pres = (1.0001 - 1.0) * config.rho_edge * config.u_envelope;  // γ ≈ 1
                p.sml = h_edge;
                p.is_ghost = true;  // CRITICAL: exclude from time integration
                p.id = particle_id++;

                envelope_particles.push_back(p);
            }
        }

        return envelope_particles;
    }

    /**
     * @brief Print envelope configuration summary
     */
    static void print_summary(const GhostEnvelopeConfig& config, int num_particles) {
        real h_edge = compute_h_edge(config.particle_mass, config.rho_edge, config.N_neighbor);
        real dr_layer = compute_layer_spacing(h_edge);

        std::cout << "\n=== Ghost Envelope Configuration (SSOT) ===" << std::endl;
        std::cout << "First principles:" << std::endl;
        std::cout << "  h_edge = (3Nm/4πρ)^(1/3) = " << h_edge << " [code]" << std::endl;
        std::cout << "  Layer spacing = 0.5h = " << dr_layer << " [code]" << std::endl;
        std::cout << "  First layer at R + 0.5h = " << config.R_cloud + 0.5*h_edge << std::endl;
        std::cout << "Layers:" << std::endl;
        for (int i = 0; i < config.num_layers; ++i) {
            real r = config.R_cloud + (i + 1.0) * dr_layer;
            std::cout << "  Layer " << i << ": r = " << r << " (q = " << (i+1)*0.5 << ")" << std::endl;
        }
        std::cout << "Created " << num_particles << " ghost particles" << std::endl;
        std::cout << "===========================================" << std::endl;
    }
};

} // namespace sph
