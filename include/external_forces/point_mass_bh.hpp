/**
 * @file point_mass_bh.hpp
 * @brief External gravitational force from a fixed point-mass black hole with sink accretion
 *
 * Implements Newtonian gravity from a massive point source (e.g., IMBH, SMBH)
 * with Plummer softening to prevent singularities near the BH.
 *
 * Physics:
 *   Softened acceleration (Plummer form):
 *     a_BH(r) = -G M_BH r / (r^2 + eps^2)^(3/2)
 *
 *   Softened potential:
 *     Φ_BH(r) = -G M_BH / sqrt(r^2 + eps^2)
 *
 * Sink/Accretion Model (Fixed BH):
 *   When a particle enters r < r_sink, check:
 *     1. Inflow check: v_rad = v · r̂ < 0 (moving inward)
 *     2. Bound check: E = 0.5|v|^2 + Φ_BH < 0 (gravitationally bound)
 *   If both conditions met, particle is removed (accreted).
 *   BH mass remains fixed (fixed-potential approximation).
 *
 * Parameters:
 *   - G: Gravitational constant [code units]
 *   - M_BH: Black hole mass [M_☉] (fixed, does not change)
 *   - r_BH: Black hole position [pc] (typically fixed at origin)
 *   - eps: Softening length [pc] (prevents singularity, ~0.001 pc typical)
 *   - r_sink: Sink radius [pc] (accretion boundary, ~0.005 pc typical)
 *
 * Time stepping note:
 *   When using external BH gravity, timestep must include gravity constraint:
 *     dt_grav = eta_t * sqrt(h / |a_grav|)
 *   where eta_t ~ 0.1-0.2 for stability near the BH.
 *
 * Use cases:
 *   - IMBH tidal disruption of molecular clouds (Oka et al. 2017)
 *   - SMBH accretion disk simulations
 *   - Binary BH + gas disk interactions
 *
 * References:
 *   - Plummer (1911): Softening model
 *   - Oka et al. (2017): CO-0.40-0.22 IMBH candidate
 */

#pragma once

#include "defines.hpp"
#include "module.hpp"
#include "vector_type.hpp"
#include <memory>

namespace sph
{

class SPHParticle;
class Simulation;

namespace external_forces
{

/**
 * @brief Parameters for point-mass black hole external force
 *
 * Default values match typical IMBH cloud-scattering simulations:
 *   - M_BH = 10^5 M_☉ (Oka et al. 2017)
 *   - eps = 0.001 pc (~200 AU) - softening prevents singularity
 *   - r_sink = 0.005 pc (~1000 AU) - accretion boundary
 */
struct PointMassBHParams
{
    vec_t position;        ///< BH position [code units] (typically origin)
    real mass;             ///< BH mass [M_☉] (fixed, does not change on accretion)
    real softening_length; ///< Gravitational softening ε [pc] (~0.001 pc typical)
    real sink_radius;      ///< Sink/accretion radius r_sink [pc] (~0.005 pc typical)
    real G_constant;       ///< Gravitational constant [code units]
    bool is_moving;        ///< Whether BH position changes with time (false for fixed BH)
    vec_t velocity;        ///< BH velocity if moving [code units/time]
    bool enable_sink;      ///< Enable sink accretion (particle removal)

    PointMassBHParams()
        : position(0.0)
        , mass(1e5)              // Default 10^5 M_☉ (IMBH)
        , softening_length(0.001) // Default 0.001 pc (~200 AU)
        , sink_radius(0.005)      // Default 0.005 pc (~1000 AU)
        , G_constant(1.0)
        , is_moving(false)        // Fixed BH by default
        , velocity(0.0)
        , enable_sink(true)       // Enable sink by default
    {}
};

/**
 * @brief External gravitational force from a fixed point-mass black hole with sink accretion
 *
 * This module computes the acceleration and potential from a single
 * massive point source, typically representing an intermediate-mass
 * or supermassive black hole fixed at the origin.
 *
 * Key features:
 *   - Plummer softening to prevent singularities: a = -GM r / (r² + ε²)^(3/2)
 *   - Sink accretion with bound+inflow criteria
 *   - Accretion diagnostics (mass removed, particle count)
 *   - Optional BH motion (for binary BH or moving frame scenarios)
 *
 * Sink accretion algorithm (apply_sink_accretion):
 *   For each particle with r < r_sink:
 *     1. Inflow check: v_rad = v · r̂ < 0 (must be moving toward BH)
 *     2. Bound check: E = 0.5|v|² + Φ_BH(r) < 0 (must be gravitationally bound)
 *   If both conditions satisfied, particle is marked for removal.
 *   BH mass is NOT updated (fixed-potential approximation).
 *
 * Usage:
 * ```cpp
 * auto bh_force = std::make_shared<PointMassBlackHole>();
 * PointMassBHParams params;
 * params.mass = 1e5;               // 10^5 M_☉
 * params.position = vec_t{0,0,0};  // Fixed at origin
 * params.softening_length = 0.001; // 0.001 pc softening
 * params.sink_radius = 0.005;      // 0.005 pc sink radius
 * params.enable_sink = true;
 * bh_force->initialize(params);
 *
 * // In timestep loop:
 * bh_force->calculation(sim);      // Apply gravity
 * bh_force->apply_sink_accretion(sim);  // Remove accreted particles
 * ```
 */
class PointMassBlackHole : public Module
{
public:
    PointMassBlackHole();
    ~PointMassBlackHole() = default;

    /**
     * @brief Initialize BH force with parameters
     * @param params Black hole configuration
     */
    void initialize(const PointMassBHParams& params);

    /**
     * @brief Initialize from SPH parameter structure
     * @param param SPH simulation parameters
     *
     * Reads from param->external_force settings
     */
    void initialize(std::shared_ptr<SPHParameters> param) override;

    /**
     * @brief Apply BH gravitational force to all particles
     * @param sim Simulation state
     *
     * Computes softened Plummer acceleration:
     *   a_BH = -G M_BH r / (r² + ε²)^(3/2)
     * and adds to particle accelerations (p[i].acc).
     *
     * NOTE: Does NOT modify p[i].phi. External BH potential is computed
     * separately in energy calculations with correct factor (no 0.5).
     */
    void calculation(std::shared_ptr<Simulation> sim) override;

    /**
     * @brief Apply sink accretion - remove particles that fall into the BH
     * @param sim Simulation state
     * @return Number of particles removed this call
     *
     * For each particle with r < r_sink, checks:
     *   1. Inflow: v_rad = v · r̂ < 0 (moving toward BH)
     *   2. Bound: E = 0.5|v|² + Φ_BH(r) < 0 (gravitationally bound)
     *
     * If both conditions satisfied, particle is removed from simulation.
     * Removed mass is tracked in m_accreted_mass (diagnostic only).
     * BH mass is NOT updated (fixed-potential approximation).
     *
     * Call this AFTER calculation() in the timestep loop.
     */
    int apply_sink_accretion(std::shared_ptr<Simulation> sim);

    /**
     * @brief Update BH position for moving black hole
     * @param dt Time step [code units]
     *
     * Only effective if is_moving == true
     */
    void update_position(real dt);

    // === Accessors ===

    vec_t get_position() const { return m_position; }
    real get_mass() const { return m_mass; }
    real get_softening() const { return m_softening_length; }
    real get_sink_radius() const { return m_sink_radius; }
    bool is_sink_enabled() const { return m_enable_sink; }

    // === Accretion diagnostics ===

    /** @brief Get total mass removed by sink accretion [code units] */
    real get_accreted_mass() const { return m_accreted_mass; }

    /** @brief Get total number of particles removed by sink accretion */
    int get_accreted_count() const { return m_accreted_count; }

    /** @brief Reset accretion counters (e.g., at simulation start) */
    void reset_accretion_counters() { m_accreted_mass = 0.0; m_accreted_count = 0; }

    // === Diagnostic methods ===

    /**
     * @brief Get total particle mass within radius r from BH
     * @param sim Simulation state
     * @param r Radius [code units]
     * @return Mass within r [code units]
     */
    real get_enclosed_mass(std::shared_ptr<Simulation> sim, real r) const;

    /**
     * @brief Count particles within specified radius
     * @param sim Simulation state
     * @param r_acc Radius [code units]
     * @return Number of particles within r_acc
     */
    int count_particles_within(std::shared_ptr<Simulation> sim, real r_acc) const;

    /**
     * @brief Compute acceleration from BH at position r
     * @param r Particle position [code units]
     * @return Acceleration vector: -G M_BH r / (r² + ε²)^(3/2)
     */
    vec_t acceleration(const vec_t& r) const;

    /**
     * @brief Compute gravitational potential from BH at position r
     * @param r Particle position [code units]
     * @return Potential: -G M_BH / sqrt(r² + ε²)
     */
    real potential(const vec_t& r) const;

private:
    // === BH parameters ===
    vec_t m_position;           ///< Current BH position [code units]
    real m_mass;                ///< BH mass [code units] (fixed)
    real m_softening_length;    ///< Plummer softening ε [code units]
    real m_softening_squared;   ///< ε² (cached for performance)
    real m_sink_radius;         ///< Sink/accretion radius [code units]
    real m_sink_radius_squared; ///< r_sink² (cached for performance)
    real m_G_constant;          ///< Gravitational constant
    bool m_is_moving;           ///< Whether BH moves
    vec_t m_velocity;           ///< BH velocity if moving
    bool m_enable_sink;         ///< Enable sink accretion
    bool m_initialized;         ///< Initialization flag

    // === Accretion tracking (diagnostics only, BH mass fixed) ===
    real m_accreted_mass;       ///< Total mass removed by sink [code units]
    int m_accreted_count;       ///< Total particles removed by sink

    /**
     * @brief Softened gravitational force magnitude factor
     * @param r_squared Distance squared from BH
     * @return Force factor: G M / (r² + ε²)^(3/2)
     */
    inline real softened_force_factor(real r_squared) const
    {
        const real denom = r_squared + m_softening_squared;
        return m_G_constant * m_mass / (denom * std::sqrt(denom));
    }
};

} // namespace external_forces
} // namespace sph
