/**
 * @file point_mass_bh.hpp
 * @brief External gravitational force from a point-mass black hole
 * 
 * Implements Newtonian gravity from a massive point source (e.g., IMBH, SMBH)
 * with optional softening to prevent singularities near the BH.
 * 
 * Physics:
 *   F = -G M_BH (r - r_BH) / (|r - r_BH|^2 + eps^2)^(3/2)
 *   Φ = -G M_BH / sqrt(|r - r_BH|^2 + eps^2)
 * 
 * Where:
 *   - G: Gravitational constant [code units]
 *   - M_BH: Black hole mass [M_☉]
 *   - r_BH: Black hole position [pc]
 *   - eps: Softening length [pc] (prevents singularity)
 * 
 * Use cases:
 *   - IMBH tidal disruption of molecular clouds
 *   - SMBH accretion disk simulations
 *   - Binary BH + gas disk interactions
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
 */
struct PointMassBHParams
{
    vec_t position;        ///< BH position [code units]
    real mass;             ///< BH mass [M_☉]
    real softening_length; ///< Gravitational softening ε [pc]
    real G_constant;       ///< Gravitational constant [code units]
    bool is_moving;        ///< Whether BH position changes with time
    vec_t velocity;        ///< BH velocity if moving [code units/time]
    
    PointMassBHParams()
        : position(0.0)
        , mass(1e5)  // Default 10^5 M_☉ (IMBH)
        , softening_length(0.01)
        , G_constant(1.0)
        , is_moving(false)
        , velocity(0.0)
    {}
};

/**
 * @brief External gravitational force from a point-mass black hole
 * 
 * This module computes the acceleration and potential from a single
 * massive point source, typically representing an intermediate-mass
 * or supermassive black hole.
 * 
 * Key features:
 *   - Plummer softening to prevent singularities
 *   - Optional BH motion (for binary BH or moving frame)
 *   - Energy diagnostics (accretion rate, tidal heating)
 * 
 * Usage:
 * ```cpp
 * auto bh_force = std::make_shared<PointMassBlackHole>();
 * PointMassBHParams params;
 * params.mass = 1e5;  // 10^5 M_☉
 * params.position = vec_t{10.0, 0.0, 0.0};  // 10 pc offset
 * params.softening_length = 0.05;  // 0.05 pc softening
 * bh_force->initialize(params);
 * bh_force->apply_force(sim);
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
     * Updates particle accelerations (p[i].acc) and potentials (p[i].phi)
     */
    void calculation(std::shared_ptr<Simulation> sim) override;
    
    /**
     * @brief Update BH position for moving black hole
     * @param dt Time step [code units]
     * 
     * Only effective if is_moving == true
     */
    void update_position(real dt);
    
    /**
     * @brief Get current BH position
     */
    vec_t get_position() const { return m_position; }
    
    /**
     * @brief Get BH mass
     */
    real get_mass() const { return m_mass; }
    
    /**
     * @brief Get total mass within radius r from BH
     * @param sim Simulation state
     * @param r Radius [code units]
     * @return Mass within r [M_☉]
     * 
     * Useful for measuring accretion rate and bound mass fraction
     */
    real get_enclosed_mass(std::shared_ptr<Simulation> sim, real r) const;
    
    /**
     * @brief Count particles within accretion radius
     * @param sim Simulation state
     * @param r_acc Accretion radius [code units]
     * @return Number of particles within r_acc
     */
    int count_accreted_particles(std::shared_ptr<Simulation> sim, real r_acc) const;
    
    /**
     * @brief Compute acceleration from BH at position r
     * @param r Particle position [code units]
     * @return Acceleration vector [code units/time²]
     */
    vec_t acceleration(const vec_t& r) const;
    
    /**
     * @brief Compute gravitational potential from BH at position r
     * @param r Particle position [code units]
     * @return Potential [code units²/time²]
     */
    real potential(const vec_t& r) const;

private:
    vec_t m_position;           ///< Current BH position [code units]
    real m_mass;                ///< BH mass [M_☉]
    real m_softening_length;    ///< Plummer softening ε [pc]
    real m_softening_squared;   ///< ε² (cached for performance)
    real m_G_constant;          ///< Gravitational constant
    bool m_is_moving;           ///< Whether BH moves
    vec_t m_velocity;           ///< BH velocity if moving
    bool m_initialized;         ///< Initialization flag
    
    /**
     * @brief Softened gravitational force magnitude
     * @param r_squared Distance squared from BH
     * @return Force magnitude factor: G M / (r² + ε²)^(3/2)
     */
    inline real softened_force_factor(real r_squared) const
    {
        const real denom = r_squared + m_softening_squared;
        return m_G_constant * m_mass / (denom * std::sqrt(denom));
    }
};

} // namespace external_forces
} // namespace sph
