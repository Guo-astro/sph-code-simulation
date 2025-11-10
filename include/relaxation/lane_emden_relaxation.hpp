#pragma once

#include <memory>
#include "relaxation/lane_emden_data.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "defines.hpp"

namespace sph
{

/**
 * @brief Parameters for Lane-Emden relaxation
 */
struct LaneEmdenRelaxationParams
{
    real alpha_scaling;  // Scaling factor: r = alpha * xi
    real rho_center;     // Central density
    real K;              // Polytropic constant
    real R;              // Physical radius
    real M_total;        // Total mass
    real G;              // Gravitational constant
    real gamma;          // Adiabatic index (5/3 for n=1.5)
};

/**
 * @brief Applies relaxation forces to drive particles toward Lane-Emden equilibrium
 * 
 * This module computes analytical forces based on the Lane-Emden solution
 * and applies them to particles to help them settle into hydrostatic equilibrium.
 * Velocities are damped to zero during relaxation.
 */
class LaneEmdenRelaxation
{
public:
    LaneEmdenRelaxation();
    
    /**
     * @brief Initialize with Lane-Emden data
     */
    void initialize(const LaneEmdenRelaxationParams& params);
    
    /**
     * @brief Compute relaxation force for a single particle
     * @param p Particle
     * @return Acceleration vector toward equilibrium
     */
    vec_t compute_relaxation_force(const SPHParticle& p) const;
    
    /**
     * @brief Apply relaxation forces to all particles and damp velocities
     * @param sim Simulation object containing particles
     * @param damping_factor Velocity damping (0-1, default 0.9 = strong damping)
     */
    void apply_relaxation(std::shared_ptr<Simulation> sim, real damping_factor = 0.9);
    
    /**
     * @brief Check if relaxation is initialized
     */
    bool is_initialized() const { return m_initialized; }
    
private:
    LaneEmdenData m_data;
    LaneEmdenRelaxationParams m_params;
    bool m_initialized;
};

} // namespace sph
