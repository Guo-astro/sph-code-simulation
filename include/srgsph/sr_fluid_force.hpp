#pragma once

#include <functional>
#include "fluid_force.hpp"

namespace sph
{
namespace srgsph
{

/**
 * Special Relativistic Godunov SPH Fluid Force
 * Based on Kitajima, Inutsuka, Seno (2025) arXiv:2510.18251v1
 * 
 * Evolves canonical momentum S and canonical energy e using Riemann solver
 * Key features:
 * - Volume-based formulation for variable smoothing length
 * - MUSCL reconstruction for 2nd order accuracy
 * - Monotonicity constraints for shock handling
 * - Primitive variable recovery from conserved variables
 */
class FluidForce : public sph::FluidForce {
    bool m_is_2nd_order;      // Enable MUSCL reconstruction
    real m_gamma;             // EOS gamma (ratio of specific heats)
    real m_c_speed;           // Speed of light (typically 1.0 in code units)
    real m_c_shock;           // Shock detection parameter (default: 3.0)
    real m_c_cd;              // Contact discontinuity parameter (default: 1.0)

    // Riemann solver: (left[4], right[4]) -> (pstar, vstar)
    // Input: [velocity, density, pressure, sound_speed]
    // Output: pstar (interface pressure), vstar (interface velocity)
    std::function<void(const real[], const real[], real & pstar, real & vstar)> m_solver;

    /**
     * Setup exact Riemann solver (Pons et al. 2000)
     * Exact solution for special relativistic hydrodynamics
     * with arbitrary tangential velocities and ideal gas EOS.
     */
    void exact_riemann_solver();
    
    /**
     * Compute post-shock velocity using Rankine-Hugoniot relations
     * @param p_a Pressure ahead of shock
     * @param rho_a Density ahead of shock
     * @param h_a Specific enthalpy ahead of shock
     * @param vx_a Normal velocity ahead of shock
     * @param W_a Lorentz factor ahead of shock
     * @param p_b Pressure behind shock
     * @param is_left True for left-propagating shock
     * @return Post-shock normal velocity
     */
    real compute_shock_velocity(const real p_a, const real rho_a, const real h_a,
                                const real vx_a, const real W_a, const real p_b,
                                const bool is_left) const;
    
    /**
     * Compute post-rarefaction velocity by integrating ODE
     * @param p_a Pressure ahead of rarefaction
     * @param rho_a Density ahead of rarefaction
     * @param h_a Specific enthalpy ahead of rarefaction
     * @param vx_a Normal velocity ahead of rarefaction
     * @param cs_a Sound speed ahead of rarefaction
     * @param p_b Pressure behind rarefaction
     * @param is_left True for left-propagating rarefaction
     * @return Post-rarefaction normal velocity
     */
    real compute_rarefaction_velocity(const real p_a, const real rho_a, const real h_a,
                                      const real vx_a, const real cs_a, const real p_b,
                                      const bool is_left) const;

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    
    /**
     * Main calculation loop
     * Computes dS/dt (canonical momentum) and de/dt (canonical energy)
     * using Riemann solver with MUSCL reconstruction
     */
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}
