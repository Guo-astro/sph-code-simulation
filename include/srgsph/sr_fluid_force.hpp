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
     * Setup iterative Riemann solver for SR-GSPH
     * Based on Marti & Mueller (1994) exact solver
     * Uses rest-frame baryon density and relativistic enthalpy
     */
    void exact_riemann_solver();

    /**
     * Setup relativistic iterative Riemann solver
     * Adapted from van Leer (1997) for special relativity
     * Uses Newton-Raphson iteration with relativistic wave impedances
     */
    void iterative_solver();

    /**
     * Compute velocity behind shock or rarefaction wave
     * Based on Marti & Mueller (1994), Eqs. 67-86 for shocks, analytic rarefaction
     * 
     * @param p Post-wave pressure
     * @param rho_a Pre-wave rest-frame density
     * @param p_a Pre-wave pressure
     * @param h_a Pre-wave specific enthalpy
     * @param cs_a Pre-wave sound speed
     * @param vel_a Pre-wave velocity
     * @param w_a Pre-wave Lorentz factor
     * @param sign Direction: -1 for left wave, +1 for right wave
     * @param[out] rho Post-wave rest-frame density
     * @param[out] h Post-wave specific enthalpy
     * @param[out] cs Post-wave sound speed
     * @param[out] vel Post-wave velocity
     */
    void get_velocity_behind_wave(
        const real p,
        const real rho_a, const real p_a, const real h_a, const real cs_a,
        const real vel_a, const real w_a,
        const real sign,
        real & rho, real & h, real & cs, real & vel) const;

    /**
     * Compute velocity difference between left and right intermediate states
     * Used in root finding for pressure
     * 
     * @param p Pressure to test
     * @param rhol Left rest-frame density
     * @param pl Left pressure
     * @param hl Left specific enthalpy
     * @param csl Left sound speed
     * @param vell Left velocity
     * @param wl Left Lorentz factor
     * @param rhor Right rest-frame density
     * @param pr Right pressure
     * @param hr Right specific enthalpy
     * @param csr Right sound speed
     * @param velr Right velocity
     * @param wr Right Lorentz factor
     * @param[out] vel_l Velocity in left intermediate state
     * @param[out] vel_r Velocity in right intermediate state
     * @return Velocity difference (vel_l - vel_r)
     */
    real get_velocity_difference(
        const real p,
        const real rhol, const real pl, const real hl, const real csl, const real vell, const real wl,
        const real rhor, const real pr, const real hr, const real csr, const real velr, const real wr,
        real & vel_l, real & vel_r) const;

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
