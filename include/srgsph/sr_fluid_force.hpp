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
 * Extended with tangential velocity support from Pons, Martí & Müller (1999)
 * "The exact solution of the Riemann problem with non-zero tangential velocities
 * in relativistic hydrodynamics" (J. Fluid Mech. 1999)
 * 
 * Evolves canonical momentum S and canonical energy e using Riemann solver
 * Key features:
 * - Volume-based formulation for variable smoothing length
 * - Primitive variable recovery from conserved variables
 * - Tangential velocity coupling through Lorentz factor
 */
class FluidForce : public sph::FluidForce {
    real m_gamma;             // EOS gamma (ratio of specific heats)
    real m_c_speed;           // Speed of light (typically 1.0 in code units)
    bool m_use_sqrt2_h;       // Whether to use √2·h factor (true for Gaussian kernel)

    // Riemann solver: (left[5], right[5]) -> (pstar, vstar)
    // Input: [normal_velocity, density, pressure, sound_speed, tangential_velocity]
    // Output: pstar (interface pressure), vstar (interface normal velocity)
    std::function<void(const real[], const real[], real & pstar, real & vstar)> m_solver;

    /**
     * Setup relativistic iterative Riemann solver
     * Extended for tangential velocities following Pons et al. (1999)
     * Uses Newton-Raphson iteration with relativistic wave impedances
     * Accounts for tangential velocity coupling through Lorentz factor
     */
    void iterative_solver();

    /**
     * Setup relativistic HLLE Riemann solver
     * More robust than iterative solver, uses wave speed estimates
     * Based on Einfeldt et al. (1991) and Mignone & Bodo (2005)
     * Slightly more dissipative but never fails
     */
    void hlle_solver();

    /**
     * Setup relativistic HLLC Riemann solver (RECOMMENDED)
     * Best balance of accuracy and robustness
     * Based on Mignone & Bodo (2005) JCP 207:545-563
     * - Resolves contact discontinuities (unlike HLL/HLLE)
     * - Always robust (unlike exact solver)
     * - Nearly as accurate as exact solver
     */
    void hllc_solver();

    /**
     * Compute velocity behind shock or rarefaction wave
     * Extended for tangential velocities following Pons et al. (1999)
     * 
     * For shocks: Uses Taub adiabat and Rankine-Hugoniot conditions
     *   Conservation: h*W*v^t/D = const (eqs. 54-55)
     * For rarefactions: Uses isentropic relations
     *   Conservation: h*W*v^t = const (eqs. 27-28)
     * 
     * @param P Post-wave pressure
     * @param na Pre-wave rest-frame density
     * @param Pa Pre-wave pressure
     * @param Ha Pre-wave specific enthalpy
     * @param csa Pre-wave sound speed
     * @param va Pre-wave normal velocity (v^x)
     * @param vta Pre-wave tangential velocity magnitude (|v^t|)
     * @param wa Pre-wave Lorentz factor (includes tangential velocity)
     * @param sign Direction: -1 for left wave, +1 for right wave
     * @param[out] n Post-wave rest-frame density
     * @param[out] H Post-wave specific enthalpy
     * @param[out] cs Post-wave sound speed
     * @param[out] vel Post-wave normal velocity
     */
    void get_velocity_behind_wave(
        const real P,
        const real na, const real Pa, const real Ha, const real csa,
        const real va, const real vta, const real wa,
        const real sign,
        real & n, real & H, real & cs, real & vel) const;

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    
    /**
     * Main calculation loop
     * Computes dS/dt (canonical momentum) and de/dt (canonical energy)
     * using iterative Riemann solver (1st order)
     */
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}
