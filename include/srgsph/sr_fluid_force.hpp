#pragma once

#include "fluid_force.hpp"
#include "particle.hpp"
#include "exception.hpp"

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
 * Variable naming follows Kitajima notation (SRGSPH paper Section 2.1):
 *   N       = lab-frame baryon number density (N = γn)
 *   n       = rest-frame baryon number density
 *   γ       = Lorentz factor
 *   H       = enthalpy per baryon
 *   P       = pressure
 *   v^x     = normal velocity component
 *   v^t     = tangential velocity component
 *   c_s     = sound speed
 *   γ_c     = ratio of specific heats
 *   S       = canonical momentum per baryon
 *   e       = canonical energy per baryon
 * 
 * Evolves canonical momentum S and canonical energy e using Riemann solver
 * Key features:
 * - Volume-based formulation for variable smoothing length
 * - Primitive variable recovery from conserved variables
 * - Tangential velocity coupling through Lorentz factor
 */
class FluidForce : public sph::FluidForce {
    real m_gamma;             // Ratio of specific heats γ_c
    real m_c_speed;           // Speed of light c (typically 1.0 in code units)
    real m_c_shock;           // Shock detection coefficient C_shock
    real m_c_cd;              // Contact discontinuity limiter C_cd
    bool m_use_muscl;         // Enable MUSCL reconstruction when gradients available

    /**
     * Solve Riemann problem at particle interface
     * 
     * @param left_state  Left state [v^x, n, P, c_s, v^t]
     * @param right_state Right state [v^x, n, P, c_s, v^t]
     * @param P_star      [output] Interface pressure P*
     * @param v_x_star    [output] Interface normal velocity v^x*
     * @param v_t_star    [output] Interface tangential velocity v^t*
     */
    void solve_interface_state(
        const real left_state[5],
        const real right_state[5],
        real & P_star,
        real & v_x_star,
        real & v_t_star) const;

public:
    /**
     * Scale SR-GSPH derivatives by baryon number so dS and de are per-baryon
     * quantities. Also keeps the legacy acc/dene aliases in sync.
     */
    static void normalize_sr_derivatives(SPHParticle & particle);

    void initialize(std::shared_ptr<SPHParameters> param) override;
    
    /**
     * Main calculation loop
     * Computes dS/dt (canonical momentum) and de/dt (canonical energy)
     * using iterative Riemann solver (1st order)
     */
    void calculation(std::shared_ptr<Simulation> sim) override;
};

inline void FluidForce::normalize_sr_derivatives(SPHParticle & particle)
{
    if(particle.nu <= 0.0) {
        THROW_ERROR("SRGSPH particle baryon number must be positive.");
    }

    const real inv_nu = 1.0 / particle.nu;
    particle.dS *= inv_nu;
    particle.de *= inv_nu;
    particle.acc = particle.dS;
    particle.dene = particle.de;
}

}
}
