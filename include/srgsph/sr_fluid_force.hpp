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
 * Evolves canonical momentum S and canonical energy e using Riemann solver
 * Key features:
 * - Volume-based formulation for variable smoothing length
 * - Primitive variable recovery from conserved variables
 * - Tangential velocity coupling through Lorentz factor
 */
class FluidForce : public sph::FluidForce {
    real m_gamma;             // EOS gamma (ratio of specific heats)
    real m_c_speed;           // Speed of light (typically 1.0 in code units)
    real m_c_shock;           // Shock detection coefficient (C_shock)
    real m_c_cd;              // Contact discontinuity limiter (C_cd)
    bool m_use_muscl;         // Enable MUSCL reconstruction when gradients are available

    void solve_interface_state(
        const real left_state[5],
        const real right_state[5],
        real & pstar,
        real & vstar,
        real & vt_star) const;

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
