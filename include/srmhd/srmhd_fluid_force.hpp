#pragma once

#include "fluid_force.hpp"
#include "particle.hpp"
#include "parameters.hpp"
#include "exception.hpp"
#include "srgsph/sr_exact_riemann.hpp"

namespace sph
{
namespace srmhd
{

/**
 * Special Relativistic Magnetohydrodynamics Godunov SPH Fluid Force
 *
 * Combines:
 * - Kitajima, Inutsuka, Seno (2025) SR-GSPH formulation
 * - Iwasaki & Inutsuka (2011) GSPMHD method for MHD
 *
 * Key features:
 * - Relativistic hydrodynamics with Lorentz factor, enthalpy
 * - MHD stress tensor with magnetic pressure and tension
 * - Riemann solver for fast magnetosonic waves (compressive)
 * - Method of Characteristics for Alfven waves (incompressible)
 * - Powell divergence cleaning for numerical stability
 *
 * Variable naming follows combined notation:
 *   N       = lab-frame baryon number density (N = gamma * n)
 *   n       = rest-frame baryon number density
 *   gamma   = Lorentz factor
 *   H       = enthalpy per baryon = 1 + u/c^2 + P/(n*c^2)
 *   P       = gas pressure
 *   P_t     = total pressure = P + B_perp^2/2
 *   B       = magnetic field (always 3D, even in 1D simulations)
 *   v       = 3-velocity
 *   S       = canonical momentum per baryon = gamma * H * v
 *   e       = canonical energy per baryon = gamma * H - P/(N*c^2)
 *   c_s     = relativistic sound speed
 *   c_A     = relativistic Alfven speed
 *   c_f     = relativistic fast magnetosonic speed
 */
class FluidForce : public sph::FluidForce {
    real m_gamma;             // Ratio of specific heats (adiabatic index)
    real m_c_speed;           // Speed of light c (typically 1.0)
    real m_c_shock;           // Shock detection coefficient
    real m_c_cd;              // Contact discontinuity limiter
    bool m_use_muscl;         // Enable MUSCL reconstruction
    bool m_use_powell;        // Enable Powell divergence cleaning
    bool m_use_mhd;           // Enable MHD (false = pure hydro with SR-GSPH solver)
    RiemannSolverType m_riemann_solver;  // Riemann solver type

    /**
     * HLL-type Riemann solver for relativistic fast magnetosonic waves
     *
     * Solves for total pressure P_t* and parallel velocity v_parallel*
     * at the particle interface. Uses relativistic wave speeds.
     *
     * @param rho_L, rho_R  Rest-frame densities
     * @param Pt_L, Pt_R    Total pressures (gas + magnetic)
     * @param v_L, v_R      Parallel velocities
     * @param c_L, c_R      Fast magnetosonic speeds
     * @param H_L, H_R      Specific enthalpies
     * @param gamma_L, gamma_R  Lorentz factors
     * @param Pt_star       [output] Interface total pressure
     * @param v_star        [output] Interface parallel velocity
     */
    void srmhd_riemann_solver(
        real rho_L, real Pt_L, real v_L, real c_L, real H_L, real gamma_L,
        real rho_R, real Pt_R, real v_R, real c_R, real H_R, real gamma_R,
        real& Pt_star, real& v_star) const;

    /**
     * Method of Characteristics for relativistic Alfven waves
     *
     * Computes transverse velocity and magnetic field at interface
     * using characteristic upwinding with relativistic impedance.
     * Impedance Z = sqrt(rho * H * gamma^2) is computed internally.
     *
     * @param vy_L, vz_L    Left transverse velocities
     * @param By_L, Bz_L    Left transverse magnetic fields
     * @param rho_L         Left rest-frame density
     * @param H_L           Left specific enthalpy
     * @param gamma_L       Left Lorentz factor
     * @param vy_R, vz_R    Right transverse velocities
     * @param By_R, Bz_R    Right transverse magnetic fields
     * @param rho_R         Right rest-frame density
     * @param H_R           Right specific enthalpy
     * @param gamma_R       Right Lorentz factor
     * @param B_parallel    Average parallel magnetic field
     * @param vy_star, vz_star  [output] Interface transverse velocities
     * @param By_star, Bz_star  [output] Interface transverse magnetic fields
     */
    void moc_alfven_relativistic(
        real vy_L, real vz_L, real By_L, real Bz_L,
        real rho_L, real H_L, real gamma_L,
        real vy_R, real vz_R, real By_R, real Bz_R,
        real rho_R, real H_R, real gamma_R,
        real B_parallel,
        real& vy_star, real& vz_star, real& By_star, real& Bz_star) const;

    /**
     * Van Leer limiter for MUSCL reconstruction
     */
    real limiter(real dq1, real dq2) const;

public:
    /**
     * Scale SRMHD derivatives by baryon number so dS, de, dB are per-baryon.
     * Also keeps the legacy acc/dene aliases in sync.
     */
    static void normalize_srmhd_derivatives(SPHParticle & particle);

    void initialize(std::shared_ptr<SPHParameters> param) override;

    /**
     * Main calculation loop
     * Computes dS/dt (canonical momentum), de/dt (canonical energy),
     * and dB/dt (magnetic field) using Riemann solver and MOC.
     */
    void calculation(std::shared_ptr<Simulation> sim) override;
};

inline void FluidForce::normalize_srmhd_derivatives(SPHParticle & particle)
{
    if(particle.nu <= 0.0) {
        THROW_ERROR("SRMHD particle baryon number must be positive.");
    }

    const real inv_nu = 1.0 / particle.nu;
    particle.dS *= inv_nu;
    particle.de *= inv_nu;
    particle.dS_t *= inv_nu;  // Tangent momentum for 1D tests
    // NOTE: dB, acc_y_mhd, acc_z_mhd are NOT normalized by 1/ν
    // because they use GSPMHD-style formulas (ν_j/N instead of V)
    // and are already correct rates, not "generalized" quantities
    particle.acc = particle.dS;
    particle.dene = particle.de;
}

}
}
