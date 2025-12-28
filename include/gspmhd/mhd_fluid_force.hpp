#pragma once

#include <functional>
#include <memory>
#include "fluid_force.hpp"

namespace sph
{
namespace gspmhd
{

/**
 * @brief GSPMHD Fluid Force module based on Iwasaki & Inutsuka (2011)
 *
 * Implements the Godunov SPH method for MHD with:
 * - Non-linear Riemann solver with magnetic pressure for compressive waves
 * - Method of Characteristics (MOC) for Alfven waves
 * - Powell 8-wave formulation for div-B correction (tensile instability suppression)
 * - MUSCL reconstruction for 2nd order spatial accuracy
 */
class FluidForce : public sph::FluidForce {
    bool m_is_2nd_order;
    real m_gamma;
    bool m_use_powell;

    // HLL-type Riemann solver for compressive waves
    // Solves for P_t* (total pressure) and v_parallel* (parallel velocity)
    void mhd_riemann_solver(
        real rho_L, real Pt_L, real v_L, real c_L,
        real rho_R, real Pt_R, real v_R, real c_R,
        real& Pt_star, real& v_star) const;

    // Method of Characteristics for Alfven waves (3D scalar version)
    // Computes B_perp* and v_perp* from characteristics using separate y,z components
    void moc_alfven_3d(
        real vy_L, real vz_L, real By_L, real Bz_L, real rho_L,
        real vy_R, real vz_R, real By_R, real Bz_R, real rho_R,
        real B_parallel,
        real& vy_star, real& vz_star, real& By_star, real& Bz_star) const;

    // Van Leer limiter for MUSCL reconstruction
    real limiter(real dq1, real dq2) const;

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}
