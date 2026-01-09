#pragma once

#include "defines.hpp"

namespace sph
{
namespace srmhd
{
namespace riemann
{

/**
 * State structure for SRMHD Riemann problem
 *
 * Contains all quantities needed for the relativistic MHD Riemann solver.
 */
struct SRMHDState {
    real rho;       // Rest-frame density n
    real P;         // Gas pressure
    real P_t;       // Total pressure P + B_perp^2/2
    real v;         // Normal velocity component
    real c_s;       // Sound speed
    real c_A;       // Alfven speed
    real c_f;       // Fast magnetosonic speed
    real H;         // Specific enthalpy
    real gamma_lor; // Lorentz factor
    real B_par;     // Parallel magnetic field
    real B_perp_sq; // B_perp^2 = By^2 + Bz^2
};

/**
 * HLL Riemann solver for relativistic fast magnetosonic waves
 *
 * Uses relativistic wave speed estimates and acoustic impedance
 * to compute interface pressure and velocity.
 *
 * @param L         Left state
 * @param R         Right state
 * @param Pt_star   [output] Interface total pressure
 * @param v_star    [output] Interface normal velocity
 */
void hll_solver(
    const SRMHDState& L,
    const SRMHDState& R,
    real& Pt_star,
    real& v_star
);

/**
 * HLLC Riemann solver for relativistic MHD
 *
 * Three-wave approximation (left, contact, right) with
 * proper treatment of the contact discontinuity.
 *
 * @param L         Left state
 * @param R         Right state
 * @param Pt_star   [output] Interface total pressure
 * @param v_star    [output] Interface normal velocity
 */
void hllc_solver(
    const SRMHDState& L,
    const SRMHDState& R,
    real& Pt_star,
    real& v_star
);

/**
 * Compute relativistic wave speeds for HLL-type solvers
 *
 * Uses relativistic velocity addition:
 *   lambda = (v +/- c_f) / (1 +/- v*c_f)
 *
 * @param v       Normal velocity
 * @param c_f     Fast magnetosonic speed
 * @return (lambda_minus, lambda_plus) characteristic speeds
 */
void compute_wave_speeds(
    real v, real c_f,
    real& lambda_minus, real& lambda_plus
);

/**
 * Compute relativistic acoustic impedance for MHD
 *
 * Z = rho * H * gamma^2 * c_f
 *
 * Used for pressure and velocity estimates.
 */
real compute_impedance(real rho, real H, real gamma_lor, real c_f);

/**
 * Acoustic approximation for interface pressure
 *
 * P* = (Z_L * P_R + Z_R * P_L + Z_L * Z_R * (v_L - v_R)) / (Z_L + Z_R)
 */
real acoustic_pressure(
    real P_L, real P_R,
    real v_L, real v_R,
    real Z_L, real Z_R
);

/**
 * Acoustic approximation for interface velocity
 *
 * v* = (Z_L * v_L + Z_R * v_R + P_L - P_R) / (Z_L + Z_R)
 */
real acoustic_velocity(
    real P_L, real P_R,
    real v_L, real v_R,
    real Z_L, real Z_R
);

}
}
}
