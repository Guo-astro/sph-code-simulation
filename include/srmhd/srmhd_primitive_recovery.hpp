#pragma once

#include "defines.hpp"
#include "vector_type.hpp"

namespace sph
{
namespace srmhd
{
namespace PrimitiveRecovery
{

/**
 * Primitive variable structure for SRMHD
 *
 * Contains all primitive variables recovered from conserved quantities.
 * Includes magnetic field contributions to wave speeds.
 */
struct PrimitiveVariables {
    vec_t vel;          // 3-velocity v (DIM components)
    real vel_t;         // Tangent velocity for 1D tests
    real density;       // Rest-frame baryon density n
    real pressure;      // Gas pressure P
    real gamma_lor;     // Lorentz factor gamma
    real enthalpy;      // Specific enthalpy H
    real sound_speed;   // Relativistic sound speed c_s
    real alfven_speed;  // Relativistic Alfven speed c_A
    real fast_speed;    // Relativistic fast magnetosonic speed c_f
};

/**
 * Solve quartic equation for Lorentz factor gamma using Newton-Raphson
 *
 * The quartic arises from combining:
 *   S = gamma * H * v
 *   e = gamma * H - P/(N*c^2)
 *   H = f(P, n) from EOS
 *   gamma = 1/sqrt(1 - v^2)
 *
 * Equation: (gamma^2 - 1)(X*e*gamma - 1)^2 - S^2*(X*gamma^2 - 1)^2 = 0
 * where X = gamma_eos / (gamma_eos - 1)
 *
 * @param S_mag     Magnitude of canonical momentum |S|
 * @param e         Canonical energy per baryon
 * @param N         Lab-frame baryon density
 * @param gamma_eos Adiabatic index
 * @param c_speed   Speed of light
 * @return Lorentz factor gamma
 */
real solve_lorentz_factor(
    const real S_mag,
    const real e,
    const real N,
    const real gamma_eos,
    const real c_speed
);

/**
 * Recover velocity from Lorentz factor
 * v = (X*gamma^2 - 1) / (gamma*(X*e*gamma - 1)) * S
 */
vec_t recover_velocity(
    const vec_t & S,
    const real gamma,
    const real e,
    const real gamma_eos
);

/**
 * Full conversion from conserved to primitive variables (no B field)
 *
 * Conserved: (S, e, N) -> Primitive: (v, n, P, gamma, H, c_s)
 */
PrimitiveVariables conserved_to_primitive(
    const vec_t & S,
    const real e,
    const real N,
    const real gamma_eos,
    const real c_speed
);

/**
 * Full conversion from conserved to primitive with tangent momentum
 *
 * For 1D simulations with transverse motion (Section 3.1.5 tests).
 * S_t = gamma * H * v_t is conserved.
 */
PrimitiveVariables conserved_to_primitive_with_tangent(
    const vec_t & S,
    const real S_t,
    const real e,
    const real N,
    const real gamma_eos,
    const real c_speed
);

/**
 * Full conversion from conserved to primitive with magnetic field
 *
 * For SRMHD, the magnetic field modifies the energy relation.
 * Also computes relativistic MHD wave speeds.
 *
 * @param S         Canonical momentum vector
 * @param S_t       Tangent momentum (for 1D)
 * @param e         Canonical energy
 * @param N         Lab-frame baryon density
 * @param B         Magnetic field vector (always 3D)
 * @param gamma_eos Adiabatic index
 * @param c_speed   Speed of light
 * @return PrimitiveVariables including MHD wave speeds
 */
PrimitiveVariables conserved_to_primitive_mhd(
    const vec_t & S,
    const real S_t,
    const real e,
    const real N,
    const vec3_t & B,
    const real gamma_eos,
    const real c_speed
);

/**
 * Convert primitive to conserved variables
 *
 * S = gamma * H * v
 * e = gamma * H - P/(N*c^2)
 */
void primitive_to_conserved(
    const PrimitiveVariables & prim,
    const real N,
    const real c_speed,
    vec_t & S_out,
    real & e_out
);

/**
 * Compute relativistic sound speed
 * c_s^2 = (gamma_eos - 1)(H - 1) / H
 */
real compute_sound_speed(real H, real gamma_eos, real c_speed);

/**
 * Compute relativistic Alfven speed
 * c_A^2 = B^2 / (rho * H * gamma^2 + B^2)
 */
real compute_alfven_speed(real B_sq, real rho, real H, real gamma_lor, real c_speed);

/**
 * Compute relativistic fast magnetosonic speed
 * c_f^2 = 0.5 * (c_s^2 + c_A^2 + sqrt((c_s^2 + c_A^2)^2 - 4*c_s^2*c_A^2*cos^2(theta)))
 * For perpendicular propagation (theta=90), c_f = sqrt(c_s^2 + c_A^2)
 */
real compute_fast_speed(real c_s, real c_A);

}
}
}
