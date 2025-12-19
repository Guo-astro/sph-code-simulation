#pragma once

#include "defines.hpp"
#include "vector_type.hpp"

namespace sph
{
namespace srgsph
{

/**
 * Primitive variables recovered from conserved variables
 */
struct PrimitiveVariables {
    vec_t vel;          // Velocity (normal components)
    real vel_t = 0.0;   // Tangent velocity (for 1D simulations with tangent)
    real density;       // Rest frame density n
    real pressure;      // Pressure P
    real sound_speed;   // Sound speed c_s
    real enthalpy;      // Specific enthalpy H
    real gamma_lor;     // Lorentz factor γ
};

/**
 * Special Relativistic Primitive Variable Recovery
 * Based on Kitajima, Inutsuka, Seno (2025) Section 2.6
 * 
 * Converts conserved variables (S, e, N) to primitive variables (v, ρ, P, ...)
 * Required because:
 * - Time evolution uses conserved variables
 * - Riemann solver needs primitive variables
 */
namespace PrimitiveRecovery
{

/**
 * Solve quartic equation for Lorentz factor γ
 * Eq. 67: (γ²-1)(Xe⟨γ⟩-1)² - S²(Xγ²-1)² = 0
 * where X = γ_c/(γ_c-1)
 * 
 * @param S Canonical momentum (magnitude)
 * @param e Canonical energy
 * @param N Baryon number density
 * @param gamma_eos EOS gamma (γ_c)
 * @param c_speed Speed of light
 * @return Lorentz factor γ
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
 * Eq. 69: v = (Xγ²-1)/(γ(Xe⟨γ⟩-1)) S
 * 
 * @param S Canonical momentum vector
 * @param gamma Lorentz factor
 * @param e Canonical energy
 * @param gamma_eos EOS gamma (γ_c)
 * @return Velocity vector
 */
vec_t recover_velocity(
    const vec_t & S,
    const real gamma,
    const real e,
    const real gamma_eos
);

/**
 * Full conversion from conserved to primitive variables
 * 
 * @param S Canonical momentum vector
 * @param e Canonical energy
 * @param N Baryon number density (lab frame)
 * @param gamma_eos EOS gamma (γ_c)
 * @param c_speed Speed of light
 * @return Complete set of primitive variables
 */
PrimitiveVariables conserved_to_primitive(
    const vec_t & S,
    const real e,
    const real N,
    const real gamma_eos,
    const real c_speed
);

/**
 * Full conversion from conserved to primitive variables with CONSERVED tangent momentum
 * For 1D simulations with tangent velocity (Section 3.1.5)
 *
 * PHYSICS: In 1D SRGSPH, the tangent MOMENTUM S_t = γHv_t is conserved
 * (dS_t/dt = 0 because there's no force in the tangent direction).
 * The tangent VELOCITY v_t is recovered from: v_t = S_t / (γH).
 * This is the physically correct approach for time evolution.
 *
 * @param S Canonical momentum vector (normal component)
 * @param S_t Tangent momentum scalar (CONSERVED quantity)
 * @param e Canonical energy
 * @param N Baryon number density (lab frame)
 * @param gamma_eos EOS gamma (γ_c)
 * @param c_speed Speed of light
 * @return Complete set of primitive variables including vel_t
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
 * Full conversion with FIXED tangent velocity (FOR INITIALIZATION ONLY)
 *
 * WARNING: This function assumes v_t is constant, which is INCORRECT for time evolution.
 * Use conserved_to_primitive_with_tangent() for normal simulation steps.
 *
 * This function is useful for:
 * - Initial condition setup where v_t is specified
 * - Testing/debugging with known v_t
 *
 * @param S Canonical momentum vector (normal component)
 * @param v_t_fixed Fixed tangent velocity (user-specified)
 * @param e Canonical energy
 * @param N Baryon number density (lab frame)
 * @param gamma_eos EOS gamma (γ_c)
 * @param c_speed Speed of light
 * @return Complete set of primitive variables with vel_t = v_t_fixed
 */
PrimitiveVariables conserved_to_primitive_fixed_vt(
    const vec_t & S,
    const real v_t_fixed,
    const real e,
    const real N,
    const real gamma_eos,
    const real c_speed
);

/**
 * Convert primitive to conserved variables
 * For initialization and testing
 * 
 * Eqs. 5-6:
 * S = γHv
 * e = γH - P/(Nc²)
 * 
 * @param prim Primitive variables
 * @param N Baryon number density (lab frame)
 * @param c_speed Speed of light
 * @param S_out Output canonical momentum
 * @param e_out Output canonical energy
 */
void primitive_to_conserved(
    const PrimitiveVariables & prim,
    const real N,
    const real c_speed,
    vec_t & S_out,
    real & e_out
);

}

}
}
