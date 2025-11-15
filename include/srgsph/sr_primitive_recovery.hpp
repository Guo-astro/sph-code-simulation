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
    vec_t vel;          // Velocity
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
