#include "srgsph/sr_primitive_recovery.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace sph
{
namespace srgsph
{
namespace PrimitiveRecovery
{

/**
 * Solve quartic equation for Lorentz factor � using Newton-Raphson
 * Based on Python's recover_primitives() in srg_sph.py (lines 132-223)
 * Equation: (��-1)(Xe�-1)� - S�(X��-1)� = 0
 * where X = �_c/(�_c-1)
 */
real solve_lorentz_factor(
    const real S_mag,
    const real e,
    const real N,
    const real gamma_eos,
    const real c_speed
)
{
    const real X = gamma_eos / (gamma_eos - 1.0);
    const real S2 = S_mag * S_mag;

    // Define quartic function and its derivative
    auto func = [&](real gamma) -> real {
        const real term1 = (gamma*gamma - 1.0) * sqr(X * e * gamma - 1.0);
        const real term2 = S2 * sqr(X * gamma*gamma - 1.0);
        return term1 - term2;
    };

    auto dfunc = [&](real gamma) -> real {
        // Derivative calculation
        const real term1_A = gamma*gamma - 1.0;
        const real term1_B = X * e * gamma - 1.0;
        const real d_term1 = 2.0*gamma * term1_B*term1_B +
                             term1_A * 2.0 * term1_B * X * e;

        const real term2_A = X * gamma*gamma - 1.0;
        const real d_term2 = S2 * 2.0 * term2_A * 2.0 * X * gamma;

        return d_term1 - d_term2;
    };

    // Initial guess for gamma (lines 160-172 in Python)
    real v_guess_mag = 0.0;
    if (e > 0.0) {
        v_guess_mag = S_mag / e;
        if (v_guess_mag >= 1.0) v_guess_mag = 0.99;
    }

    real gamma_guess = 1.0 / std::sqrt(1.0 - v_guess_mag * v_guess_mag);

    // Newton-Raphson iteration
    real gamma = gamma_guess;
    constexpr int max_iter = 50;
    constexpr real tol = 1e-10;

    for (int iter = 0; iter < max_iter; ++iter) {
        const real f = func(gamma);
        const real df = dfunc(gamma);

        if (std::abs(df) < 1e-15) {
            break; // Avoid division by zero
        }

        const real delta = f / df;
        gamma -= delta;

        // Safety bounds
        if (gamma < 1.0) gamma = 1.0 + 1e-10;

        if (std::abs(delta) < tol) {
            break;
        }
    }

    // Safety check
    if (gamma < 1.0) gamma = 1.0;

    return gamma;
}

/**
 * Recover velocity from Lorentz factor
 * Based on Python lines 183-193
 * Equation: v = (Xγ²-1)/(γ(Xeγ-1)) * S
 */
vec_t recover_velocity(
    const vec_t & S,
    const real gamma,
    const real e,
    const real gamma_eos
)
{
    const real X = gamma_eos / (gamma_eos - 1.0);

    const real denom = gamma * (X * e * gamma - 1.0);
    const real num = X * gamma * gamma - 1.0;

    if (std::abs(denom) < 1e-10) {
        return vec_t(0.0);
    }

    const real factor = num / denom;
    return S * factor;
}

/**
 * Full conversion from conserved to primitive variables
 * Based on Python's recover_primitives() function (lines 132-223)
 */
PrimitiveVariables conserved_to_primitive(
    const vec_t & S,
    const real e,
    const real N,
    const real gamma_eos,
    const real c_speed
)
{
    PrimitiveVariables prim;

    // Calculate S magnitude
    const real S_mag = std::abs(S);

    // Solve for Lorentz factor
    prim.gamma_lor = solve_lorentz_factor(S_mag, e, N, gamma_eos, c_speed);

    // Recover velocity
    prim.vel = recover_velocity(S, prim.gamma_lor, e, gamma_eos);

    // Recover rest-frame density (line 205)
    // N = � / Vp,  n = N / �
    prim.density = N / prim.gamma_lor;

    // Recover enthalpy (lines 207-216)
    // H(X�� - 1) = X�e - 1
    const real X = gamma_eos / (gamma_eos - 1.0);
    const real denom_H = X * prim.gamma_lor * prim.gamma_lor - 1.0;

    if (std::abs(denom_H) < 1e-10) {
        prim.enthalpy = 1.0 + 1e-8;
    } else {
        prim.enthalpy = (X * e * prim.gamma_lor - 1.0) / denom_H;
    }

    // Safety: enthalpy must be > 1
    if (prim.enthalpy < 1.0 + 1e-8) {
        prim.enthalpy = 1.0 + 1e-8;
    }

    // Recover pressure (line 218)
    // P = n(H-1)(�-1)/�
    prim.pressure = prim.density * (prim.enthalpy - 1.0) *
                    (gamma_eos - 1.0) / gamma_eos;

    // Sound speed (lines 220-222)
    // cs^2 = (gamma-1)(H-1)/H
    const real cs2 = (gamma_eos - 1.0) * (prim.enthalpy - 1.0) / prim.enthalpy;
    prim.sound_speed = std::sqrt(std::max(cs2, 0.0)) * c_speed;

    return prim;
}

/**
 * Convert primitive to conserved variables
 * Based on equations in Python (implicit in initialization)
 * Eqs. 5-6 from paper:
 * S = �Hv
 * e = �H - P/(Nc�)
 */
void primitive_to_conserved(
    const PrimitiveVariables & prim,
    const real N,
    const real c_speed,
    vec_t & S_out,
    real & e_out
)
{
    // S = �Hv
    S_out = prim.vel * (prim.gamma_lor * prim.enthalpy);

    // e = �H - P/(Nc�)
    e_out = prim.gamma_lor * prim.enthalpy - prim.pressure / (N * c_speed * c_speed);
}

} // namespace PrimitiveRecovery
} // namespace srgsph
} // namespace sph
