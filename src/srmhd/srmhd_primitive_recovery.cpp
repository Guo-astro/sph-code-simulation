/**
 * SRMHD Primitive Variable Recovery
 *
 * Converts conserved variables (S, e, N, B) to primitive (v, n, P, gamma, H).
 *
 * Based on:
 * - Kitajima et al. (2025) for relativistic hydrodynamics
 * - Noble et al. (2006) for relativistic MHD primitive recovery
 *
 * Key equations:
 *   S = gamma * H * v                    (canonical momentum)
 *   e = gamma * H - P/(N*c^2)            (canonical energy)
 *   gamma = 1/sqrt(1 - v^2/c^2)          (Lorentz factor)
 *   H = 1 + u/c^2 + P/(n*c^2)            (specific enthalpy)
 *   n = N/gamma                          (rest-frame density)
 */

#include "srmhd/srmhd_primitive_recovery.hpp"
#include <cmath>
#include <algorithm>

namespace sph
{
namespace srmhd
{
namespace PrimitiveRecovery
{

real compute_sound_speed(real H, real gamma_eos, real c_speed)
{
    // c_s^2 = (gamma - 1)(H - 1) / H
    const real cs2 = (gamma_eos - 1.0) * (H - 1.0) / H;
    return std::sqrt(std::max(cs2, 0.0)) * c_speed;
}

real compute_alfven_speed(real B_sq, real rho, real H, real gamma_lor, real c_speed)
{
    // For relativistic MHD, Alfven speed is:
    // c_A^2 = b^2 / (rho*h + b^2)
    // where b^2 = B^2/gamma^2 + (v.B)^2 (comoving magnetic field squared)
    // For simplicity, assume v.B = 0 (perpendicular), so b^2 = B^2/gamma^2
    //
    // Alternatively, use lab-frame formula:
    // c_A = |B_parallel| / sqrt(rho*H*gamma^2 + B^2)

    const real rho_H_gamma2 = rho * H * gamma_lor * gamma_lor;
    const real denominator = rho_H_gamma2 + B_sq;

    if (denominator < 1e-30) return 0.0;

    // Return physical Alfven speed (same units as sound speed)
    return std::sqrt(B_sq / denominator) * c_speed;
}

real compute_fast_speed(real c_s, real c_A)
{
    // Fast magnetosonic speed (perpendicular propagation, theta=90):
    // c_f^2 = c_s^2 + c_A^2 - c_s^2 * c_A^2
    // This is the relativistic addition formula for perpendicular waves.
    //
    // For general angle theta:
    // c_f^2 = 0.5*(c_s^2 + c_A^2 + sqrt((c_s^2+c_A^2)^2 - 4*c_s^2*c_A^2*cos^2(theta)))
    //
    // We use the maximum (theta=90) for safety:
    const real c_s2 = c_s * c_s;
    const real c_A2 = c_A * c_A;
    const real c_f2 = c_s2 + c_A2 - c_s2 * c_A2;  // Relativistic addition

    return std::sqrt(std::max(c_f2, 0.0));
}

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
    // f(gamma) = (gamma^2 - 1)(X*e*gamma - 1)^2 - S^2*(X*gamma^2 - 1)^2 = 0
    auto func = [&](real gamma) -> real {
        const real term1 = (gamma * gamma - 1.0) * sqr(X * e * gamma - 1.0);
        const real term2 = S2 * sqr(X * gamma * gamma - 1.0);
        return term1 - term2;
    };

    auto dfunc = [&](real gamma) -> real {
        const real A = gamma * gamma - 1.0;
        const real B = X * e * gamma - 1.0;
        const real d_term1 = 2.0 * gamma * B * B + A * 2.0 * B * X * e;

        const real C = X * gamma * gamma - 1.0;
        const real d_term2 = S2 * 2.0 * C * 2.0 * X * gamma;

        return d_term1 - d_term2;
    };

    // Initial guess for gamma
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
            break;
        }

        const real delta = f / df;
        gamma -= delta;

        // Safety bounds
        if (gamma < 1.0) gamma = 1.0 + 1e-10;

        if (std::abs(delta) < tol) {
            break;
        }
    }

    if (gamma < 1.0) gamma = 1.0;

    return gamma;
}

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

PrimitiveVariables conserved_to_primitive(
    const vec_t & S,
    const real e,
    const real N,
    const real gamma_eos,
    const real c_speed
)
{
    return conserved_to_primitive_with_tangent(S, 0.0, e, N, gamma_eos, c_speed);
}

PrimitiveVariables conserved_to_primitive_with_tangent(
    const vec_t & S,
    const real S_t,
    const real e,
    const real N,
    const real gamma_eos,
    const real c_speed
)
{
    PrimitiveVariables prim;
    prim.alfven_speed = 0.0;
    prim.fast_speed = 0.0;

    const real S_normal_mag = std::abs(S);
    const real S_t_abs = std::abs(S_t);
    const real X = gamma_eos / (gamma_eos - 1.0);

    // Get initial estimate using total momentum magnitude
    const real S_total = std::sqrt(S_normal_mag * S_normal_mag + S_t_abs * S_t_abs);
    real gamma_lor = solve_lorentz_factor(S_total, e, N, gamma_eos, c_speed);

    // Compute H from current gamma
    real denom_H = X * gamma_lor * gamma_lor - 1.0;
    real H = (denom_H > 1e-10) ? (X * e * gamma_lor - 1.0) / denom_H : 1.0 + 1e-8;
    H = std::max(H, 1.0 + 1e-8);

    // Extract velocities from S = gamma*H*v
    real gamma_H = gamma_lor * H;
    real v_x = (gamma_H > 1e-10) ? S_normal_mag / gamma_H : 0.0;
    real v_t = (gamma_H > 1e-10) ? S_t / gamma_H : 0.0;
    real v_t2 = v_t * v_t;
    real v_x2 = v_x * v_x;

    // Iterative refinement
    for (int iter = 0; iter < 10; ++iter) {
        real v2 = v_x2 + v_t2;

        if (v2 >= 0.9999) {
            real scale = std::sqrt(0.99 / v2);
            v_x *= scale;
            v_t *= scale;
            v_x2 = v_x * v_x;
            v_t2 = v_t * v_t;
            v2 = v_x2 + v_t2;
        }

        real gamma_new = 1.0 / std::sqrt(1.0 - v2);

        denom_H = X * gamma_new * gamma_new - 1.0;
        if (denom_H < 1e-10) denom_H = 1e-10;
        real H_new = (X * e * gamma_new - 1.0) / denom_H;
        if (H_new < 1.0 + 1e-8) H_new = 1.0 + 1e-8;

        gamma_H = gamma_new * H_new;
        real v_x_new = S_normal_mag / gamma_H;
        real v_t_new = S_t / gamma_H;

        if (std::abs(gamma_new - gamma_lor) / gamma_lor < 1e-8) {
            gamma_lor = gamma_new;
            H = H_new;
            v_x = v_x_new;
            v_t = v_t_new;
            break;
        }

        gamma_lor = gamma_new;
        H = H_new;
        v_x = v_x_new;
        v_t = v_t_new;
        v_x2 = v_x * v_x;
        v_t2 = v_t * v_t;
    }

    prim.gamma_lor = gamma_lor;
    prim.vel = S / (gamma_lor * H);
    prim.vel_t = v_t;
    prim.density = N / gamma_lor;
    prim.enthalpy = H;
    prim.pressure = prim.density * (H - 1.0) * (gamma_eos - 1.0) / gamma_eos;

    // Safety: ensure subluminal velocity
    // This can happen when the iterative solver doesn't fully converge
    const real v_total2 = std::abs(prim.vel) * std::abs(prim.vel) + prim.vel_t * prim.vel_t;
    if (v_total2 >= 0.9999) {
        const real scale = std::sqrt(0.99 / v_total2);
        prim.vel = prim.vel * scale;
        prim.vel_t *= scale;
        // Recompute gamma from limited velocity
        prim.gamma_lor = 1.0 / std::sqrt(1.0 - (std::abs(prim.vel) * std::abs(prim.vel) + prim.vel_t * prim.vel_t));
    }

    // Floor pressure to prevent cold relativistic states
    // Use same floor as SR-GSPH (1e-6) for stability
    if (prim.pressure < 1e-6) {
        prim.pressure = 1e-6;
    }

    prim.sound_speed = compute_sound_speed(H, gamma_eos, c_speed);

    return prim;
}

PrimitiveVariables conserved_to_primitive_mhd(
    const vec_t & S,
    const real S_t,
    const real e,
    const real N,
    const vec3_t & B,
    const real gamma_eos,
    const real c_speed
)
{
    // For ideal MHD, the magnetic field doesn't directly enter the
    // primitive recovery (it's separately evolved). However, we need
    // to compute MHD wave speeds.
    //
    // First, recover hydrodynamic primitives:
    PrimitiveVariables prim = conserved_to_primitive_with_tangent(
        S, S_t, e, N, gamma_eos, c_speed);

    // Compute magnetic field magnitude
    const real B_sq = inner_product(B, B);

    // Compute Alfven speed
    prim.alfven_speed = compute_alfven_speed(
        B_sq, prim.density, prim.enthalpy, prim.gamma_lor, c_speed);

    // Compute fast magnetosonic speed
    prim.fast_speed = compute_fast_speed(prim.sound_speed, prim.alfven_speed);

    return prim;
}

void primitive_to_conserved(
    const PrimitiveVariables & prim,
    const real N,
    const real c_speed,
    vec_t & S_out,
    real & e_out
)
{
    // S = gamma * H * v
    S_out = prim.vel * (prim.gamma_lor * prim.enthalpy);

    // e = gamma * H - P/(N*c^2)
    e_out = prim.gamma_lor * prim.enthalpy - prim.pressure / (N * c_speed * c_speed);
}

}
}
}
