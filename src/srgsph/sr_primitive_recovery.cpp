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
    // Default: no tangent momentum
    return conserved_to_primitive_with_tangent(S, 0.0, e, N, gamma_eos, c_speed);
}

/**
 * Full conversion from conserved to primitive variables with tangent momentum
 * For 1D simulations with tangent velocity (Section 3.1.5)
 * 
 * IMPORTANT: In 1D tangent velocity tests, v_t is a constant of the motion.
 * We solve for γ correctly accounting for v_t being nonzero using iteration.
 */
PrimitiveVariables conserved_to_primitive_with_tangent(
    const vec_t & S,
    const real S_t,      // Used to extract v_t for the constraint
    const real e,
    const real N,
    const real gamma_eos,
    const real c_speed
)
{
    PrimitiveVariables prim;

    const real S_normal_mag = std::abs(S);
    const real S_t_abs = std::abs(S_t);
    const real X = gamma_eos / (gamma_eos - 1.0);
    const real c2 = c_speed * c_speed;
    
    // For tangent velocity problems, we need to solve the coupled system:
    // γ = 1/√(1 - v_x² - v_t²)
    // S = γ H v_x
    // S_t = γ H v_t
    // e = γH - P/(Nc²)
    //
    // Key insight: from S and S_t, we have v_x/v_t = S/S_t (when S_t ≠ 0)
    // And v_x² + v_t² = 1 - 1/γ²
    
    // Special case: when S_t is significant, use iterative approach
    // that properly accounts for tangent velocity
    
    // Get initial estimate using total momentum magnitude
    const real S_total = std::sqrt(S_normal_mag * S_normal_mag + S_t_abs * S_t_abs);
    real gamma_lor = solve_lorentz_factor(S_total, e, N, gamma_eos, c_speed);
    
    // Compute H from current γ
    real denom_H = X * gamma_lor * gamma_lor - 1.0;
    real H = (denom_H > 1e-10) ? (X * e * gamma_lor - 1.0) / denom_H : 1.0 + 1e-8;
    H = std::max(H, 1.0 + 1e-8);
    
    // Extract velocities from S = γHv and S_t = γHv_t
    real gamma_H = gamma_lor * H;
    real v_x = (gamma_H > 1e-10) ? S_normal_mag / gamma_H : 0.0;
    real v_t = (gamma_H > 1e-10) ? S_t / gamma_H : 0.0;  // Note: preserves sign
    real v_t2 = v_t * v_t;
    real v_x2 = v_x * v_x;
    
    // Iterative refinement: ensure consistency with γ constraint
    for (int iter = 0; iter < 10; ++iter) {
        // Compute the total velocity squared
        real v2 = v_x2 + v_t2;
        
        // Safety check for superluminal
        if (v2 >= 0.9999) {
            // Scale down both velocities proportionally
            real scale = std::sqrt(0.99 / v2);
            v_x *= scale;
            v_t *= scale;
            v_x2 = v_x * v_x;
            v_t2 = v_t * v_t;
            v2 = v_x2 + v_t2;
        }
        
        // Compute γ from the constraint
        real gamma_new = 1.0 / std::sqrt(1.0 - v2);
        
        // Compute H from new γ
        denom_H = X * gamma_new * gamma_new - 1.0;
        if (denom_H < 1e-10) denom_H = 1e-10;
        real H_new = (X * e * gamma_new - 1.0) / denom_H;
        if (H_new < 1.0 + 1e-8) H_new = 1.0 + 1e-8;
        
        // Re-extract velocities with new γ and H
        gamma_H = gamma_new * H_new;
        real v_x_new = S_normal_mag / gamma_H;
        real v_t_new = S_t / gamma_H;
        
        // Check convergence
        if (std::abs(gamma_new - gamma_lor) / gamma_lor < 1e-8 &&
            std::abs(v_t_new - v_t) < 1e-8) {
            gamma_lor = gamma_new;
            H = H_new;
            v_x = v_x_new;
            v_t = v_t_new;
            break;
        }
        
        // Update for next iteration
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
    
    // Recover rest-frame density
    prim.density = N / prim.gamma_lor;
    
    // Store enthalpy
    prim.enthalpy = H;

    // Safety: enthalpy must be > 1
    if (prim.enthalpy < 1.0 + 1e-8) {
        prim.enthalpy = 1.0 + 1e-8;
    }

    // Recover pressure: P = n(H-1)(γ-1)/γ
    prim.pressure = prim.density * (prim.enthalpy - 1.0) *
                    (gamma_eos - 1.0) / gamma_eos;

    // Floor pressure to prevent unphysical near-vacuum states
    if (prim.pressure < 1e-6) {
        prim.pressure = 1e-6;
    }

    // Sound speed: cs² = (γ-1)(H-1)/H
    const real cs2 = (gamma_eos - 1.0) * (prim.enthalpy - 1.0) / prim.enthalpy;
    prim.sound_speed = std::sqrt(std::max(cs2, 0.0)) * c_speed;

    return prim;
}

/**
 * Full conversion with FIXED tangent velocity
 * For 1D simulations where v_t is known to be constant (Section 3.1.5)
 * 
 * The algorithm:
 * 1. We know v_t is fixed, so v^2 = v_x^2 + v_t^2 and gamma = 1/sqrt(1 - v^2)
 * 2. From S = gamma*H*v_x, we get v_x = S/(gamma*H)
 * 3. From e = gamma*H - P/(Nc^2) and H = 1 + u + P/(nc^2), we solve for gamma iteratively
 * 4. The constraint gamma = 1/sqrt(1 - v_x^2 - v_t^2) must be satisfied
 */
PrimitiveVariables conserved_to_primitive_fixed_vt(
    const vec_t & S,
    const real v_t_fixed,
    const real e,
    const real N,
    const real gamma_eos,
    const real c_speed
)
{
    PrimitiveVariables prim;
    
    const real S_normal_mag = std::abs(S);
    const real X = gamma_eos / (gamma_eos - 1.0);
    const real c2 = c_speed * c_speed;
    const real v_t2 = v_t_fixed * v_t_fixed;
    
    // Max allowed v_x for subluminal motion with fixed v_t
    const real v_x_max_sq = std::max(0.9999 - v_t2, 0.0);
    const real v_x_max = std::sqrt(v_x_max_sq);
    
    // Initial guess: use standard solve without v_t constraint
    real gamma_lor = solve_lorentz_factor(S_normal_mag, e, N, gamma_eos, c_speed);
    
    // Iterative refinement to satisfy v_t constraint
    for (int iter = 0; iter < 20; ++iter) {
        // Compute H from current gamma
        real denom_H = X * gamma_lor * gamma_lor - 1.0;
        if (denom_H < 1e-10) denom_H = 1e-10;
        real H = (X * e * gamma_lor - 1.0) / denom_H;
        if (H < 1.0 + 1e-8) H = 1.0 + 1e-8;
        
        // Extract v_x from S = gamma*H*v_x
        real gamma_H = gamma_lor * H;
        real v_x = (gamma_H > 1e-10) ? S_normal_mag / gamma_H : 0.0;
        
        // Clamp v_x to subluminal
        if (std::abs(v_x) > v_x_max) {
            v_x = std::copysign(v_x_max * 0.99, v_x);
        }
        
        // Compute total velocity magnitude
        real v2 = v_x * v_x + v_t2;
        
        // Safety: ensure v^2 < 1
        if (v2 >= 0.9999) {
            // Scale v_x down (v_t is fixed)
            real v_x_new = std::copysign(std::sqrt(0.99 - v_t2), v_x);
            v_x = v_x_new;
            v2 = v_x * v_x + v_t2;
        }
        
        // Compute gamma from constraint
        real gamma_new = 1.0 / std::sqrt(1.0 - v2);
        
        // Check convergence
        if (std::abs(gamma_new - gamma_lor) / gamma_lor < 1e-8) {
            gamma_lor = gamma_new;
            
            // Recompute final H
            denom_H = X * gamma_lor * gamma_lor - 1.0;
            if (denom_H < 1e-10) denom_H = 1e-10;
            H = (X * e * gamma_lor - 1.0) / denom_H;
            if (H < 1.0 + 1e-8) H = 1.0 + 1e-8;
            
            prim.gamma_lor = gamma_lor;
            prim.enthalpy = H;
            prim.vel = S / (gamma_lor * H);
            prim.vel_t = v_t_fixed;
            break;
        }
        
        gamma_lor = gamma_new;
        
        // On last iteration, store values
        if (iter == 19) {
            denom_H = X * gamma_lor * gamma_lor - 1.0;
            if (denom_H < 1e-10) denom_H = 1e-10;
            H = (X * e * gamma_lor - 1.0) / denom_H;
            if (H < 1.0 + 1e-8) H = 1.0 + 1e-8;
            
            prim.gamma_lor = gamma_lor;
            prim.enthalpy = H;
            prim.vel = S / (gamma_lor * H);
            prim.vel_t = v_t_fixed;
        }
    }
    
    // Recover rest-frame density
    prim.density = N / prim.gamma_lor;
    
    // Recover pressure: P = n(H-1)(gamma-1)/gamma
    prim.pressure = prim.density * (prim.enthalpy - 1.0) *
                    (gamma_eos - 1.0) / gamma_eos;
    
    // Floor pressure to prevent unphysical near-vacuum states
    // that cause Riemann solver failures
    if (prim.pressure < 1e-6) {
        prim.pressure = 1e-6;
    }
    
    // Sound speed: cs^2 = (gamma-1)(H-1)/H
    const real cs2 = (gamma_eos - 1.0) * (prim.enthalpy - 1.0) / prim.enthalpy;
    prim.sound_speed = std::sqrt(std::max(cs2, 0.0)) * c_speed;
    
    return prim;
}

/**
 * Convert primitive to conserved variables
 * Based on equations in Python (implicit in initialization)
 * Eqs. 5-6 from paper:
 * S = γHv
 * e = γH - P/(Nc²)
 */
void primitive_to_conserved(
    const PrimitiveVariables & prim,
    const real N,
    const real c_speed,
    vec_t & S_out,
    real & e_out
)
{
    // S = γHv
    S_out = prim.vel * (prim.gamma_lor * prim.enthalpy);

    // e = γH - P/(Nc²)
    e_out = prim.gamma_lor * prim.enthalpy - prim.pressure / (N * c_speed * c_speed);
}

} // namespace PrimitiveRecovery
} // namespace srgsph
} // namespace sph
