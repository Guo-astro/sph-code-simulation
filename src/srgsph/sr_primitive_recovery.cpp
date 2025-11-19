/**
 * SR-GSPH Primitive Variable Recovery
 * 
 * Solves quartic equation to recover primitive variables (γ, v, P, n) from
 * conserved variables (S, e, N) using Kitajima et al. (2025) formulation.
 * 
 * Key Equations:
 * - Conserved variables:
 *   S = γ H v                      (canonical momentum per baryon)
 *   e = γ H - P/(N c²)             (canonical energy per baryon)
 *   N = ν/V_p                      (volume-based number density, Eq. 239)
 * 
 * - Primitive variables:
 *   γ = 1/√(1 - v²/c²)             (Lorentz factor)
 *   H = 1 + u/c² + P/(n·c²)        (specific enthalpy)
 *   P = (Γ - 1) n u                (ideal gas EOS)
 *   n = N/γ                        (rest-frame density)
 * 
 * - Quartic Equation (solved for γ):
 *   (γ² - 1)(X e γ - 1)² - S²(X γ² - 1)² = 0
 *   where X = Γ/(Γ - 1)
 * 
 * Algorithm:
 * 1. solve_lorentz_factor(): Newton-Raphson iteration for γ from quartic
 * 2. recover_velocity(): Compute v = [(Xγ² - 1) / (γ(Xeγ - 1))] S
 * 3. conserved_to_primitive(): Full recovery chain γ → v → n → H → P → c_s
 * 
 * References:
 * - Kitajima et al. (2025): Variable-h volume-based formulation
 * - Takami et al. (2014): Original quartic solver for conservative GSPH
 * - Pons et al. (1999): Relativistic hydrodynamics with tangential velocities
 */
#include "srgsph/sr_primitive_recovery.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>  // For debug output

namespace sph
{
namespace srgsph
{
namespace PrimitiveRecovery
{

/**
 * Solve quartic equation for Lorentz factor γ
 * 
 * Quartic Equation (derived from S and e definitions):
 *   (γ² - 1)(X e γ - 1)² - S²(X γ² - 1)² = 0
 * where:
 *   X = Γ/(Γ - 1) = gamma_eos / (gamma_eos - 1)
 *   S = |S| = magnitude of canonical momentum
 * 
 * Physical Constraints:
 * - γ ≥ 1 (causality)
 * - γ → 1 as v → 0 (non-relativistic limit)
 * - γ → ∞ as v → c (ultra-relativistic limit)
 * 
 * Newton-Raphson Iteration:
 * - f(γ) = (γ² - 1)(Xeγ - 1)² - S²(Xγ² - 1)²
 * - f'(γ) = 2γ(Xeγ-1)² + 2(γ²-1)(Xeγ-1)Xe - 4S²Xγ(Xγ²-1)
 * - γ_{n+1} = γ_n - f(γ_n)/f'(γ_n)
 * 
 * Initial Guess:
 * - Weak-field: γ ≈ 1 + ½(S²/e²) (valid for S << e)
 * - Fallback: γ = 1.0 (rest frame)
 * 
 * Convergence:
 * - tol = 1e-10 (relative tolerance)
 * - max_iter = 100 (typically converges in < 10 iterations)
 */
real solve_lorentz_factor(
    const real S_mag,
    const real e,
    const real N,
    const real gamma_eos,
    const real c_speed)
{
    // X = γ_c / (γ_c - 1)
    const real X = gamma_eos / (gamma_eos - 1.0);
    
    // Quartic equation: (γ²-1)(Xeγ-1)² - S²(Xγ²-1)² = 0
    // Let w = γ, then:
    // (w² - 1)(Xew - 1)² - S²(Xw² - 1)² = 0
    
    // Expand and solve using Newton-Raphson
    // Initial guess: non-relativistic limit or γ ~ 1 + v²/2
    real gamma = 1.0;
    if (S_mag > 0.0 && e > 0.0)
    {
        // Better initial guess from weak-field approximation
        gamma = 1.0 + 0.5 * S_mag * S_mag / (e * e);
    }
    
    const real tol = 1e-10;
    const int max_iter = 1000;
    bool converged = false;
    int final_iter = 0;

    for (int iter = 0; iter < max_iter; ++iter)
    {
        final_iter = iter;
        const real w = gamma;
        const real w2 = w * w;

        // f(w) = (w²-1)(Xew-1)² - S²(Xw²-1)²
        const real term1 = X * e * w - 1.0;
        const real term2 = X * w2 - 1.0;
        const real f = (w2 - 1.0) * term1 * term1 - S_mag * S_mag * term2 * term2;

        // f'(w)
        const real df_dw = 2.0 * w * term1 * term1
                         + 2.0 * (w2 - 1.0) * term1 * X * e
                         - 2.0 * S_mag * S_mag * term2 * 2.0 * X * w;

        if (std::abs(df_dw) < 1e-20)
        {
            break; // Avoid division by zero
        }

        const real delta = f / df_dw;
        gamma -= delta;

        // Ensure γ >= 1
        gamma = std::max(1.0, gamma);

        if (std::abs(delta) < tol * gamma)
        {
            converged = true;
            break;
        }
    }

    // Diagnostic: warn if primitive recovery didn't converge
    static int warn_count = 0;
    if (!converged && warn_count < 10)
    {
        std::cerr << "WARNING: Primitive recovery did not converge after " << final_iter
                  << " iterations. S=" << S_mag << ", e=" << e << ", gamma=" << gamma << std::endl;
        warn_count++;
        if (warn_count == 10)
        {
            std::cerr << "(Further warnings suppressed)" << std::endl;
        }
    }

    return gamma;
}

/**
 * Recover velocity from canonical momentum
 * 
 * Formula (derived from S = γ H v and quartic solution):
 *   v = [(X γ² - 1) / (γ (X e γ - 1))] S
 * 
 * Physical Interpretation:
 * - Numerator (Xγ² - 1): Relativistic correction to momentum
 * - Denominator γ(Xeγ - 1): Enthalpy-weighted Lorentz factor
 * - Returns vectorial velocity v with same direction as S
 * 
 * Edge Cases:
 * - If denominator → 0: return v = 0 (rest frame)
 * - v is guaranteed to satisfy |v| < c if γ is correct
 */
vec_t recover_velocity(
    const vec_t &S,
    const real gamma,
    const real e,
    const real gamma_eos)
{
    // X = γ_c / (γ_c - 1)
    const real X = gamma_eos / (gamma_eos - 1.0);
    
    // v = (Xγ²-1) / (γ(Xeγ-1)) S
    const real gamma2 = gamma * gamma;
    const real numerator = X * gamma2 - 1.0;
    const real denominator = gamma * (X * e * gamma - 1.0);
    
    if (std::abs(denominator) < 1e-20)
    {
        return vec_t(0.0);
    }
    
    const real factor = numerator / denominator;
    return S * factor;
}

/**
 * Full primitive variable recovery from conserved variables
 * 
 * Input (Conserved):
 * - S: canonical momentum per baryon (vector)
 * - e: canonical energy per baryon (scalar)
 * - N: volume-based number density (from V_p, Kitajima Eq. 239)
 * 
 * Output (Primitive):
 * - γ: Lorentz factor
 * - v: velocity (vector)
 * - n: rest-frame density
 * - H: specific enthalpy
 * - P: pressure
 * - c_s: sound speed
 * 
 * Recovery Chain:
 * 1. Solve quartic for γ = γ(S, e, N)
 * 2. Recover v = v(S, γ, e)
 * 3. Compute n = N/γ (relativistic density relation)
 * 4. Compute H = (Xeγ - 1)/(Xγ² - 1) from quartic solution
 * 5. Compute u = c²(H - 1)/Γ (internal energy)
 * 6. Compute P = (Γ - 1) n u (ideal gas EOS)
 * 7. Compute c_s = √[(Γ - 1)(H - 1)/H] (relativistic sound speed)
 * 
 * Physical Consistency:
 * - Enforces γ ≥ 1, |v| < c, P ≥ 0, c_s² ≥ 0
 * - H > 1 guaranteed by quartic solution (enthalpy includes rest mass)
 */
PrimitiveVariables conserved_to_primitive(
    const vec_t &S,
    const real e,
    const real N,
    const real gamma_eos,
    const real c_speed)
{
    PrimitiveVariables prim;
    
    // Solve for Lorentz factor
    const real S_mag = std::abs(S);
    prim.gamma_lor = solve_lorentz_factor(S_mag, e, N, gamma_eos, c_speed);
    
    // Recover velocity
    prim.vel = recover_velocity(S, prim.gamma_lor, e, gamma_eos);
    const real v_mag = std::abs(prim.vel);
    
    // Rest-frame density: n = N / γ
    prim.density = N / prim.gamma_lor;
    
    // X = γ_c / (γ_c - 1)
    const real X = gamma_eos / (gamma_eos - 1.0);
    
    // From e = γH - P/(Nc²), we get:
    // H = (Xeγ - 1) / (Xγ² - 1)  [from the quartic derivation]
    const real gamma2 = prim.gamma_lor * prim.gamma_lor;
    prim.enthalpy = (X * e * prim.gamma_lor - 1.0) / (X * gamma2 - 1.0);
    
    // From H = 1 + u/c² + P/(nc²), and ideal gas P = (γ_c - 1) n u:
    // H = 1 + u/c²(1 + (γ_c-1)) = 1 + γ_c u/c²
    // So u = c²(H - 1)/γ_c
    const real u = c_speed * c_speed * (prim.enthalpy - 1.0) / gamma_eos;
    
    // P = (γ_c - 1) n u
    prim.pressure = (gamma_eos - 1.0) * prim.density * u;
    
    // Sound speed: c_s² = (γ_c - 1)(H - 1) / H
    const real cs_sq = (gamma_eos - 1.0) * (prim.enthalpy - 1.0) / prim.enthalpy;
    prim.sound_speed = std::sqrt(std::max(0.0, cs_sq));
    
    return prim;
}

/**
 * Convert primitive variables back to conserved variables
 * 
 * Inverse transformation for verification or initialization:
 * - S = γ H v    (canonical momentum)
 * - e = γ H - P/(N c²)  (canonical energy)
 * 
 * Used for:
 * - Initializing conserved variables from IC primitives
 * - Round-trip verification of primitive recovery
 * - Debugging and testing
 */
void primitive_to_conserved(
    const PrimitiveVariables &prim,
    const real N,
    const real c_speed,
    vec_t &S_out,
    real &e_out)
{
    // S = γHv
    S_out = prim.vel * (prim.gamma_lor * prim.enthalpy);
    
    // e = γH - P/(Nc²)
    e_out = prim.gamma_lor * prim.enthalpy - prim.pressure / (N * c_speed * c_speed);
}

} // namespace PrimitiveRecovery
} // namespace srgsph
} // namespace sph
