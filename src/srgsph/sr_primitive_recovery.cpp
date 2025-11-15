#include "defines.hpp"
#include "srgsph/sr_primitive_recovery.hpp"
#include <cmath>
#include <stdexcept>

namespace sph
{
namespace srgsph
{
namespace PrimitiveRecovery
{

/**
 * Solve quartic equation for Lorentz factor γ
 * Eq. 67: (γ²-1)(Xeγ-1)² - S²(Xγ²-1)² = 0
 * where X = γ_c/(γ_c-1)
 * 
 * This is derived from:
 * S = γHv  and  e = γH - P/(Nc²)
 * with EOS: P = (γ_c-1)nu
 * 
 * Uses Newton-Raphson iteration for robustness
 */
real solve_lorentz_factor(
    const real S_mag,
    const real e,
    const real N,
    const real gamma_eos,
    const real c_speed
)
{
    // X = γ_c/(γ_c-1)
    const real X = gamma_eos / (gamma_eos - 1.0);
    const real c2 = c_speed * c_speed;
    
    // Normalized variables for better numerics
    const real S2 = S_mag * S_mag;
    
    // Initial guess for γ
    // For weak fields: γ ≈ 1 + v²/(2c²)
    // For strong fields: use safe estimate
    real gamma = 1.0 + std::sqrt(1.0 + S2 / c2) / X;
    
    // Newton-Raphson iteration
    const int max_iter = 100;
    const real tol = 1.0e-10;
    
    for(int iter = 0; iter < max_iter; ++iter) {
        const real gamma2 = gamma * gamma;
        const real gamma2_m1 = gamma2 - 1.0;
        const real Xe_gamma_m1 = X * e * gamma - 1.0;
        const real X_gamma2_m1 = X * gamma2 - 1.0;
        
        // f(γ) = (γ²-1)(Xeγ-1)² - S²(Xγ²-1)²
        const real f = gamma2_m1 * Xe_gamma_m1 * Xe_gamma_m1
                     - S2 * X_gamma2_m1 * X_gamma2_m1 / c2;
        
        // f'(γ) derivative
        const real df = 2.0 * gamma * Xe_gamma_m1 * Xe_gamma_m1
                      + 2.0 * gamma2_m1 * Xe_gamma_m1 * X * e
                      - 2.0 * S2 * X_gamma2_m1 * 2.0 * X * gamma / c2;
        
        // Newton step
        const real delta = f / df;
        gamma -= delta;
        
        // Ensure γ >= 1
        if(gamma < 1.0) {
            gamma = 1.0 + 1.0e-10;
        }
        
        // Check convergence
        if(std::abs(delta) < tol * gamma) {
            return gamma;
        }
    }
    
    // If didn't converge, throw error
    throw std::runtime_error("Lorentz factor solver did not converge");
}

/**
 * Recover velocity from Lorentz factor
 * Eq. 69: v = (Xγ²-1)/(γ(Xeγ-1)) S
 */
vec_t recover_velocity(
    const vec_t & S,
    const real gamma,
    const real e,
    const real gamma_eos
)
{
    const real X = gamma_eos / (gamma_eos - 1.0);
    const real gamma2 = gamma * gamma;
    
    // Numerator: Xγ² - 1
    const real num = X * gamma2 - 1.0;
    
    // Denominator: γ(Xeγ - 1)
    const real denom = gamma * (X * e * gamma - 1.0);
    
    // v = (num/denom) * S
    const real factor = num / denom;
    return S * factor;
}

/**
 * Full conversion from conserved to primitive variables
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
    
    // 1. Solve for Lorentz factor
    const real S_mag = std::abs(S);
    prim.gamma_lor = solve_lorentz_factor(S_mag, e, N, gamma_eos, c_speed);
    
    // 2. Recover velocity
    prim.vel = recover_velocity(S, prim.gamma_lor, e, gamma_eos);
    
    // 3. Compute rest frame density
    // N = γn  =>  n = N/γ
    prim.density = N / prim.gamma_lor;
    
    // 4. Compute enthalpy H from energy equation
    // e = γH - P/(Nc²)
    // H = (e + P/(Nc²))/γ
    // But we don't know P yet, so use momentum equation:
    // S = γHv  =>  H = |S|/(γ|v|)
    const real v_mag = std::abs(prim.vel);
    if(v_mag > 1.0e-10) {
        prim.enthalpy = S_mag / (prim.gamma_lor * v_mag);
    } else {
        // At rest: H = 1 + u/c²
        const real c2 = c_speed * c_speed;
        prim.enthalpy = 1.0 + (e - 1.0) / c2;  // Approximate
    }
    
    // 5. Compute pressure from enthalpy and EOS
    // H = 1 + u/c² + P/(nc²)
    // P = (γ_c - 1)nu
    // Combining: H = 1 + P/[(γ_c-1)nc²] + P/(nc²)
    //              = 1 + γ_c P/[(γ_c-1)nc²]
    // Solve for P:
    const real c2 = c_speed * c_speed;
    prim.pressure = (gamma_eos - 1.0) * prim.density * c2 * (prim.enthalpy - 1.0) / gamma_eos;
    
    // 6. Compute sound speed
    // c_s² = γ_c P / (ρh)  where ρ = γn, h = H (relativistic enthalpy)
    const real rho = N;  // Lab frame density
    prim.sound_speed = std::sqrt(gamma_eos * prim.pressure / (rho * prim.enthalpy));
    
    return prim;
}

/**
 * Convert primitive to conserved variables
 * For initialization and testing
 */
void primitive_to_conserved(
    const PrimitiveVariables & prim,
    const real N,
    const real c_speed,
    vec_t & S_out,
    real & e_out
)
{
    // S = γHv  (Eq. 5)
    S_out = prim.vel * (prim.gamma_lor * prim.enthalpy);
    
    // e = γH - P/(Nc²)  (Eq. 6)
    const real c2 = c_speed * c_speed;
    e_out = prim.gamma_lor * prim.enthalpy - prim.pressure / (N * c2);
}

}
}
}
