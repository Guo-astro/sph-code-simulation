/**
 * @file inoue_inutsuka_cooling.cpp
 * @brief Implementation of ISM Cooling from Inoue & Inutsuka (2008)
 * 
 * Reference: "Two-Fluid MHD Simulations of Converging HI Flows in the ISM"
 * Inoue, T. & Inutsuka, S. (2008)
 * 
 * Implements the simplified cooling function (Eq. 7-9) that fits the
 * detailed Koyama & Inutsuka (2000) cooling/heating rates.
 * 
 * For IMBH-cloud scattering simulations, only single-fluid HD with this
 * cooling function is needed (two-fluid effects are negligible).
 */

#include "thermal/inoue_inutsuka_cooling.hpp"
#include "logger.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace sph {
namespace thermal {

// Namespace alias for convenience
namespace constants = inoue_inutsuka_constants;

InoueInutsukaCooling::InoueInutsukaCooling(real gamma)
    : m_gamma(gamma)
    , m_gamma_m1(gamma - 1.0)
    , m_energy_factor(constants::k_B / ((gamma - 1.0) * constants::m_n))
{
    WRITE_LOG << "======================================================================";
    WRITE_LOG << "Inoue & Inutsuka (2008) ISM Cooling";
    WRITE_LOG << "======================================================================";
    WRITE_LOG << "Based on: Two-Fluid MHD Simulations of Converging HI Flows";
    WRITE_LOG << "         in the Interstellar Medium. I: Methodology and Basic Results";
    WRITE_LOG << "";
    WRITE_LOG << "Cooling function (Eq. 7-9, corrected):";
    WRITE_LOG << "  rho_n L = n_n (-Gamma + n_n Lambda)  [erg cm^-3 s^-1]";
    WRITE_LOG << "  Gamma = 2e-26 erg s^-1";
    WRITE_LOG << "  Lambda/Gamma = 10^7 exp(-114800/(T+1000)) + 0.014 sqrt(T) exp(-92/T)";
    WRITE_LOG << "  (Note: Exponent corrected from 118400 to 114800)";
    WRITE_LOG << "";
    WRITE_LOG << "Physical parameters:";
    WRITE_LOG << "  Mean particle mass: m_n = " << constants::m_n / constants::m_proton << " m_p";
    WRITE_LOG << "  Adiabatic index: gamma = " << m_gamma;
    WRITE_LOG << "";
    WRITE_LOG << "Equilibrium phases:";
    WRITE_LOG << "  WNM: n ~ 0.57 cm^-3, T ~ 6000 K, P/k ~ 3500 K cm^-3";
    WRITE_LOG << "  CNM: n ~ 30 cm^-3,   T ~ 100 K,  P/k ~ 3000 K cm^-3";
    WRITE_LOG << "======================================================================";
}

real InoueInutsukaCooling::equilibrium_temperature(real n_H) const
{
    // At equilibrium: Γ = n Λ
    // Need to solve: Γ = n Γ (Λ/Γ)(T)
    // => 1 = n (Λ/Γ)(T)
    // => (Λ/Γ)(T) = 1/n
    
    if (n_H <= 0.0) return T_max;
    
    // Initial guess based on typical ISM phases
    real T_guess;
    if (n_H < 0.1) {
        T_guess = 8000.0;  // WNM
    } else if (n_H < 10.0) {
        T_guess = 1000.0;  // Transition
    } else {
        T_guess = 100.0;   // CNM
    }
    
    return solve_equilibrium_temperature(n_H, T_guess);
}

real InoueInutsukaCooling::solve_equilibrium_temperature(
    real n_H, real T_guess, real tol, int max_iter) const
{
    // Solve: f(T) = n (Λ/Γ)(T) - 1 = 0
    // Using Newton-Raphson with bisection fallback
    
    real T = T_guess;
    real T_lo = T_min;
    real T_hi = T_max;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        const real ratio = cooling_coefficient_ratio(T);
        const real f = n_H * ratio - 1.0;
        
        // Check convergence
        if (std::abs(f) < tol) {
            return T;
        }
        
        // Numerical derivative df/dT
        const real dT = T * 1e-6;
        const real ratio_plus = cooling_coefficient_ratio(T + dT);
        const real df = n_H * (ratio_plus - ratio) / dT;
        
        // Newton step with bounds
        real T_new = T - f / df;
        
        // Bisection fallback if Newton goes out of bounds
        if (T_new < T_lo || T_new > T_hi || !std::isfinite(T_new)) {
            T_new = std::sqrt(T_lo * T_hi);  // Geometric mean
        }
        
        // Update bounds
        if (f > 0) {
            T_hi = T;  // f(T) > 0 means T too low
        } else {
            T_lo = T;  // f(T) < 0 means T too high
        }
        
        T = T_new;
    }
    
    // Return best guess if not converged
    return T;
}

real InoueInutsukaCooling::equilibrium_pressure(real n_H) const
{
    // P = n k_B T => P/k_B = n T
    const real T_eq = equilibrium_temperature(n_H);
    return n_H * T_eq;
}

bool InoueInutsukaCooling::is_thermally_unstable(real n_H, real T) const
{
    // Balbus criterion: ∂(L/T)/∂T|_P < 0
    // At constant pressure P = n k T, we have dn/dT = -n/T
    // 
    // L/T = (-Γ + n Λ) / T
    // d(L/T)/dT|_P = d(L/T)/dT + (∂(L/T)/∂n) (dn/dT)
    //             = d(L/T)/dT + (Λ/T) (-n/T)
    //             = (-Γ + n Λ)(-1/T^2) + (n/T) dΛ/dT - n Λ / T^2
    //             = Γ/T^2 - 2 n Λ/T^2 + (n/T) dΛ/dT
    
    const real Gamma = constants::Gamma_0;
    const real Lambda = cooling_coefficient(T);
    
    // Numerical derivative dΛ/dT
    const real dT = T * 1e-6;
    const real Lambda_plus = cooling_coefficient(T + dT);
    const real dLambda_dT = (Lambda_plus - Lambda) / dT;
    
    const real T2 = T * T;
    const real deriv = Gamma / T2 - 2.0 * n_H * Lambda / T2 + (n_H / T) * dLambda_dT;
    
    return deriv < 0.0;
}

real InoueInutsukaCooling::cooling_timescale(real n_H, real T) const
{
    // t_cool = k_B T / (m_n |L| (γ-1))
    // where L is the net cooling rate per H nucleus [erg s^-1]
    
    const real L = std::abs(net_cooling_rate(n_H, T));
    
    if (L < 1e-40) {
        return 1e30;  // Near equilibrium
    }
    
    return constants::k_B * T / (constants::m_n * L * m_gamma_m1);
}

real InoueInutsukaCooling::cooling_rate_sph(
    real rho, real u, real dt,
    real n_to_cgs, real u_to_cgs, real t_to_cgs) const
{
    // Convert to CGS
    const real n_H = rho * n_to_cgs;
    const real u_cgs = u * u_to_cgs;
    const real dt_cgs = dt * t_to_cgs;
    
    // Get current temperature
    real T = temperature_from_energy(u_cgs);
    T = std::max(T_min, std::min(T, T_max));
    
    // Get equilibrium temperature at this density
    const real T_eq = equilibrium_temperature(n_H);
    
    // Get cooling timescale
    const real t_cool = cooling_timescale(n_H, T);
    
    // If very close to equilibrium, no cooling needed
    if (std::abs(T - T_eq) < 1e-3 * T_eq) {
        return 0.0;
    }
    
    // Use implicit relaxation for stability
    // du/dt = (u_eq - u) / t_cool
    // For stiff cooling, use exponential approach:
    // u(t+dt) = u_eq + (u - u_eq) exp(-dt/t_cool)
    
    const real u_eq = energy_from_temperature(T_eq);
    const real exp_factor = std::exp(-dt_cgs / t_cool);
    const real u_new = u_eq + (u_cgs - u_eq) * exp_factor;
    
    // Return du/dt in code units
    const real dudt_cgs = (u_new - u_cgs) / dt_cgs;
    return dudt_cgs / u_to_cgs * t_to_cgs;
}

} // namespace thermal
} // namespace sph
