/**
 * @file inoue_inutsuka_cooling.hpp
 * @brief ISM Cooling following Inoue & Inutsuka (2008)
 * 
 * Implements the simplified cooling function from:
 * "Two-Fluid MHD Simulations of Converging HI Flows in the Interstellar Medium.
 *  I: Methodology and Basic Results"
 * Inoue, T. & Inutsuka, S. (2008), ApJ
 * 
 * This is a fitting function to the detailed Koyama & Inutsuka (2000) cooling rates.
 * 
 * The net cooling function is:
 *   ρ_n L = n_n (-Γ + n_n Λ) [erg cm^-3 s^-1]
 * 
 * Where Γ (heating) and Λ (cooling) are analytic functions of temperature.
 * 
 * Note: Corrected coefficient from original paper (8↔4 typo fix in exponent):
 *   Λ/Γ = 10^7 exp(-114800/(T+1000)) + 0.014 √T exp(-92/T)
 * 
 * For IMBH-cloud scattering simulations:
 *   - Two-fluid effects are NOT needed (τ_AD >> τ_dyn)
 *   - Single-fluid HD with this cooling function is sufficient
 *   - See docs/imbh-sim/theoretical_derivation_imbh_cloud_scattering.md
 */

#pragma once

#include "defines.hpp"
#include <cmath>
#include <algorithm>

namespace sph {
namespace thermal {

/**
 * @brief Physical constants in CGS units for ISM cooling
 */
namespace inoue_inutsuka_constants {
    /// Boltzmann constant [erg K^-1]
    constexpr real k_B = 1.380649e-16;
    
    /// Proton mass [g]
    constexpr real m_proton = 1.6726219e-24;
    
    /// Mean neutral particle mass [g] (91% H + 9% He, from paper Section 2.1)
    constexpr real m_n = 1.27 * m_proton;
    
    /// Heating rate constant Γ [erg s^-1] (Eq. 8)
    constexpr real Gamma_0 = 2.0e-26;
    
    /// Adiabatic index
    constexpr real gamma = 5.0 / 3.0;
    
} // namespace inoue_inutsuka_constants

/**
 * @brief ISM Cooling following Inoue & Inutsuka (2008)
 * 
 * Key physics:
 * - Simplified analytic cooling function (Eq. 7-9)
 * - Suitable for single-fluid hydrodynamics
 * 
 * Equilibrium phases:
 * - Warm Neutral Medium (WNM): n ~ 0.57 cm^-3, T ~ 6000 K
 * - Cold Neutral Medium (CNM): n ~ 30 cm^-3, T ~ 100 K
 */
class InoueInutsukaCooling {
public:
    /**
     * @brief Constructor with optional parameters
     * @param gamma Adiabatic index (default 5/3)
     */
    explicit InoueInutsukaCooling(real gamma = inoue_inutsuka_constants::gamma);
    
    // ========== Net Cooling Function (Eq. 7-9) ==========
    
    /**
     * @brief Heating rate per H nucleus
     * @return Γ [erg s^-1] constant heating rate
     */
    static constexpr real heating_rate_constant() { return inoue_inutsuka_constants::Gamma_0; }
    
    /**
     * @brief Cooling coefficient Λ/Γ as function of temperature (Eq. 9)
     * @param T Temperature [K]
     * @return Λ/Γ [cm^3]
     * 
     * Formula: Λ/Γ = 10^7 exp(-114800/(T+1000)) + 0.014 √T exp(-92/T)
     * (Corrected exponent: 118400 → 114800)
     */
    real cooling_coefficient_ratio(real T) const;
    
    /**
     * @brief Cooling coefficient Λ(T)
     * @param T Temperature [K]
     * @return Λ [erg cm^3 s^-1]
     */
    real cooling_coefficient(real T) const;
    
    /**
     * @brief Net cooling function L = -Γ + n Λ
     * @param n_H Number density [cm^-3]
     * @param T Temperature [K]
     * @return L [erg s^-1] per H nucleus (positive = cooling, negative = heating)
     */
    real net_cooling_rate(real n_H, real T) const;
    
    /**
     * @brief Net heating-cooling rate for energy equation
     * @param n_H Number density [cm^-3]
     * @param T Temperature [K]
     * @return ρ_n L = n_n (-Γ + n_n Λ) [erg cm^-3 s^-1]
     */
    real volumetric_cooling_rate(real n_H, real T) const;
    
    // ========== Thermal Equilibrium ==========
    
    /**
     * @brief Compute equilibrium temperature at given density
     * @param n_H Number density [cm^-3]
     * @return T_eq [K] where L(n, T_eq) = 0
     * 
     * Finds temperature where Γ = n Λ (thermal equilibrium)
     */
    real equilibrium_temperature(real n_H) const;
    
    /**
     * @brief Compute equilibrium pressure at given density
     * @param n_H Number density [cm^-3]
     * @return P/k_B [K cm^-3]
     */
    real equilibrium_pressure(real n_H) const;
    
    /**
     * @brief Check if point is thermally unstable (Balbus criterion)
     * @param n_H Number density [cm^-3]
     * @param T Temperature [K]
     * @return true if ∂(L/T)/∂T|_P < 0 (thermally unstable)
     */
    bool is_thermally_unstable(real n_H, real T) const;
    
    // ========== Cooling Timescale ==========
    
    /**
     * @brief Cooling/heating timescale
     * @param n_H Number density [cm^-3]
     * @param T Temperature [K]
     * @return t_cool = k_B T / (m_n |L| (γ-1)) [s]
     */
    real cooling_timescale(real n_H, real T) const;
    
    // ========== SPH Integration ==========
    
    /**
     * @brief Compute energy change rate for SPH particle
     * @param rho Density [code units]
     * @param u Specific internal energy [code units]  
     * @param dt Timestep [code units]
     * @param n_to_cgs Density conversion factor to cm^-3
     * @param u_to_cgs Energy conversion factor to erg/g
     * @param t_to_cgs Time conversion factor to seconds
     * @return du/dt [code units]
     * 
     * Uses implicit subcycling to handle stiff cooling
     */
    real cooling_rate_sph(real rho, real u, real dt,
                          real n_to_cgs, real u_to_cgs, real t_to_cgs) const;
    
    /**
     * @brief Get temperature from specific internal energy
     * @param u Specific internal energy [erg g^-1]
     * @return T [K]
     */
    real temperature_from_energy(real u) const;
    
    /**
     * @brief Get specific internal energy from temperature
     * @param T Temperature [K]
     * @return u [erg g^-1]
     */
    real energy_from_temperature(real T) const;

private:
    real m_gamma;           ///< Adiabatic index
    real m_gamma_m1;        ///< γ - 1
    real m_energy_factor;   ///< k_B / ((γ-1) m_n)
    
    /// Minimum temperature for stability [K]
    static constexpr real T_min = 10.0;
    
    /// Maximum temperature for stability [K]
    static constexpr real T_max = 1.0e8;
    
    /**
     * @brief Newton-Raphson solver for equilibrium temperature
     * @param n_H Number density [cm^-3]
     * @param T_guess Initial guess [K]
     * @param tol Relative tolerance
     * @param max_iter Maximum iterations
     * @return T_eq [K]
     */
    real solve_equilibrium_temperature(real n_H, real T_guess, 
                                        real tol = 1e-8, int max_iter = 100) const;
};

// ========== Inline implementations ==========

inline real InoueInutsukaCooling::cooling_coefficient_ratio(real T) const
{
    // Equation (9) from the paper (corrected: 118400 → 114800)
    // Λ/Γ = 10^7 exp(-114800/(T+1000)) + 0.014 √T exp(-92/T)
    
    T = std::max(T, T_min);
    
    const real term1 = 1.0e7 * std::exp(-114800.0 / (T + 1000.0));
    const real term2 = 1.4e-2 * std::sqrt(T) * std::exp(-92.0 / T);
    
    return term1 + term2;
}

inline real InoueInutsukaCooling::cooling_coefficient(real T) const
{
    return inoue_inutsuka_constants::Gamma_0 * cooling_coefficient_ratio(T);
}

inline real InoueInutsukaCooling::net_cooling_rate(real n_H, real T) const
{
    // L = (-Γ + n Λ) per H nucleus [erg s^-1]
    const real Lambda = cooling_coefficient(T);
    return (-inoue_inutsuka_constants::Gamma_0 + n_H * Lambda);
}

inline real InoueInutsukaCooling::volumetric_cooling_rate(real n_H, real T) const
{
    // ρ_n L = n_n (-Γ + n_n Λ) [erg cm^-3 s^-1]
    return n_H * net_cooling_rate(n_H, T);
}

inline real InoueInutsukaCooling::temperature_from_energy(real u) const
{
    // u = k_B T / ((γ-1) m_n)
    // T = u (γ-1) m_n / k_B
    return u / m_energy_factor;
}

inline real InoueInutsukaCooling::energy_from_temperature(real T) const
{
    // u = k_B T / ((γ-1) m_n)
    return m_energy_factor * T;
}

} // namespace thermal
} // namespace sph
