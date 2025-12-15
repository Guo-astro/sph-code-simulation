/**
 * @file koyama_inutsuka_cooling.hpp
 * @brief Full ISM Cooling following Koyama & Inutsuka (2000) - SSOT with Embedded Data
 * 
 * This is the authoritative Single Source of Truth (SSOT) C++ implementation.
 * All tabulated data is embedded in koyama_inutsuka_data.hpp - NO external files needed.
 * 
 * Reference: Koyama, H. & Inutsuka, S. (2000), ApJ 532, 980
 */

#pragma once

#include "defines.hpp"
#include "thermal/koyama_inutsuka_data.hpp"
#include <cmath>
#include <array>

namespace sph {
namespace thermal {

namespace koyama_inutsuka_constants {
    constexpr real k_B = 1.380649e-16;
    constexpr real m_proton = 1.6726219e-24;
    constexpr real m_n = 1.27 * m_proton;
    constexpr real gamma = 5.0 / 3.0;
    constexpr real G0_default = 1.7;
    constexpr real zeta_CR_default = 1.0e-17;
}

/**
 * @brief Full ISM Cooling with embedded K&I 2000 data
 */
class KoyamaInutsukaCooling {
public:
    /**
     * @brief Constructor - uses embedded data (no external files)
     */
    explicit KoyamaInutsukaCooling(
        real gamma = koyama_inutsuka_constants::gamma,
        real G0 = koyama_inutsuka_constants::G0_default,
        real zeta_CR = koyama_inutsuka_constants::zeta_CR_default
    );
    
    // ========== Equilibrium Properties ==========
    
    real equilibrium_temperature(real n_H, real N_H) const;
    real equilibrium_pressure(real n_H, real N_H) const;
    real h2_fraction(real n_H, real N_H) const;
    real co_fraction(real n_H, real N_H) const;
    real electron_fraction(real n_H, real N_H) const;
    
    // ========== Cooling Timescale ==========
    
    real cooling_timescale(real n_H, real T, real N_H) const;
    
    // ========== SPH Integration ==========
    
    real cooling_rate_sph(
        real rho, real u, real N_H, real dt,
        real n_to_cgs, real u_to_cgs, real t_to_cgs) const;
    
    // ========== Temperature/Energy Conversion ==========
    
    real temperature_from_energy(real u) const;
    real energy_from_temperature(real T) const;
    
private:
    real m_gamma;
    real m_gamma_m1;
    real m_energy_factor;
    real m_G0;
    real m_zeta_CR;
    
    static constexpr real T_min = 5.0;
    static constexpr real T_max = 1.0e8;
    
    real interpolate_embedded(
        real n_H, real N_H,
        const real* val_n19, size_t size_n19,
        const real* val_n20, size_t size_n20) const;
    
    real interpolate_1d(const real* x, const real* y, size_t n, real x0) const;
};

// ========== Inline implementations ==========

inline real KoyamaInutsukaCooling::temperature_from_energy(real u) const
{
    return u / m_energy_factor;
}

inline real KoyamaInutsukaCooling::energy_from_temperature(real T) const
{
    return T * m_energy_factor;
}

} // namespace thermal
} // namespace sph
