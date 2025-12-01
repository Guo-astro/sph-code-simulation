/**
 * @file koyama_inutsuka_cooling.cpp
 * @brief ISM cooling/chemistry using hardcoded interpolation data
 */

#include "../include/thermal/koyama_inutsuka_cooling.hpp"
#include "../include/thermal/koyama_inutsuka_data.hpp"
#include "../include/logger.hpp"
#include <algorithm>

namespace sph {
namespace thermal {

using namespace koyama_inutsuka_data;

KoyamaInutsukaCooling::KoyamaInutsukaCooling(const std::string& /*data_dir*/, real N_H_column)
    : m_N_H_column(N_H_column)
{
    WRITE_LOG << "======================================================================";
    WRITE_LOG << "Koyama-Inutsuka (2000) ISM Cooling (Hardcoded Data)";
    WRITE_LOG << "======================================================================";
    
    // Select column density
    if (N_H_column >= 5e19) {
        WRITE_LOG << "Column density: N_H = 1e20 cm^-2 (high shielding)";
        m_use_high_column = true;
    } else {
        WRITE_LOG << "Column density: N_H = 1e19 cm^-2 (low shielding)";
        m_use_high_column = false;
    }
    
    // Set density range from temperature data
    if (m_use_high_column) {
        m_n_min = n_T_1e20[0];
        m_n_max = n_T_1e20[N_T_1e20 - 1];
        
        // Copy data to member vectors for interpolation
        m_log_n.assign(N_T_1e20, 0.0);
        m_log_T.assign(N_T_1e20, 0.0);
        for (size_t i = 0; i < N_T_1e20; ++i) {
            m_log_n[i] = std::log10(n_T_1e20[i]);
            m_log_T[i] = std::log10(val_T_1e20[i]);
        }
        
        // Interpolate pressure onto same grid
        interpolate_to_master_grid(n_P_1e20, val_P_1e20, N_P_1e20, m_log_P);
        interpolate_to_master_grid(n_xe_1e20, val_xe_1e20, N_xe_1e20, m_log_x_e);
        
    } else {
        m_n_min = n_T_1e19[0];
        m_n_max = n_T_1e19[N_T_1e19 - 1];
        
        // Copy data to member vectors
        m_log_n.assign(N_T_1e19, 0.0);
        m_log_T.assign(N_T_1e19, 0.0);
        for (size_t i = 0; i < N_T_1e19; ++i) {
            m_log_n[i] = std::log10(n_T_1e19[i]);
            m_log_T[i] = std::log10(val_T_1e19[i]);
        }
        
        // Interpolate other quantities onto same grid
        interpolate_to_master_grid(n_P_1e19, val_P_1e19, N_P_1e19, m_log_P);
        interpolate_to_master_grid(n_xe_1e19, val_xe_1e19, N_xe_1e19, m_log_x_e);
    }
    
    // Chemical fractions (same for both column densities)
    interpolate_to_master_grid(n_xH2, val_xH2, N_xH2, m_log_x_H2);
    interpolate_to_master_grid(n_xCO, val_xCO, N_xCO, m_log_x_CO);
    
    // Heating rates
    interpolate_to_master_grid(n_Gamma_PE, val_Gamma_PE, N_Gamma_PE, m_log_Gamma_PE);
    interpolate_to_master_grid(n_Gamma_CR, val_Gamma_CR, N_Gamma_CR, m_log_Gamma_CR);
    interpolate_to_master_grid(n_Gamma_H2, val_Gamma_H2, N_Gamma_H2, m_log_Gamma_H2);
    
    // Cooling rates
    interpolate_to_master_grid(n_Lambda_CII, val_Lambda_CII, N_Lambda_CII, m_log_Lambda_CII);
    interpolate_to_master_grid(n_Lambda_OI, val_Lambda_OI, N_Lambda_OI, m_log_Lambda_OI);
    interpolate_to_master_grid(n_Lambda_Lya, val_Lambda_Lya, N_Lambda_Lya, m_log_Lambda_Lya);
    interpolate_to_master_grid(n_Lambda_CO, val_Lambda_CO, N_Lambda_CO, m_log_Lambda_CO);
    
    // Timescales
    interpolate_to_master_grid(n_t_cool, val_t_cool, N_t_cool, m_log_t_cool);
    interpolate_to_master_grid(n_t_rec, val_t_rec, N_t_rec, m_log_t_rec);
    interpolate_to_master_grid(n_t_ff, val_t_ff, N_t_ff, m_log_t_ff);
    interpolate_to_master_grid(n_t_H2, val_t_H2, N_t_H2, m_log_t_H2);
    
    WRITE_LOG << "Density range: [" << m_n_min << ", " << m_n_max << "] cm^-3";
    WRITE_LOG << "Number of points: " << m_log_n.size();
    WRITE_LOG << "✓ All interpolation tables initialized";
    WRITE_LOG << "======================================================================";
}

void KoyamaInutsukaCooling::interpolate_to_master_grid(
    const real* n_data, const real* val_data, size_t n_points,
    std::vector<real>& log_val_out)
{
    log_val_out.resize(m_log_n.size());
    
    for (size_t i = 0; i < m_log_n.size(); ++i) {
        real n = std::pow(10.0, m_log_n[i]);
        
        // Find bracketing points in source data
        size_t j = 0;
        while (j < n_points - 1 && n_data[j+1] < n) {
            ++j;
        }
        
        if (j >= n_points - 1) {
            // Extrapolate at high end
            log_val_out[i] = std::log10(val_data[n_points - 1]);
        } else if (j == 0 && n < n_data[0]) {
            // Extrapolate at low end
            log_val_out[i] = std::log10(val_data[0]);
        } else {
            // Interpolate in log-log space
            real log_n_0 = std::log10(n_data[j]);
            real log_n_1 = std::log10(n_data[j+1]);
            real log_val_0 = std::log10(val_data[j]);
            real log_val_1 = std::log10(val_data[j+1]);
            
            real frac = (m_log_n[i] - log_n_0) / (log_n_1 - log_n_0);
            log_val_out[i] = log_val_0 + frac * (log_val_1 - log_val_0);
        }
    }
}

real KoyamaInutsukaCooling::interpolate_log(real log_n_query,
                                             const std::vector<real>& log_val_data) const
{
    // Clamp to valid range
    if (log_n_query <= m_log_n.front()) {
        return std::pow(10.0, log_val_data.front());
    }
    if (log_n_query >= m_log_n.back()) {
        return std::pow(10.0, log_val_data.back());
    }
    
    // Binary search for interval
    size_t left = 0;
    size_t right = m_log_n.size() - 1;
    
    while (right - left > 1) {
        size_t mid = (left + right) / 2;
        if (log_n_query < m_log_n[mid]) {
            right = mid;
        } else {
            left = mid;
        }
    }
    
    // Linear interpolation in log space
    real frac = (log_n_query - m_log_n[left]) / (m_log_n[right] - m_log_n[left]);
    real log_val = log_val_data[left] + frac * (log_val_data[right] - log_val_data[left]);
    
    return std::pow(10.0, log_val);
}

// ========== Thermodynamic State ==========

real KoyamaInutsukaCooling::temperature(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_T);
}

real KoyamaInutsukaCooling::pressure(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_P);
}

// ========== Chemical Fractions ==========

real KoyamaInutsukaCooling::electron_fraction(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_x_e);
}

real KoyamaInutsukaCooling::h2_fraction(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_x_H2);
}

real KoyamaInutsukaCooling::co_fraction(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_x_CO);
}

// ========== Heating Rates ==========

real KoyamaInutsukaCooling::photoelectric_heating(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_Gamma_PE);
}

real KoyamaInutsukaCooling::cosmic_ray_heating(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_Gamma_CR);
}

real KoyamaInutsukaCooling::h2_heating(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_Gamma_H2);
}

real KoyamaInutsukaCooling::total_heating(real n_H) const
{
    return photoelectric_heating(n_H) + cosmic_ray_heating(n_H) + h2_heating(n_H);
}

// ========== Cooling Rates ==========

real KoyamaInutsukaCooling::cii_cooling(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_Lambda_CII);
}

real KoyamaInutsukaCooling::oi_cooling(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_Lambda_OI);
}

real KoyamaInutsukaCooling::lya_cooling(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_Lambda_Lya);
}

real KoyamaInutsukaCooling::co_cooling(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_Lambda_CO);
}

real KoyamaInutsukaCooling::total_cooling(real n_H) const
{
    return cii_cooling(n_H) + oi_cooling(n_H) + lya_cooling(n_H) + co_cooling(n_H);
}

real KoyamaInutsukaCooling::net_heating_cooling(real n_H) const
{
    return total_heating(n_H) - total_cooling(n_H);
}

// ========== Timescales ==========

real KoyamaInutsukaCooling::cooling_timescale(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_t_cool);
}

real KoyamaInutsukaCooling::recombination_timescale(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_t_rec);
}

real KoyamaInutsukaCooling::freefall_timescale(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_t_ff);
}

real KoyamaInutsukaCooling::h2_formation_timescale(real n_H) const
{
    if (n_H <= 0.0) return 0.0;
    return interpolate_log(std::log10(n_H), m_log_t_H2);
}

// ========== SPH Integration ==========

real KoyamaInutsukaCooling::cooling_rate(real n_H, real T_current, real t_relax) const
{
    if (t_relax <= 0.0) return 0.0;
    
    // Get equilibrium temperature
    real T_eq = temperature(n_H);
    
    // Compute specific internal energy (assuming ideal gas)
    // u = (k_B * T) / ((gamma - 1) * mu * m_H)
    constexpr real gamma = 5.0 / 3.0;
    constexpr real mu = 1.0;  // Neutral H
    constexpr real factor = k_B_cgs / ((gamma - 1.0) * mu * m_H_cgs);
    
    real u_current = factor * T_current;
    real u_eq = factor * T_eq;
    
    // Relaxation: du/dt = (u_eq - u_current) / t_relax
    return (u_eq - u_current) / t_relax;
}

} // namespace thermal
} // namespace sph
