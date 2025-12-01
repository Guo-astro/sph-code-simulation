/**
 * @file koyama_inutsuka_cooling.hpp
 * @brief Interpolation-based cooling/heating from Koyama & Inutsuka (2000) Figure 1
 * 
 * Uses pixel-perfect digitized data from PostScript files.
 * No chemistry solving - direct table lookup for thermal equilibrium curves.
 */

#pragma once

#include "../defines.hpp"
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace sph {
namespace thermal {

/**
 * @brief Complete thermal and chemical equilibrium data from Koyama & Inutsuka (2000)
 * 
 * Provides equilibrium state via interpolation of digitized Figure 1 data:
 * - Temperature T(n) and Pressure P(n) 
 * - Chemical fractions: x_e(n), x_H2(n), x_CO(n)
 * - Heating rates: Γ_PE, Γ_CR, Γ_H2
 * - Cooling rates: Λ_CII, Λ_OI, Λ_Lya, Λ_CO
 * - Timescales: t_cool, t_rec, t_ff, t_H2
 */
class KoyamaInutsukaCooling {
public:
    /**
     * @brief Construct from digitized data files
     * @param data_dir Directory containing f1{a,b,c,d}_curve_*.txt files
     * @param N_H_column Column density [cm^-2]: use 1e19 or 1e20
     */
    explicit KoyamaInutsukaCooling(const std::string& data_dir, real N_H_column = 1e19);
    
    // ========== Thermodynamic State ==========
    
    /**
     * @brief Get equilibrium temperature at given density
     * @param n_H Number density [cm^-3]
     * @return Temperature [K]
     */
    real temperature(real n_H) const;
    
    /**
     * @brief Get equilibrium pressure at given density
     * @param n_H Number density [cm^-3]
     * @return Pressure P/k_B [K cm^-3]
     */
    real pressure(real n_H) const;
    
    // ========== Chemical Fractions ==========
    
    /**
     * @brief Get electron fraction
     * @param n_H Number density [cm^-3]
     * @return x_e = n_e/n_H (dimensionless)
     */
    real electron_fraction(real n_H) const;
    
    /**
     * @brief Get H2 molecular fraction
     * @param n_H Number density [cm^-3]
     * @return x_H2 = 2*n_H2/n_H (dimensionless)
     */
    real h2_fraction(real n_H) const;
    
    /**
     * @brief Get CO molecular fraction
     * @param n_H Number density [cm^-3]
     * @return x_CO = n_CO/n_H (dimensionless)
     */
    real co_fraction(real n_H) const;
    
    // ========== Heating Rates ==========
    
    /**
     * @brief Get photoelectric heating rate
     * @param n_H Number density [cm^-3]
     * @return Γ_PE [erg s^-1 per H nucleus]
     */
    real photoelectric_heating(real n_H) const;
    
    /**
     * @brief Get cosmic ray / X-ray heating rate
     * @param n_H Number density [cm^-3]
     * @return Γ_CR [erg s^-1 per H nucleus]
     */
    real cosmic_ray_heating(real n_H) const;
    
    /**
     * @brief Get H2 formation/dissociation heating rate
     * @param n_H Number density [cm^-3]
     * @return Γ_H2 [erg s^-1 per H nucleus]
     */
    real h2_heating(real n_H) const;
    
    /**
     * @brief Get total heating rate
     * @param n_H Number density [cm^-3]
     * @return Γ_total [erg s^-1 per H nucleus]
     */
    real total_heating(real n_H) const;
    
    // ========== Cooling Rates ==========
    
    /**
     * @brief Get C II fine-structure cooling rate
     * @param n_H Number density [cm^-3]
     * @return Λ_CII [erg s^-1 per H nucleus]
     */
    real cii_cooling(real n_H) const;
    
    /**
     * @brief Get O I fine-structure cooling rate
     * @param n_H Number density [cm^-3]
     * @return Λ_OI [erg s^-1 per H nucleus]
     */
    real oi_cooling(real n_H) const;
    
    /**
     * @brief Get Ly-alpha cooling rate
     * @param n_H Number density [cm^-3]
     * @return Λ_Lya [erg s^-1 per H nucleus]
     */
    real lya_cooling(real n_H) const;
    
    /**
     * @brief Get CO rotational/vibrational cooling rate
     * @param n_H Number density [cm^-3]
     * @return Λ_CO [erg s^-1 per H nucleus]
     */
    real co_cooling(real n_H) const;
    
    /**
     * @brief Get total cooling rate
     * @param n_H Number density [cm^-3]
     * @return Λ_total [erg s^-1 per H nucleus]
     */
    real total_cooling(real n_H) const;
    
    /**
     * @brief Get net heating-cooling rate
     * @param n_H Number density [cm^-3]
     * @return Γ - Λ [erg s^-1 per H nucleus] (positive = heating, negative = cooling)
     */
    real net_heating_cooling(real n_H) const;
    
    // ========== Timescales ==========
    
    /**
     * @brief Get cooling timescale
     * @param n_H Number density [cm^-3]
     * @return t_cool [years]
     */
    real cooling_timescale(real n_H) const;
    
    /**
     * @brief Get recombination timescale
     * @param n_H Number density [cm^-3]
     * @return t_rec [years]
     */
    real recombination_timescale(real n_H) const;
    
    /**
     * @brief Get free-fall timescale
     * @param n_H Number density [cm^-3]
     * @return t_ff [years]
     */
    real freefall_timescale(real n_H) const;
    
    /**
     * @brief Get H2 formation timescale
     * @param n_H Number density [cm^-3]
     * @return t_H2 [years]
     */
    real h2_formation_timescale(real n_H) const;
    
    // ========== SPH Integration ==========
    
    /**
     * @brief Compute cooling/heating rate to drive toward equilibrium
     * @param n_H Current density [cm^-3]
     * @param T_current Current temperature [K]
     * @param t_relax Relaxation timescale [code time units]
     * @return du/dt [code energy units / code time]
     * 
     * Uses simple relaxation: du/dt = (u_eq - u) / t_relax
     */
    real cooling_rate(real n_H, real T_current, real t_relax) const;
    
    /**
     * @brief Check if density is within valid range
     */
    bool is_valid_density(real n_H) const {
        return n_H >= m_n_min && n_H <= m_n_max;
    }
    
    /**
     * @brief Get density range
     */
    void get_density_range(real& n_min, real& n_max) const {
        n_min = m_n_min;
        n_max = m_n_max;
    }

private:
    // ========== Interpolation Data (log10 space) ==========
    
    // Master density grid and thermodynamics
    std::vector<real> m_log_n;      ///< log10(density) [cm^-3] - master grid
    std::vector<real> m_log_T;      ///< log10(temperature) [K]
    std::vector<real> m_log_P;      ///< log10(pressure/k_B) [K cm^-3]
    
    // Chemical fractions
    std::vector<real> m_log_x_e;    ///< log10(electron fraction)
    std::vector<real> m_log_x_H2;   ///< log10(H2 fraction)
    std::vector<real> m_log_x_CO;   ///< log10(CO fraction)
    
    // Heating rates [erg/s/H]
    std::vector<real> m_log_Gamma_PE;   ///< log10(photoelectric heating)
    std::vector<real> m_log_Gamma_CR;   ///< log10(cosmic ray heating)
    std::vector<real> m_log_Gamma_H2;   ///< log10(H2 heating)
    
    // Cooling rates [erg/s/H]
    std::vector<real> m_log_Lambda_CII; ///< log10(C II cooling)
    std::vector<real> m_log_Lambda_OI;  ///< log10(O I cooling)
    std::vector<real> m_log_Lambda_Lya; ///< log10(Ly-alpha cooling)
    std::vector<real> m_log_Lambda_CO;  ///< log10(CO cooling)
    
    // Timescales [years]
    std::vector<real> m_log_t_cool;  ///< log10(cooling timescale)
    std::vector<real> m_log_t_rec;   ///< log10(recombination timescale)
    std::vector<real> m_log_t_ff;    ///< log10(free-fall timescale)
    std::vector<real> m_log_t_H2;    ///< log10(H2 formation timescale)
    
    real m_n_min;                    ///< Minimum density
    real m_n_max;                    ///< Maximum density
    real m_N_H_column;               ///< Column density used
    bool m_use_high_column;          ///< True if N_H >= 5e19
    
    // Physical constants (CGS)
    static constexpr real k_B_cgs = 1.380649e-16;  ///< Boltzmann constant [erg/K]
    static constexpr real m_H_cgs = 1.6737236e-24; ///< Hydrogen mass [g]
    
    /**
     * @brief Interpolate data onto master density grid
     */
    void interpolate_to_master_grid(const real* n_data, const real* val_data, 
                                    size_t n_points, std::vector<real>& log_val_out);
    
    /**
     * @brief Linear interpolation in log-log space
     */
    real interpolate_log(real log_n_query, const std::vector<real>& log_val_data) const;
};

} // namespace thermal
} // namespace sph
