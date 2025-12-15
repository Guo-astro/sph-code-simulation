/**
 * @file koyama_inutsuka_cooling.cpp
 * @brief Implementation of full ISM cooling from Koyama & Inutsuka (2000)
 * 
 * SINGLE SOURCE OF TRUTH (SSOT):
 * This C++ implementation is the authoritative source for K&I 2000 cooling
 * physics in SPH simulations. It uses tabulated data from data/ki2000_extracted/
 * which was extracted from the original paper figures.
 * 
 * Visualization:
 * - Use scripts/plot_ki2000_fig1_exact.py to reproduce original Figure 1
 * - The Python script uses the same tabulated data files
 * - All other K&I plotting scripts have been consolidated into this SSOT script
 * 
 * Data Source:
 * - data/ki2000_extracted/ contains extracted equilibrium curves from K&I 2000
 * - This C++ implementation loads and interpolates these tables
 * - Python visualization uses same data for validation
 */

#include "thermal/koyama_inutsuka_cooling.hpp"
#include "logger.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace sph {
namespace thermal {

namespace constants = koyama_inutsuka_constants;

KoyamaInutsukaCooling::KoyamaInutsukaCooling(
    real gamma, real G0, real zeta_CR)
    : m_gamma(gamma)
    , m_gamma_m1(gamma - 1.0)
    , m_energy_factor(constants::k_B / ((gamma - 1.0) * constants::m_n))
    , m_G0(G0)
    , m_zeta_CR(zeta_CR)
{
    WRITE_LOG << "======================================================================";
    WRITE_LOG << "Koyama & Inutsuka (2000) Full ISM Cooling - SSOT Implementation";
    WRITE_LOG << "======================================================================";
    WRITE_LOG << "Reference: 'Two-dimensional Simulations of Thermal Instability in the ISM'";
    WRITE_LOG << "           Koyama, H. & Inutsuka, S. (2000), ApJ 532, 980";
    WRITE_LOG << "";
    WRITE_LOG << "Single Source of Truth:";
    WRITE_LOG << "  C++ Implementation: src/thermal/koyama_inutsuka_cooling.{hpp,cpp}";
    WRITE_LOG << "  Embedded Data: include/thermal/koyama_inutsuka_data.hpp";
    WRITE_LOG << "  Python Visualization: scripts/plot_ki2000_fig1_exact.py";
    WRITE_LOG << "  NO external data files - all data embedded as constexpr arrays";
    WRITE_LOG << "";
    WRITE_LOG << "Physical processes:";
    WRITE_LOG << "  Cooling:";
    WRITE_LOG << "    - [CII] 158 micron fine-structure";
    WRITE_LOG << "    - [OI] 63/145 micron fine-structure";
    WRITE_LOG << "    - Lyman-alpha";
    WRITE_LOG << "    - H2 rovibrational";
    WRITE_LOG << "    - CO rotational";
    WRITE_LOG << "    - Gas-grain collisions";
    WRITE_LOG << "  Heating:";
    WRITE_LOG << "    - Photoelectric from dust/PAHs";
    WRITE_LOG << "    - Cosmic ray ionization";
    WRITE_LOG << "    - X-ray heating";
    WRITE_LOG << "    - H2 formation heating";
    WRITE_LOG << "";
    WRITE_LOG << "Parameters:";
    WRITE_LOG << "  FUV field: G0 = " << m_G0 << " (Habing units)";
    WRITE_LOG << "  Cosmic ray rate: zeta = " << m_zeta_CR << " s^-1";
    WRITE_LOG << "  Adiabatic index: gamma = " << m_gamma;
    WRITE_LOG << "";
    WRITE_LOG << "Embedded data tables:";
    WRITE_LOG << "  N_H = 10^19 cm^-2: " << ki2000_data::N_POINTS_N19 << " density points";
    WRITE_LOG << "  N_H = 10^20 cm^-2: " << ki2000_data::N_POINTS_N20 << " density points";
    WRITE_LOG << "";
    WRITE_LOG << "Valid density range: 10^-2 - 10^6 cm^-3";
    WRITE_LOG << "======================================================================";
}

real KoyamaInutsukaCooling::interpolate_1d(
    const real* x, const real* y, size_t n, real x0) const
{
    if (n == 0) return 0.0;
    
    // Clamp to bounds
    if (x0 <= x[0]) return y[0];
    if (x0 >= x[n-1]) return y[n-1];
    
    // Find bracketing points (log-log interpolation)
    real log_x0 = std::log10(x0);
    
    for (size_t i = 0; i < n - 1; ++i) {
        real log_x1 = std::log10(x[i]);
        real log_x2 = std::log10(x[i + 1]);
        
        if (log_x0 >= log_x1 && log_x0 <= log_x2) {
            // Log-log interpolation
            real log_y1 = std::log10(std::max(y[i], 1e-30));
            real log_y2 = std::log10(std::max(y[i + 1], 1e-30));
            
            real alpha = (log_x0 - log_x1) / (log_x2 - log_x1);
            real log_y0 = log_y1 + alpha * (log_y2 - log_y1);
            
            return std::pow(10.0, log_y0);
        }
    }
    
    return y[n-1];
}

real KoyamaInutsukaCooling::interpolate_embedded(
    real n_H, real N_H,
    const real* val_n19, size_t size_n19,
    const real* val_n20, size_t size_n20) const
{
    using namespace ki2000_data;
    
    // Use N_H to decide which table(s) to use
    constexpr real N_H_19 = 1.0e19;
    constexpr real N_H_20 = 1.0e20;
    
    if (N_H <= N_H_19) {
        // Use only N_H = 10^19 table
        return interpolate_1d(N19::n.data(), val_n19, size_n19, n_H);
    } 
    else if (N_H >= N_H_20) {
        // Use only N_H = 10^20 table
        return interpolate_1d(N20::n.data(), val_n20, size_n20, n_H);
    } 
    else {
        // Interpolate between two column densities
        real val1 = interpolate_1d(N19::n.data(), val_n19, size_n19, n_H);
        real val2 = interpolate_1d(N20::n.data(), val_n20, size_n20, n_H);
        
        // Linear interpolation in log(N_H)
        real log_N1 = std::log10(N_H_19);
        real log_N2 = std::log10(N_H_20);
        real log_N0 = std::log10(N_H);
        
        real alpha = (log_N0 - log_N1) / (log_N2 - log_N1);
        
        // Linear interpolation in log space
        real log_val1 = std::log10(std::max(val1, 1e-30));
        real log_val2 = std::log10(std::max(val2, 1e-30));
        
        return std::pow(10.0, log_val1 + alpha * (log_val2 - log_val1));
    }
}

real KoyamaInutsukaCooling::equilibrium_temperature(real n_H, real N_H) const
{
    using namespace ki2000_data;
    return interpolate_embedded(n_H, N_H,
                               N19::T_eq.data(), N_POINTS_N19,
                               N20::T_eq.data(), N_POINTS_N20);
}

real KoyamaInutsukaCooling::equilibrium_pressure(real n_H, real N_H) const
{
    using namespace ki2000_data;
    return interpolate_embedded(n_H, N_H,
                               N19::P_eq.data(), N_POINTS_N19,
                               N20::P_eq.data(), N_POINTS_N20);
}

real KoyamaInutsukaCooling::h2_fraction(real n_H, real N_H) const
{
    using namespace ki2000_data;
    return interpolate_embedded(n_H, N_H,
                               N19::x_H2.data(), N19::x_H2.size(),
                               N20::x_H2.data(), N20::x_H2.size());
}

real KoyamaInutsukaCooling::co_fraction(real n_H, real N_H) const
{
    using namespace ki2000_data;
    return interpolate_embedded(n_H, N_H,
                               N19::x_CO.data(), N19::x_CO.size(),
                               N20::x_CO.data(), N20::x_CO.size());
}

real KoyamaInutsukaCooling::electron_fraction(real n_H, real N_H) const
{
    using namespace ki2000_data;
    return interpolate_embedded(n_H, N_H,
                               N19::x_e.data(), N19::x_e.size(),
                               N20::x_e.data(), N20::x_e.size());
}

real KoyamaInutsukaCooling::cooling_timescale(real n_H, real T, real N_H) const
{
    // Get equilibrium temperature
    const real T_eq = equilibrium_temperature(n_H, N_H);
    
    // Estimate cooling rate from distance to equilibrium
    // Lambda - Gamma ~ (T - T_eq) / t_cool
    // t_cool ~ k_B T / (n |T - T_eq| * typical_rate)
    
    const real dT = std::abs(T - T_eq);
    
    if (dT < 1e-3 * T_eq) {
        return 1e30;  // Near equilibrium
    }
    
    // Typical cooling rate ~ 1e-26 erg cm^3 s^-1
    const real typical_rate = 1e-26;
    
    return constants::k_B * T / (n_H * typical_rate * m_gamma_m1);
}

real KoyamaInutsukaCooling::cooling_rate_sph(
    real rho, real u, real N_H, real dt,
    real n_to_cgs, real u_to_cgs, real t_to_cgs) const
{
    // Convert to CGS
    const real n_H = rho * n_to_cgs;
    const real u_cgs = u * u_to_cgs;
    const real dt_cgs = dt * t_to_cgs;
    
    // Get current temperature
    real T = temperature_from_energy(u_cgs);
    T = std::max(T_min, std::min(T, T_max));
    
    // Get equilibrium temperature
    const real T_eq = equilibrium_temperature(n_H, N_H);
    
    // Get cooling timescale
    const real t_cool = cooling_timescale(n_H, T, N_H);
    
    // If very close to equilibrium, no cooling needed
    if (std::abs(T - T_eq) < 1e-3 * T_eq) {
        return 0.0;
    }
    
    // Use implicit relaxation: u(t+dt) = u_eq + (u - u_eq) exp(-dt/t_cool)
    const real u_eq = energy_from_temperature(T_eq);
    const real exp_factor = std::exp(-dt_cgs / t_cool);
    const real u_new = u_eq + (u_cgs - u_eq) * exp_factor;
    
    // Return du/dt in code units
    const real dudt_cgs = (u_new - u_cgs) / dt_cgs;
    return dudt_cgs / u_to_cgs * t_to_cgs;
}

} // namespace thermal
} // namespace sph
