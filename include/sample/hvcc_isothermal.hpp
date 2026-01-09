/**
 * @file hvcc_isothermal.hpp
 * @brief High-Velocity Compact Cloud IC with isothermal equilibrium (SSOT)
 *
 * Creates initial conditions for HVCC-IMBH interaction (Oka 2017 scenario)
 * using K&I 2000 thermal equilibrium at T = 10K.
 *
 * Physical model:
 *   - Isothermal sphere at T = 10K (on K&I 2000 equilibrium curve)
 *   - Uniform or centrally concentrated density profile
 *   - Ghost envelope for external pressure confinement
 *   - Point-mass IMBH at configurable position
 *
 * Reference:
 *   - Ballone et al. (2018), MNRAS - Oka 2017 IMBH candidate (CO-0.40-0.22)
 *   - Koyama & Inutsuka (2000), ApJ 532, 980 - ISM thermal equilibrium
 */

#pragma once

#include "defines.hpp"
#include "vector_type.hpp"
#include <vector>

namespace sph {

class SPHParticle;

namespace sample {

/**
 * @brief Configuration for isothermal HVCC initial conditions
 *
 * Default values match Oka 2017 cloud parameters but with
 * more physical 10K temperature on K&I 2000 curve.
 */
struct HVCCIsothermalConfig {
    // Cloud physical parameters
    real M_cloud;           ///< Cloud mass [M_sun] (default: 40)
    real R_cloud;           ///< Cloud radius [pc] (default: 0.14)
    real T_cloud;           ///< Cloud temperature [K] (default: 10)
    real n_center;          ///< Central density [cm^-3] (default: 1000)
    real N_H;               ///< Column density [cm^-2] (default: 1e19)

    /// Density profile type for cloud structure
    enum class ProfileType {
        UNIFORM,            ///< Constant density (Oka 2017-like)
        CENTRAL_CONC,       ///< Centrally concentrated (r^-2 falloff)
        KI2000_BAROTROPIC   ///< Use K&I 2000 T_eq(n) barotropic profile
    };
    ProfileType profile_type;

    // Ghost envelope parameters
    bool use_envelope;      ///< Enable ghost envelope for pressure confinement
    int envelope_layers;    ///< Number of ghost layers (default: 4)
    real P_ext;             ///< External pressure P/k_B [K cm^-3] (default: 250)

    // IMBH parameters
    real M_BH;              ///< IMBH mass [M_sun] (default: 5e4)
    real b_impact;          ///< Impact parameter / initial distance [pc] (default: 5)
    real v_rel;             ///< Relative velocity toward BH [km/s] (default: 0)
    real epsilon_BH;        ///< BH softening length [pc] (default: 0.001)

    // Numerical parameters
    int N_particles;        ///< Number of cloud particles
    int N_neighbor;         ///< Target neighbor number (default: 50)

    /**
     * @brief Default constructor with Oka 2017-like values at 10K
     *
     * K&I 2000 equilibrium at T = 10K (from first principles analysis):
     *   N_H = 10^19 cm^-2: n_eq ~ 34 cm^-3, P/k_B ~ 920 K cm^-3
     *   N_H = 10^20 cm^-2: n_eq ~ 26 cm^-3, P/k_B ~ 710 K cm^-3
     */
    HVCCIsothermalConfig()
        : M_cloud(40.0)
        , R_cloud(0.14)
        , T_cloud(10.0)
        , n_center(1000.0)
        , N_H(1.0e20)                  // Use N20 for better shielding
        , profile_type(ProfileType::UNIFORM)
        , use_envelope(true)
        , envelope_layers(4)
        , P_ext(710.0)                 // K&I 2000: P/k_B at 10K for N_H=10^20
        , M_BH(5.0e4)
        , b_impact(5.0)
        , v_rel(0.0)
        , epsilon_BH(0.001)
        , N_particles(10000)
        , N_neighbor(50)
    {}
};

/**
 * @brief Find density where T_eq = T_target on K&I 2000 curve
 *
 * Uses bisection to find n_H where equilibrium_temperature(n_H, N_H) = T_target.
 * The K&I 2000 curve is monotonically decreasing with density in the cold phase.
 *
 * @param T_target Target temperature [K]
 * @param N_H Column density [cm^-2] (1e19 or 1e20)
 * @return Density n_H [cm^-3] where T_eq = T_target
 */
real find_equilibrium_density(real T_target, real N_H);

/**
 * @brief Compute isothermal sound speed at temperature T
 *
 * c_s = sqrt(k_B T / (mu m_H))
 *
 * @param T Temperature [K]
 * @param mu Mean molecular weight (default: 1.27 for molecular gas)
 * @return Sound speed [km/s]
 */
real isothermal_sound_speed(real T, real mu = 1.27);

/**
 * @brief Convert number density to code density
 *
 * rho [M_sun/pc^3] = n_H [cm^-3] / density_to_n
 * where density_to_n = (M_sun/pc^3) / (mu * m_H)
 *
 * @param n_H Number density [cm^-3]
 * @return Mass density [M_sun/pc^3]
 */
real n_to_code_density(real n_H);

/**
 * @brief Convert code density to number density
 *
 * n_H [cm^-3] = rho [M_sun/pc^3] * density_to_n
 *
 * @param rho Mass density [M_sun/pc^3]
 * @return Number density [cm^-3]
 */
real code_density_to_n(real rho);

} // namespace sample
} // namespace sph
