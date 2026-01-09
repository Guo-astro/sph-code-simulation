/**
 * @file koyama_inutsuka_relaxation.hpp
 * @brief Analytical relaxation for pressure-truncated Bonnor-Ebert spheres
 *        using Koyama & Inutsuka (2000) barotropic EOS
 *
 * This module computes equilibrium pressure gradients from the K&I 2000
 * thermal equilibrium curve P(ρ) = n k_B T_eq(n) and subtracts them from
 * SPH accelerations to drive particles toward hydrostatic equilibrium.
 */

#pragma once

#include <memory>
#include <vector>
#include "simulation.hpp"
#include "particle.hpp"
#include "defines.hpp"
#include "thermal/koyama_inutsuka_cooling.hpp"

namespace sph {

/**
 * @brief Parameters for K&I 2000 Bonnor-Ebert relaxation
 */
struct KIRelaxationParams {
    real R_cloud;           // Cloud truncation radius [code units]
    real M_cloud;           // Total cloud mass [code units]
    real P_ext;             // External pressure [code units]
    real rho_center;        // Central density [code units]
    real N_H;               // Column density [cm^-2] for K&I tables
    real G;                 // Gravitational constant [code units]
    real gamma;             // Adiabatic index

    // Unit conversions (code → CGS)
    real density_to_n;      // code density → n_H [cm^-3]
    real pressure_to_cgs;   // code pressure → dyn/cm^2
};

/**
 * @brief Analytical relaxation using K&I 2000 barotropic EOS
 *
 * Solves the hydrostatic equilibrium ODE:
 *   dP/dr = -ρ G M(r) / r²
 * with the barotropic EOS:
 *   P(ρ) = n k_B T_eq(n)
 *
 * where T_eq(n) comes from the embedded K&I 2000 tables.
 */
class KoyamaInutsukaRelaxation {
public:
    KoyamaInutsukaRelaxation();

    /**
     * @brief Initialize with parameters and compute equilibrium profile
     */
    void initialize(const KIRelaxationParams& params);

    /**
     * @brief Compute analytical relaxation force for a particle
     * @param p Particle to compute force for
     * @return Acceleration vector (analytical pressure gradient)
     */
    vec_t compute_relaxation_force(const SPHParticle& p) const;

    /**
     * @brief Apply relaxation to all particles
     * Subtracts analytical pressure gradient from SPH acceleration
     */
    void apply_relaxation(std::shared_ptr<Simulation> sim, real damping_factor = 0.0);

    /**
     * @brief Check if relaxation is initialized
     */
    bool is_initialized() const { return m_initialized; }

    /**
     * @brief Save computed profile to file for sharing with relaxation
     * @param filename Path to output file
     */
    void save_profile_to_file(const std::string& filename) const;

    /**
     * @brief Load profile from file (computed by IC generator)
     * @param filename Path to profile file
     * @return true if loaded successfully
     */
    bool load_profile_from_file(const std::string& filename);

    /**
     * @brief Get equilibrium density at radius r
     */
    real get_rho_eq(real r) const;

    /**
     * @brief Get equilibrium pressure at radius r
     */
    real get_P_eq(real r) const;

    /**
     * @brief Get equilibrium temperature at radius r
     */
    real get_T_eq(real r) const;

    /**
     * @brief Get enclosed mass at radius r
     */
    real get_M_enclosed(real r) const;

    /**
     * @brief Get profile radius array for iteration
     */
    const std::vector<real>& get_r_table() const { return m_r_table; }

    /**
     * @brief Get cloud radius from profile truncation
     */
    real get_R_cloud() const { return m_R_actual; }

    /**
     * @brief Get total cloud mass from integration
     */
    real get_M_cloud() const { return m_M_actual; }

private:
    KIRelaxationParams m_params;
    bool m_initialized;

    // Pre-computed equilibrium profile tables
    std::vector<real> m_r_table;        // Radius grid [code units]
    std::vector<real> m_rho_table;      // Density ρ(r) [code units]
    std::vector<real> m_P_table;        // Pressure P(r) [code units]
    std::vector<real> m_M_table;        // Enclosed mass M(r) [code units]
    std::vector<real> m_T_table;        // Temperature T(r) [K]
    std::vector<real> m_drho_dr_table;  // Density gradient dρ/dr

    // Actual profile properties (from ODE integration)
    real m_R_actual;    // Truncation radius where P = P_ext
    real m_M_actual;    // Total mass at truncation

    // K&I 2000 cooling instance for T_eq(n) lookup
    thermal::KoyamaInutsukaCooling m_ki_cooling;

    // Physical constants
    static constexpr real k_B = 1.380649e-16;      // erg/K
    static constexpr real m_proton = 1.6726219e-24; // g
    static constexpr real mu = 2.33;                // Mean molecular weight
    static constexpr real m_n = mu * m_proton;      // Mean particle mass

    /**
     * @brief Compute the Bonnor-Ebert equilibrium profile
     * Integrates the hydrostatic ODE from center until P = P_ext
     */
    void compute_equilibrium_profile();

    /**
     * @brief Scale the profile to match target R_cloud and M_cloud
     * Preserves the dimensionless profile shape while adjusting
     * radius and density to match user-specified targets.
     */
    void scale_profile_to_targets();

    /**
     * @brief Get effective sound speed squared c_eff² = dP/dρ
     * Includes the log-derivative of T_eq(n)
     */
    real get_c_eff_squared(real rho) const;

    /**
     * @brief Get equilibrium temperature from K&I 2000 tables
     */
    real get_T_eq_from_rho(real rho) const;

    /**
     * @brief Interpolate tabulated value at radius r
     */
    real interpolate(const std::vector<real>& table, real r) const;

    /**
     * @brief Compute density gradient at radius r
     */
    real get_drho_dr(real r) const;
};

} // namespace sph
