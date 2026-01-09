/**
 * @file isothermal_relaxation.hpp
 * @brief Analytical relaxation for TRUE Bonnor-Ebert spheres
 *
 * Solves the isothermal Lane-Emden equation:
 *   (1/ξ²) d/dξ(ξ² dψ/dξ) = exp(-ψ)
 *
 * Key equations:
 *   ξ = r / r_0                        [dimensionless radius]
 *   r_0 = c_s / sqrt(4πG ρ_c)          [scale length]
 *   ρ(ξ) = ρ_c × exp(-ψ(ξ))            [TRUE BE density profile]
 *   a_eq = (c_s²/r_0) × (dψ/dξ) × r̂    [equilibrium pressure gradient]
 */

#pragma once

#include <memory>
#include <vector>
#include "simulation.hpp"
#include "particle.hpp"
#include "defines.hpp"

namespace sph {

/**
 * @brief Parameters for isothermal sphere relaxation
 */
struct IsothermalRelaxationParams {
    real T_cloud;           // Temperature [K]
    real rho_center;        // Central density [code units]
    real r_0;               // Scale length r_0 = c_s/sqrt(4πGρ_c) [code units]
    real R_cloud;           // Cloud truncation radius [code units]
    real xi_s;              // Dimensionless truncation radius (R_cloud/r_0)
    real P_ext;             // External pressure [code units or K cm^-3]
    real G;                 // Gravitational constant [code units]

    // Unit conversions (code → CGS)
    real density_to_n;      // code density → n_H [cm^-3]

    // Mean molecular weight
    real mu = 1.27;         // For neutral ISM
};

/**
 * @brief Analytical relaxation for TRUE Bonnor-Ebert sphere
 *
 * Solves the isothermal Lane-Emden equation and uses:
 *   ρ(ξ) = ρ_c × exp(-ψ(ξ))     [TRUE BE density profile]
 *
 * The equilibrium pressure gradient force is:
 *   a_eq = (c_s²/r_0) × (dψ/dξ) × r̂    [outward, pointing away from center]
 *
 * Relaxation subtracts this from SPH accelerations to drive particles
 * toward positions where SPH matches analytical equilibrium.
 */
class IsothermalRelaxation {
public:
    IsothermalRelaxation();

    /**
     * @brief Initialize with parameters
     */
    void initialize(const IsothermalRelaxationParams& params);

    /**
     * @brief Compute analytical relaxation force for a particle
     * @param p Particle to compute force for
     * @return Acceleration vector (analytical pressure gradient, outward)
     */
    vec_t compute_relaxation_force(const SPHParticle& p) const;

    /**
     * @brief Apply relaxation to all particles
     * Subtracts analytical pressure gradient from SPH acceleration
     * @param sim Simulation containing particles
     * @param damping_factor Velocity damping factor (0 = no damping)
     */
    void apply_relaxation(std::shared_ptr<Simulation> sim, real damping_factor = 0.0);

    /**
     * @brief Check if relaxation is initialized
     */
    bool is_initialized() const { return m_initialized; }

    /**
     * @brief Get equilibrium density at radius r
     * ρ(r) = ρ_c × exp(-ψ(r/r_0))
     */
    real get_rho_eq(real r) const;

    /**
     * @brief Get equilibrium pressure at radius r
     * P = rho * c_s^2
     */
    real get_P_eq(real r) const;

    /**
     * @brief Get sound speed squared (constant for isothermal)
     */
    real get_c_s_squared() const { return m_c_s_sq; }

    /**
     * @brief Get scale length r_0
     */
    real get_r_0() const { return m_r_0; }

    /**
     * @brief Get cloud radius
     */
    real get_R_cloud() const { return m_params.R_cloud; }

    /**
     * @brief Update analytical profile to better match SPH density
     *
     * Measures the actual SPH density profile in radial bins and
     * fits new rho_center and r_c parameters to minimize mismatch.
     * This allows iterative convergence toward true equilibrium.
     *
     * @param sim Simulation containing particles
     * @return true if profile was updated, false if no update needed
     */
    bool update_profile_from_sph(std::shared_ptr<Simulation> sim);

    /**
     * @brief Get current central density parameter
     */
    real get_rho_center() const { return m_params.rho_center; }

private:
    IsothermalRelaxationParams m_params;
    bool m_initialized;

    // Pre-computed constants
    real m_c_s_sq;          // Sound speed squared [code units: (km/s)^2]
    real m_r_0;             // Scale length [code units]
    real m_rho_edge;        // Density at cloud edge (constant for r > R_cloud)

    // Lane-Emden solution storage
    int m_n_profile;                // Number of profile points
    std::vector<real> m_xi_arr;     // Dimensionless radius array
    std::vector<real> m_psi_arr;    // Dimensionless potential ψ(ξ)
    std::vector<real> m_dpsi_arr;   // Derivative dψ/dξ

    // Physical constants
    static constexpr real k_B_cgs = 1.380649e-16;       // erg/K
    static constexpr real m_proton_cgs = 1.6726219e-24; // g
    static constexpr real kms_cgs = 1.0e5;              // cm/s

    /**
     * @brief Interpolate ψ(ξ) from Lane-Emden solution
     */
    real interpolate_psi(real xi) const;

    /**
     * @brief Interpolate dψ/dξ from Lane-Emden solution
     */
    real interpolate_dpsi(real xi) const;
};

} // namespace sph
