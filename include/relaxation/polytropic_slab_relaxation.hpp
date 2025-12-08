#pragma once

#include <memory>
#include <vector>
#include "simulation.hpp"
#include "particle.hpp"
#include "defines.hpp"

namespace sph
{

/**
 * @brief Parameters for 1D polytropic slab relaxation
 */
struct PolytropicSlabRelaxationParams
{
    real alpha_scaling;  // Scaling factor: x = alpha * xi
    real rho_center;     // Central density ρ_c
    real K;              // Polytropic constant (P = K ρ^γ)
    real x_surface;      // Physical half-width of slab
    real gamma;          // Adiabatic index
    real n;              // Polytropic index n = 1/(γ-1)
};

/**
 * @brief Stores 1D planar Lane-Emden solution data
 * 
 * Solves d²θ/dξ² = -θ^n with θ(0)=1, θ'(0)=0
 * for 1D planar geometry.
 */
class PolytropicSlabData
{
public:
    PolytropicSlabData();
    
    /**
     * @brief Generate Lane-Emden solution for given polytropic index
     * @param n Polytropic index (n = 1/(γ-1))
     */
    void generate_solution(real n);
    
    /**
     * @brief Load solution from file
     * @param filename Path to data file
     */
    void load_from_file(const std::string& filename);
    
    /**
     * @brief Save solution to file
     * @param filename Path to data file
     */
    void save_to_file(const std::string& filename) const;
    
    /**
     * @brief Get θ(ξ) by linear interpolation
     */
    real get_theta(real xi) const;
    
    /**
     * @brief Get dθ/dξ by linear interpolation
     */
    real dtheta_dxi(real xi) const;
    
    /**
     * @brief Get first zero of θ (dimensionless surface)
     */
    real get_xi_1() const { return m_xi_1; }
    
    /**
     * @brief Get polytropic index
     */
    real get_n() const { return m_n; }
    
    bool is_loaded() const { return m_loaded; }
    
private:
    std::vector<real> m_xi_array;
    std::vector<real> m_theta_array;
    std::vector<real> m_dtheta_array;
    
    real m_xi_1;    // First zero of θ
    real m_n;       // Polytropic index
    bool m_loaded;
};

/**
 * @brief Applies relaxation forces to drive particles toward 1D polytropic slab equilibrium
 * 
 * For a self-gravitating polytropic slab:
 *   - Hydrostatic equilibrium: (1/ρ) dP/dx = -g(x) where g(x) is gravitational acceleration
 *   - The Lane-Emden solution gives the equilibrium density profile
 *   - During relaxation (no gravity), we apply analytical pressure gradient forces
 *     to maintain the density structure while SPH errors are damped out
 * 
 * The relaxation force is the analytical pressure gradient:
 *   a_relax = -(1/ρ) dP/dx = -K γ ρ^(γ-2) dρ/dx
 * 
 * This balances against SPH pressure forces, driving particles to correct positions.
 */
class PolytropicSlabRelaxation
{
public:
    PolytropicSlabRelaxation();
    
    /**
     * @brief Initialize with slab parameters
     */
    void initialize(const PolytropicSlabRelaxationParams& params);
    
    /**
     * @brief Compute relaxation force for a single particle
     * @param p Particle
     * @return Acceleration vector toward equilibrium
     */
    vec_t compute_relaxation_force(const SPHParticle& p) const;
    
    /**
     * @brief Apply relaxation forces to all particles
     * @param sim Simulation object
     * @param damping_factor Velocity damping (0-1)
     */
    void apply_relaxation(std::shared_ptr<Simulation> sim, real damping_factor = 0.9);
    
    /**
     * @brief Check if initialized
     */
    bool is_initialized() const { return m_initialized; }
    
    /**
     * @brief Get analytic density at position x
     */
    real get_analytic_density(real x) const;
    
    /**
     * @brief Get analytic pressure at position x
     */
    real get_analytic_pressure(real x) const;
    
private:
    PolytropicSlabData m_data;
    PolytropicSlabRelaxationParams m_params;
    bool m_initialized;
};

} // namespace sph
