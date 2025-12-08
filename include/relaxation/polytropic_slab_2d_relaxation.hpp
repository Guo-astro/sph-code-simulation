#pragma once

#include <memory>
#include <vector>
#include "simulation.hpp"
#include "particle.hpp"
#include "defines.hpp"

namespace sph
{

/**
 * @brief Parameters for 2D polytropic slab relaxation
 */
struct PolytropicSlab2DRelaxationParams
{
    real alpha_scaling;  // Scaling factor: y = alpha * xi
    real rho_center;     // Central density ρ_c
    real K;              // Polytropic constant (P = K ρ^γ)
    real y_surface;      // Physical half-width of slab (y-direction)
    real gamma;          // Adiabatic index
    real n;              // Polytropic index n = 1/(γ-1)
    real L_x;            // Domain width in x-direction
};

/**
 * @brief Stores planar Lane-Emden solution data for 2D slab
 * 
 * Solves d²θ/dξ² = -θ^n with θ(0)=1, θ'(0)=0
 * for 1D planar geometry (gravity acts only in y-direction).
 */
class PolytropicSlab2DData
{
public:
    PolytropicSlab2DData();
    
    /**
     * @brief Generate Lane-Emden solution for given polytropic index
     * @param n Polytropic index (n = 1/(γ-1))
     */
    void generate_solution(real n);
    
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
 * @brief Applies relaxation forces to drive particles toward 2D polytropic slab equilibrium
 * 
 * For a 2D self-gravitating polytropic slab:
 *   - Density varies in y-direction only (uniform in x)
 *   - Hydrostatic equilibrium: (1/ρ) ∂P/∂y = -g(y) where g(y) is gravity in y-direction
 *   - The planar Lane-Emden solution gives the equilibrium density profile
 *   - During relaxation, we apply analytical pressure gradient forces in y-direction
 *     to maintain the density structure while SPH errors are damped out
 * 
 * The relaxation force (y-direction only):
 *   a_relax_y = -(1/ρ) ∂P/∂y = -K γ ρ^(γ-2) ∂ρ/∂y
 * 
 * This balances against SPH pressure forces, driving particles to correct positions.
 */
class PolytropicSlab2DRelaxation
{
public:
    PolytropicSlab2DRelaxation();
    
    /**
     * @brief Initialize with slab parameters
     */
    void initialize(const PolytropicSlab2DRelaxationParams& params);
    
    /**
     * @brief Compute relaxation force for a single particle
     * @param p Particle
     * @return Acceleration vector toward equilibrium (only y-component is non-zero)
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
     * @brief Get analytic density at position y
     */
    real get_analytic_density(real y) const;
    
    /**
     * @brief Get analytic pressure at position y
     */
    real get_analytic_pressure(real y) const;
    
    /**
     * @brief Get the slab surface location in y
     */
    real get_y_surface() const { return m_params.y_surface; }
    
    /**
     * @brief Get polytropic constant K
     */
    real get_K() const { return m_params.K; }
    
    /**
     * @brief Get adiabatic index gamma
     */
    real get_gamma() const { return m_params.gamma; }
    
private:
    PolytropicSlab2DData m_data;
    PolytropicSlab2DRelaxationParams m_params;
    bool m_initialized;
};

} // namespace sph
