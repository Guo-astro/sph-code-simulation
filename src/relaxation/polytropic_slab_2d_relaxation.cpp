/**
 * @file polytropic_slab_2d_relaxation.cpp
 * @brief 2D Polytropic Slab Relaxation Implementation
 * 
 * Implements relaxation for 2D self-gravitating polytropic slab using
 * the planar Lane-Emden solution: d²θ/dξ² = -θ^n
 * 
 * The 2D slab is uniform in x and varies in y according to the Lane-Emden solution.
 * Gravity acts only in the y-direction.
 */

#include "relaxation/polytropic_slab_2d_relaxation.hpp"
#include "exception.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sph
{

// ============================================================================
// PolytropicSlab2DData Implementation
// ============================================================================

PolytropicSlab2DData::PolytropicSlab2DData()
    : m_xi_1(0.0)
    , m_n(0.0)
    , m_loaded(false)
{
}

void PolytropicSlab2DData::generate_solution(real n)
{
    m_n = n;
    m_xi_array.clear();
    m_theta_array.clear();
    m_dtheta_array.clear();
    
    // Solve planar Lane-Emden: d²θ/dξ² = -θ^n
    // with θ(0) = 1, θ'(0) = 0
    
    const real dxi = 0.0001;  // High resolution
    const real xi_max = 5.0;  // Search for first zero
    
    real xi = 0.0;
    real theta = 1.0;
    real dtheta = 0.0;
    
    m_xi_array.push_back(xi);
    m_theta_array.push_back(theta);
    m_dtheta_array.push_back(dtheta);
    
    while (xi < xi_max && theta > 1e-10) {
        // RK4 integration
        real k1_t = dtheta;
        real k1_dt = (theta > 0) ? -std::pow(theta, n) : 0.0;
        
        real t2 = theta + 0.5 * dxi * k1_t;
        real dt2 = dtheta + 0.5 * dxi * k1_dt;
        real k2_t = dt2;
        real k2_dt = (t2 > 0) ? -std::pow(t2, n) : 0.0;
        
        real t3 = theta + 0.5 * dxi * k2_t;
        real dt3 = dtheta + 0.5 * dxi * k2_dt;
        real k3_t = dt3;
        real k3_dt = (t3 > 0) ? -std::pow(t3, n) : 0.0;
        
        real t4 = theta + dxi * k3_t;
        real dt4 = dtheta + dxi * k3_dt;
        real k4_t = dt4;
        real k4_dt = (t4 > 0) ? -std::pow(t4, n) : 0.0;
        
        theta += dxi * (k1_t + 2*k2_t + 2*k3_t + k4_t) / 6.0;
        dtheta += dxi * (k1_dt + 2*k2_dt + 2*k3_dt + k4_dt) / 6.0;
        xi += dxi;
        
        if (theta > 0) {
            m_xi_array.push_back(xi);
            m_theta_array.push_back(theta);
            m_dtheta_array.push_back(dtheta);
        }
    }
    
    // Find first zero by interpolation
    if (m_theta_array.size() >= 2) {
        size_t last = m_theta_array.size() - 1;
        real theta_prev = m_theta_array[last - 1];
        real theta_curr = m_theta_array[last];
        real xi_prev = m_xi_array[last - 1];
        real xi_curr = m_xi_array[last];
        
        // Linear interpolation to find exact zero
        if (theta_prev > 0 && theta_curr <= 0) {
            m_xi_1 = xi_prev - theta_prev * (xi_curr - xi_prev) / (theta_curr - theta_prev);
        } else {
            m_xi_1 = m_xi_array.back();
        }
    }
    
    m_loaded = true;
    
    std::cout << "PolytropicSlab2DData: Generated solution for n=" << n << std::endl;
    std::cout << "  ξ₁ = " << m_xi_1 << " (" << m_xi_array.size() << " points)" << std::endl;
}

real PolytropicSlab2DData::get_theta(real xi) const
{
    if (!m_loaded || m_xi_array.empty()) {
        return 0.0;
    }
    
    xi = std::abs(xi);  // Symmetric about y=0
    
    if (xi >= m_xi_1) {
        return 0.0;
    }
    
    if (xi <= m_xi_array.front()) {
        return m_theta_array.front();
    }
    
    // Binary search for interval
    auto it = std::lower_bound(m_xi_array.begin(), m_xi_array.end(), xi);
    if (it == m_xi_array.end()) {
        return 0.0;
    }
    
    size_t idx = std::distance(m_xi_array.begin(), it);
    if (idx == 0) {
        return m_theta_array[0];
    }
    
    // Linear interpolation
    real xi0 = m_xi_array[idx - 1];
    real xi1 = m_xi_array[idx];
    real t0 = m_theta_array[idx - 1];
    real t1 = m_theta_array[idx];
    
    real t = (xi - xi0) / (xi1 - xi0);
    return t0 + t * (t1 - t0);
}

real PolytropicSlab2DData::dtheta_dxi(real xi) const
{
    if (!m_loaded || m_xi_array.empty()) {
        return 0.0;
    }
    
    real sign = (xi >= 0) ? 1.0 : -1.0;  // Antisymmetric derivative
    xi = std::abs(xi);
    
    if (xi >= m_xi_1) {
        return 0.0;
    }
    
    if (xi <= m_xi_array.front()) {
        return 0.0;  // dθ/dξ = 0 at center
    }
    
    // Binary search
    auto it = std::lower_bound(m_xi_array.begin(), m_xi_array.end(), xi);
    if (it == m_xi_array.end()) {
        return 0.0;
    }
    
    size_t idx = std::distance(m_xi_array.begin(), it);
    if (idx == 0) {
        return 0.0;
    }
    
    // Linear interpolation
    real xi0 = m_xi_array[idx - 1];
    real xi1 = m_xi_array[idx];
    real dt0 = m_dtheta_array[idx - 1];
    real dt1 = m_dtheta_array[idx];
    
    real t = (xi - xi0) / (xi1 - xi0);
    return sign * (dt0 + t * (dt1 - dt0));
}

// ============================================================================
// PolytropicSlab2DRelaxation Implementation
// ============================================================================

PolytropicSlab2DRelaxation::PolytropicSlab2DRelaxation()
    : m_initialized(false)
{
}

void PolytropicSlab2DRelaxation::initialize(const PolytropicSlab2DRelaxationParams& params)
{
    m_params = params;
    
    // Generate Lane-Emden solution for this n
    m_data.generate_solution(params.n);
    
    // Compute y_surface from α and ξ₁ if not provided
    if (m_params.y_surface <= 0) {
        m_params.y_surface = m_params.alpha_scaling * m_data.get_xi_1();
    }
    
    m_initialized = true;
    
    std::cout << "PolytropicSlab2DRelaxation: Initialized with:" << std::endl;
    std::cout << "  γ = " << params.gamma << ", n = " << params.n << std::endl;
    std::cout << "  ρ_c = " << params.rho_center << ", K = " << params.K << std::endl;
    std::cout << "  α = " << params.alpha_scaling << std::endl;
    std::cout << "  y_surface = " << m_params.y_surface << std::endl;
    std::cout << "  L_x = " << m_params.L_x << std::endl;
    std::cout << "  ξ₁ = " << m_data.get_xi_1() << std::endl;
}

real PolytropicSlab2DRelaxation::get_analytic_density(real y) const
{
    if (!m_initialized) return 0.0;
    
    real xi = y / m_params.alpha_scaling;
    real theta = m_data.get_theta(xi);
    
    if (theta <= 0) return 0.0;
    
    // ρ = ρ_c θ^n
    return m_params.rho_center * std::pow(theta, m_params.n);
}

real PolytropicSlab2DRelaxation::get_analytic_pressure(real y) const
{
    real rho = get_analytic_density(y);
    if (rho <= 0) return 0.0;
    
    // P = K ρ^γ
    return m_params.K * std::pow(rho, m_params.gamma);
}

vec_t PolytropicSlab2DRelaxation::compute_relaxation_force(const SPHParticle& p) const
{
    vec_t force(0.0);
    
    if (!m_initialized) {
        THROW_ERROR("PolytropicSlab2DRelaxation not initialized");
    }
    
#if DIM != 2
    // This relaxation is only for 2D
    return force;
#else
    
    // Get y-position (density varies in y, uniform in x)
    real y = p.pos[1];
    real xi = y / m_params.alpha_scaling;
    real xi_abs = std::abs(xi);
    
    if (xi_abs >= m_data.get_xi_1()) {
        return force;  // Outside slab
    }
    
    if (xi_abs < 1e-12) {
        return force;  // At center, no force (dθ/dξ = 0)
    }
    
    // Get Lane-Emden derivatives
    real theta = m_data.get_theta(xi);
    real dtheta = m_data.dtheta_dxi(xi);  // Already handles sign for y < 0
    
    if (theta <= 0) {
        return force;
    }
    
    // Compute analytical pressure gradient acceleration in y-direction
    // a_y = -(1/ρ) ∂P/∂y
    // 
    // P = K ρ^γ
    // ∂P/∂y = K γ ρ^(γ-1) ∂ρ/∂y
    // 
    // ρ = ρ_c θ^n
    // ∂ρ/∂y = ρ_c n θ^(n-1) ∂θ/∂y = ρ_c n θ^(n-1) (1/α) dθ/dξ
    // 
    // Therefore:
    // a_y = -(K γ ρ^(γ-2) ∂ρ/∂y)
    //     = -(K γ (ρ_c θ^n)^(γ-2) ρ_c n θ^(n-1) (1/α) dθ/dξ)
    //     = -(K γ ρ_c^(γ-1) n / α) θ^(n(γ-2) + n - 1) dθ/dξ
    //     = -(K γ ρ_c^(γ-1) n / α) θ^(nγ - n - 1) dθ/dξ
    //
    // For γ = 1 + 1/n: γ - 1 = 1/n, so γ = (n+1)/n
    // nγ - n - 1 = n(n+1)/n - n - 1 = (n+1) - n - 1 = 0
    //
    // So the exponent simplifies to 0:
    // a_y = -(K γ ρ_c^(γ-1) n / α) dθ/dξ
    
    real gamma = m_params.gamma;
    real n = m_params.n;
    real alpha = m_params.alpha_scaling;
    real rho_c = m_params.rho_center;
    real K = m_params.K;
    
    real rho_c_pow = std::pow(rho_c, gamma - 1.0);
    real exponent = n * gamma - n - 1.0;
    real theta_factor = (std::abs(exponent) < 1e-10) ? 1.0 : std::pow(theta, exponent);
    
    real prefactor = K * gamma * rho_c_pow * n / alpha;
    real a_y = -prefactor * theta_factor * dtheta;
    
    // Force acts only in y-direction (uniform in x)
    force[0] = 0.0;  // No x-force
    force[1] = a_y;  // y-force from analytical pressure gradient
    
    return force;
#endif
}

void PolytropicSlab2DRelaxation::apply_relaxation(std::shared_ptr<Simulation> sim, real damping_factor)
{
    if (!m_initialized) {
        THROW_ERROR("PolytropicSlab2DRelaxation not initialized");
    }
    
#if DIM != 2
    return;  // Only for 2D
#else
    
    auto& particles = sim->get_particles();
    const int num_p = sim->get_particle_num();
    
    // For relaxation, we DON'T subtract analytical forces.
    // Instead, we let SPH forces + gravity determine the equilibrium.
    // The damping removes kinetic energy so particles settle.
    // The isentropic EOS reset (done in solver) maintains the correct thermodynamics.
    
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < num_p; ++i) {
        // Apply velocity damping to accelerate convergence
        // With damping, particles gradually settle toward SPH equilibrium
        particles[i].vel[0] *= damping_factor;
        particles[i].vel[1] *= damping_factor;
    }
#endif
}

} // namespace sph
