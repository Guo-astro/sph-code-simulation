/**
 * @file polytropic_slab_relaxation.cpp
 * @brief 1D Polytropic Slab Relaxation Implementation
 * 
 * Implements relaxation for 1D self-gravitating polytropic slab using
 * the planar Lane-Emden solution: d²θ/dξ² = -θ^n
 */

#include "relaxation/polytropic_slab_relaxation.hpp"
#include "exception.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sph
{

// ============================================================================
// PolytropicSlabData Implementation
// ============================================================================

PolytropicSlabData::PolytropicSlabData()
    : m_xi_1(0.0)
    , m_n(0.0)
    , m_loaded(false)
{
}

void PolytropicSlabData::generate_solution(real n)
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
    
    std::cout << "PolytropicSlabData: Generated solution for n=" << n << std::endl;
    std::cout << "  ξ₁ = " << m_xi_1 << " (" << m_xi_array.size() << " points)" << std::endl;
}

void PolytropicSlabData::load_from_file(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        THROW_ERROR("Cannot open file: " + filename);
    }
    
    m_xi_array.clear();
    m_theta_array.clear();
    m_dtheta_array.clear();
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        real xi, theta, dtheta;
        if (iss >> xi >> theta >> dtheta) {
            m_xi_array.push_back(xi);
            m_theta_array.push_back(theta);
            m_dtheta_array.push_back(dtheta);
        }
    }
    
    // Extract n and xi_1 from header or calculate
    if (!m_xi_array.empty()) {
        m_xi_1 = m_xi_array.back();
        m_loaded = true;
    }
}

void PolytropicSlabData::save_to_file(const std::string& filename) const
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        THROW_ERROR("Cannot open file for writing: " + filename);
    }
    
    file << "# 1D Planar Lane-Emden Solution\n";
    file << "# n = " << m_n << "\n";
    file << "# xi_1 = " << m_xi_1 << "\n";
    file << "# xi  theta  dtheta/dxi\n";
    
    file << std::scientific << std::setprecision(12);
    for (size_t i = 0; i < m_xi_array.size(); ++i) {
        file << m_xi_array[i] << " " 
             << m_theta_array[i] << " " 
             << m_dtheta_array[i] << "\n";
    }
}

real PolytropicSlabData::get_theta(real xi) const
{
    if (!m_loaded || m_xi_array.empty()) {
        return 0.0;
    }
    
    xi = std::abs(xi);  // Symmetric about x=0
    
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

real PolytropicSlabData::dtheta_dxi(real xi) const
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
// PolytropicSlabRelaxation Implementation
// ============================================================================

PolytropicSlabRelaxation::PolytropicSlabRelaxation()
    : m_initialized(false)
{
}

void PolytropicSlabRelaxation::initialize(const PolytropicSlabRelaxationParams& params)
{
    m_params = params;
    
    // Generate Lane-Emden solution for this n
    m_data.generate_solution(params.n);
    
    // Compute x_surface from α and ξ₁ if not provided
    if (m_params.x_surface <= 0) {
        m_params.x_surface = m_params.alpha_scaling * m_data.get_xi_1();
    }
    
    m_initialized = true;
    
    std::cout << "PolytropicSlabRelaxation: Initialized with:" << std::endl;
    std::cout << "  γ = " << params.gamma << ", n = " << params.n << std::endl;
    std::cout << "  ρ_c = " << params.rho_center << ", K = " << params.K << std::endl;
    std::cout << "  α = " << params.alpha_scaling << std::endl;
    std::cout << "  x_surface = " << m_params.x_surface << std::endl;
    std::cout << "  ξ₁ = " << m_data.get_xi_1() << std::endl;
}

real PolytropicSlabRelaxation::get_analytic_density(real x) const
{
    if (!m_initialized) return 0.0;
    
    real xi = x / m_params.alpha_scaling;
    real theta = m_data.get_theta(xi);
    
    if (theta <= 0) return 0.0;
    
    // ρ = ρ_c θ^n
    return m_params.rho_center * std::pow(theta, m_params.n);
}

real PolytropicSlabRelaxation::get_analytic_pressure(real x) const
{
    real rho = get_analytic_density(x);
    if (rho <= 0) return 0.0;
    
    // P = K ρ^γ
    return m_params.K * std::pow(rho, m_params.gamma);
}

vec_t PolytropicSlabRelaxation::compute_relaxation_force(const SPHParticle& p) const
{
    vec_t force(0.0);
    
    if (!m_initialized) {
        THROW_ERROR("PolytropicSlabRelaxation not initialized");
    }
    
#if DIM != 1
    // This relaxation is only for 1D
    return force;
#endif
    
    real x = p.pos[0];
    real xi = x / m_params.alpha_scaling;
    real xi_abs = std::abs(xi);
    
    if (xi_abs >= m_data.get_xi_1()) {
        return force;  // Outside slab
    }
    
    if (xi_abs < 1e-12) {
        return force;  // At center, no force
    }
    
    // Get Lane-Emden derivatives
    real theta = m_data.get_theta(xi);
    real dtheta = m_data.dtheta_dxi(xi);  // Already handles sign
    
    if (theta <= 0) {
        return force;
    }
    
    // Compute analytical pressure gradient acceleration
    // a = -(1/ρ) dP/dx
    // 
    // P = K ρ^γ
    // dP/dx = K γ ρ^(γ-1) dρ/dx
    // 
    // ρ = ρ_c θ^n
    // dρ/dx = ρ_c n θ^(n-1) dθ/dx = ρ_c n θ^(n-1) (1/α) dθ/dξ
    // 
    // Therefore:
    // a = -(K γ ρ^(γ-2) dρ/dx)
    //   = -(K γ (ρ_c θ^n)^(γ-2) ρ_c n θ^(n-1) (1/α) dθ/dξ)
    //   = -(K γ ρ_c^(γ-1) n / α) θ^(n(γ-2) + n - 1) dθ/dξ
    //   = -(K γ ρ_c^(γ-1) n / α) θ^(nγ - n - 1) dθ/dξ
    //
    // For γ = 1 + 1/n: γ - 1 = 1/n, so γ = (n+1)/n
    // nγ - n - 1 = n(n+1)/n - n - 1 = (n+1) - n - 1 = 0
    //
    // So the exponent simplifies to 0:
    // a = -(K γ ρ_c^(γ-1) n / α) dθ/dξ
    
    real gamma = m_params.gamma;
    real n = m_params.n;
    real alpha = m_params.alpha_scaling;
    real rho_c = m_params.rho_center;
    real K = m_params.K;
    
    real rho_c_pow = std::pow(rho_c, gamma - 1.0);
    real exponent = n * gamma - n - 1.0;
    real theta_factor = (std::abs(exponent) < 1e-10) ? 1.0 : std::pow(theta, exponent);
    
    real prefactor = K * gamma * rho_c_pow * n / alpha;
    real a_x = -prefactor * theta_factor * dtheta;
    
    force[0] = a_x;
    
    return force;
}

void PolytropicSlabRelaxation::apply_relaxation(std::shared_ptr<Simulation> sim, real damping_factor)
{
    if (!m_initialized) {
        THROW_ERROR("PolytropicSlabRelaxation not initialized");
    }
    
    auto& particles = sim->get_particles();
    const int num_p = sim->get_particle_num();
    
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < num_p; ++i) {
        // Compute relaxation force (analytical pressure gradient)
        vec_t relax_acc = compute_relaxation_force(particles[i]);
        
        // SUBTRACT analytical pressure gradient from SPH acceleration
        // When SPH pressure gradient matches analytical → net force = 0 → equilibrium
        particles[i].acc[0] -= relax_acc[0];
#if DIM >= 2
        particles[i].acc[1] -= relax_acc[1];
#endif
#if DIM == 3
        particles[i].acc[2] -= relax_acc[2];
#endif
    }
}

} // namespace sph
