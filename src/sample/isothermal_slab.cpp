#include <cmath>
#include <vector>
#include <iostream>

#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"

namespace sph
{

// 1D Isothermal Self-Gravitating Slab
// For studying diffusive instability in GSPH without grad-h correction
//
// Density profile: ρ(x) = ρ_0 * sech²(x/H)
// where H is the scale height determined by hydrostatic equilibrium
//
// In hydrostatic equilibrium:
//   dP/dx = -ρ * g(x)
//   g(x) = 2πG * Σ(x) where Σ(x) = ∫_{-∞}^{x} ρ dx
//
// For isothermal equation of state P = c_s² ρ:
//   c_s² dρ/dx = -2πG ρ Σ
//
// Solution: ρ(x) = ρ_0 sech²(x/H) with H = c_s / √(2πG ρ_0)

void Solver::make_isothermal_slab()
{
#if DIM != 1
    THROW_ERROR("Isothermal slab requires DIM == 1");
#else

    // Read parameters
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real rho_center = boost::any_cast<double>(m_sample_parameters["rho_center"]);
    const real c_s = boost::any_cast<double>(m_sample_parameters["sound_speed"]);
    const real H_input = boost::any_cast<double>(m_sample_parameters["scale_height"]);
    
    const real G = m_param->gravity.constant;
    const real gamma = m_param->physics.gamma;
    
    std::cout << "\n=== Creating 1D Isothermal Self-Gravitating Slab ===" << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "  N = " << N << " particles" << std::endl;
    std::cout << "  ρ_center = " << rho_center << std::endl;
    std::cout << "  c_s = " << c_s << std::endl;
    std::cout << "  H (input) = " << H_input << std::endl;
    std::cout << "  G = " << G << std::endl;
    std::cout << "  γ = " << gamma << std::endl;
    
    // For self-consistent hydrostatic equilibrium:
    // H = c_s / sqrt(2πG ρ_0)
    const real H_equilibrium = c_s / std::sqrt(2.0 * M_PI * G * rho_center);
    std::cout << "  H (equilibrium) = " << H_equilibrium << std::endl;
    
    // Use input H (may not be exactly in equilibrium, which is fine for testing)
    const real H = H_input;
    
    // Total mass of the slab (integrate ρ from -∞ to ∞)
    // M = ∫ ρ_0 sech²(x/H) dx = 2 ρ_0 H
    const real M_total = 2.0 * rho_center * H;
    const real mass_per_particle = M_total / N;
    
    std::cout << "  M_total = " << M_total << std::endl;
    std::cout << "  mass per particle = " << mass_per_particle << std::endl;
    
    // ========================================================================
    // Place particles with equal mass using cumulative mass function
    // M(x) = ρ_0 H [tanh(x/H) + 1] for x from -∞ to x
    // Normalized: F(x) = M(x)/M_total = 0.5 * [tanh(x/H) + 1]
    // Inverse: x = H * arctanh(2F - 1)
    // ========================================================================
    
    std::vector<SPHParticle> particles(N);
    
    for (int i = 0; i < N; ++i) {
        auto & p_i = particles[i];
        
        // Equal-mass placement: F_i = (i + 0.5) / N
        const real F_i = (i + 0.5) / N;
        
        // Invert cumulative mass function
        // x = H * arctanh(2F - 1)
        // But arctanh is only defined for |arg| < 1, so we need to be careful
        // at the boundaries. Use a slightly smaller range.
        const real arg = 2.0 * F_i - 1.0;
        const real x = H * std::atanh(arg);
        
        // Position
        p_i.pos[0] = x;
        
        // Velocity (start at rest)
        p_i.vel[0] = 0.0;
        
        // Density from sech² profile
        const real sech_val = 1.0 / std::cosh(x / H);
        p_i.dens = rho_center * sech_val * sech_val;
        
        // ====================================================================
        // For adiabatic gas with γ ≠ 1, the equilibrium pressure is:
        //
        // Hydrostatic equilibrium: dP/dx = -ρ g(x)
        // For sech² profile: g(x) = -2πG ρ₀ H tanh(x/H)
        //
        // Integrating: P(x) = P₀ - 2πG ρ₀² H² [sech²(x/H) - 1]
        //                   = P₀ - 2πG ρ₀² H² sech²(x/H) + 2πG ρ₀² H²
        //
        // At center (x=0): P(0) = P₀
        // At infinity: P(∞) → P₀ + 2πG ρ₀² H²
        //
        // For c_s² = γP/ρ at center to give scale height H:
        //   H = c_s / √(2πGρ₀)  →  c_s² = 2πGρ₀H²
        //   So P₀ = ρ₀ c_s² / γ = 2πGρ₀²H² / γ
        //
        // Full pressure profile:
        //   P(x) = (2πGρ₀²H²/γ) + 2πGρ₀²H²[1 - sech²(x/H)]
        //        = 2πGρ₀²H² [1/γ + 1 - ρ(x)/ρ₀]
        // ====================================================================
        
        const real P0 = 2.0 * M_PI * G * rho_center * rho_center * H * H / gamma;
        const real P_offset = 2.0 * M_PI * G * rho_center * rho_center * H * H;
        p_i.pres = P0 + P_offset * (1.0 - p_i.dens / rho_center);
        
        // Internal energy for adiabatic EOS: P = (γ-1)ρu
        p_i.ene = p_i.pres / ((gamma - 1.0) * p_i.dens);
        
        // Sound speed: c_s = √(γP/ρ)
        p_i.sound = std::sqrt(gamma * p_i.pres / p_i.dens);
        
        // Mass
        p_i.mass = mass_per_particle;
        
        // ID
        p_i.id = i;
        
        // Initial smoothing length estimate
        // For 1D: h ≈ N_ngb * m / (2 * ρ)
        const int n_ngb = m_param->physics.neighbor_number;
        p_i.sml = n_ngb * p_i.mass / (2.0 * p_i.dens);
    }
    
    // ========================================================================
    // Compute initial gravitational acceleration for well-balanced Riemann solver
    // g(x) = 2πG Σ(x) = 2πG ρ_0 H tanh(x/H)
    // ========================================================================
    
    for (int i = 0; i < N; ++i) {
        auto & p_i = particles[i];
        const real x = p_i.pos[0];
        
        // Surface density from center to x
        const real Sigma_x = rho_center * H * std::tanh(x / H);
        
        // Gravitational acceleration (points toward center)
        // g = -2πG Σ * sign(x) = -2πG ρ_0 H tanh(x/H)
        p_i.grav_acc[0] = -2.0 * M_PI * G * Sigma_x;
        
        // Gravitational potential (for diagnostics)
        // φ(x) = 2πG ρ_0 H² [ln(cosh(x/H)) + const]
        p_i.phi = 2.0 * M_PI * G * rho_center * H * H * std::log(std::cosh(x / H));
    }
    
    // ========================================================================
    // Report setup
    // ========================================================================
    
    std::cout << "\nParticle distribution:" << std::endl;
    std::cout << "  x range: [" << particles[0].pos[0] << ", " << particles[N-1].pos[0] << "]" << std::endl;
    std::cout << "  ρ range: [" << particles[N/2].dens << " (center), " 
              << particles[0].dens << " (edge)]" << std::endl;
    
    // Estimate dynamical time
    // The effective sound speed at center is c_eff = sqrt(gamma * P0 / rho0)
    const real P0_center = 2.0 * M_PI * G * rho_center * rho_center * H * H / gamma;
    const real c_eff = std::sqrt(gamma * P0_center / rho_center);
    const real t_dyn = H / c_eff;
    std::cout << "\nTimescales:" << std::endl;
    std::cout << "  c_eff (at center) = " << c_eff << std::endl;
    std::cout << "  t_dyn = H/c_eff = " << t_dyn << std::endl;
    
    // Estimate expected growth rate (theory)
    // Γ = ε * c_s * h / L² where L ~ H
    const real h_avg = particles[N/2].sml;
    const real epsilon = 0.4; // Expected force error without grad-h
    const real Gamma_theory = epsilon * c_eff * h_avg / (H * H);
    std::cout << "  Expected Γ (without grad-h) ≈ " << Gamma_theory << std::endl;
    std::cout << "  Collapse time ≈ " << 1.0 / Gamma_theory << " = " 
              << (1.0 / Gamma_theory) / t_dyn << " t_dyn" << std::endl;
    
    std::cout << "\n=== Isothermal Slab Setup Complete ===" << std::endl;
    
    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());
    
#endif
}

}
