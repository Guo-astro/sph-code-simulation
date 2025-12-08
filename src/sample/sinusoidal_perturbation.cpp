/**
 * @file sinusoidal_perturbation.cpp
 * @brief Periodic box with sinusoidal density perturbation for diffusion studies
 *
 * This initial condition creates a 1D periodic box with a sinusoidal density
 * perturbation. It is designed to study numerical diffusion in GSPH without
 * the grad-h correction.
 *
 * Theory:
 *   - Initial: ρ(x) = ρ_mean * (1 + A*sin(2πx/λ))
 *   - Evolution: ρ(x,t) = ρ_mean * (1 + A*exp(-Γt)*sin(2πx/λ))
 *   - Growth rate: Γ = D_eff * k² where k = 2π/λ
 *   - For GSPH without grad-h: D_eff ≈ ε*c_s*h with ε ≈ 0.35-0.4
 *
 * Advantages over slab tests:
 *   - Periodic boundaries eliminate edge effects
 *   - No gravity needed
 *   - Pure pressure-driven diffusion
 *   - Known analytical solution
 *
 * @author SPH Code Team
 */

#include "solver.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "particle.hpp"
#include "exception.hpp"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>

namespace sph {

void Solver::make_sinusoidal_perturbation()
{
    // ============================================================
    // Parameters from JSON config
    // ============================================================
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real wavelength = boost::any_cast<real>(m_sample_parameters["wavelength"]);
    const real amplitude = boost::any_cast<real>(m_sample_parameters["amplitude"]);
    const real rho_mean = boost::any_cast<real>(m_sample_parameters["rho_mean"]);
    const real pressure = boost::any_cast<real>(m_sample_parameters["pressure"]);
    const real gamma = m_param->physics.gamma;

    std::cout << "\n=== Sinusoidal Perturbation Initial Condition ===" << std::endl;
    std::cout << "  N particles     = " << N << std::endl;
    std::cout << "  Wavelength λ    = " << wavelength << std::endl;
    std::cout << "  Amplitude A     = " << amplitude << std::endl;
    std::cout << "  Mean density    = " << rho_mean << std::endl;
    std::cout << "  Pressure        = " << pressure << std::endl;
    std::cout << "  Gamma           = " << gamma << std::endl;

    // Domain: one wavelength with periodic boundaries
    const real L = wavelength;
    const real x_min = 0.0;
    const real x_max = L;

    // Compute derived quantities
    const real c_s = std::sqrt(gamma * pressure / rho_mean);  // Sound speed
    const real k = 2.0 * M_PI / wavelength;                    // Wavenumber
    const real u = pressure / ((gamma - 1.0) * rho_mean);      // Internal energy per mass

    std::cout << "  Sound speed c_s = " << c_s << std::endl;
    std::cout << "  Wavenumber k    = " << k << std::endl;
    std::cout << "  Internal energy = " << u << std::endl;

    // ============================================================
    // Equal-mass particle placement using CDF inversion
    // ============================================================
    // For sinusoidal density: ρ(x) = ρ_mean * (1 + A*sin(2πx/λ))
    // CDF: F(x) = ∫₀ˣ ρ(x')dx' / ∫₀ᴸ ρ(x')dx'
    //     = [x - (Aλ/2π)(cos(2πx/λ) - 1)] / L
    //     = x/L - (A/2π)(cos(2πx/λ) - 1)
    //
    // For N particles, we want F(x_i) = (i + 0.5)/N
    // We solve this numerically using Newton-Raphson

    const real total_mass = rho_mean * L;  // Total mass in domain
    const real m_particle = total_mass / N;  // Mass per particle

    std::vector<real> x_positions(N);
    
    // Newton-Raphson to solve F(x) = target for each particle
    const real a_coeff = amplitude / (2.0 * M_PI);  // A/(2π)
    
    auto compute_cdf = [&](real x) -> real {
        // F(x) = x/L - (A/2π)(cos(2πx/λ) - 1)
        return x / L - a_coeff * (std::cos(2.0 * M_PI * x / L) - 1.0);
    };
    
    auto compute_pdf = [&](real x) -> real {
        // f(x) = dF/dx = 1/L + (A/L)*sin(2πx/λ)
        return (1.0 + amplitude * std::sin(2.0 * M_PI * x / L)) / L;
    };

    for (int i = 0; i < N; ++i) {
        real target = (static_cast<real>(i) + 0.5) / N;
        
        // Initial guess: uniform placement
        real x = target * L;
        
        // Newton-Raphson iteration
        for (int iter = 0; iter < 50; ++iter) {
            real F_x = compute_cdf(x);
            real f_x = compute_pdf(x);
            real dx = (target - F_x) / f_x;
            x += dx;
            
            // Clamp to domain
            x = std::max(x_min + 1e-10, std::min(x_max - 1e-10, x));
            
            if (std::abs(dx) < 1e-12) break;
        }
        
        x_positions[i] = x;
    }

    // ============================================================
    // Compute local density at particle positions
    // ============================================================
    std::vector<real> rho_local(N);
    for (int i = 0; i < N; ++i) {
        rho_local[i] = rho_mean * (1.0 + amplitude * std::sin(2.0 * M_PI * x_positions[i] / L));
    }

    // ============================================================
    // Create particles
    // ============================================================
    auto& particles = m_sim->get_particles();
    particles.resize(N);
    m_sim->set_particle_num(N);

    // Initial smoothing length estimate (based on neighbor count)
    const int n_neighbor = m_param->physics.neighbor_number;
    const real h_init = L / N * static_cast<real>(n_neighbor) / 2.0;

    for (int i = 0; i < N; ++i) {
        particles[i].id = i;
        particles[i].mass = m_particle;
        
        // Position
        particles[i].pos = vec_t(0.0);
        particles[i].pos[0] = x_positions[i];
        
        // Velocity (initially zero - stationary perturbation)
        particles[i].vel = vec_t(0.0);
        
        // Thermodynamics: uniform pressure, varying density
        particles[i].dens = rho_local[i];
        particles[i].pres = pressure;
        particles[i].ene = pressure / ((gamma - 1.0) * rho_local[i]);  // u = P/[(γ-1)ρ]
        particles[i].sound = std::sqrt(gamma * pressure / rho_local[i]);
        
        // Smoothing length
        particles[i].sml = h_init;
        
        // Initialize other fields
        particles[i].acc = vec_t(0.0);
        particles[i].dene = 0.0;
        particles[i].alpha = m_param->av.alpha;
        particles[i].balsara = 1.0;
        particles[i].is_ghost = false;
    }

    // ============================================================
    // Set periodic boundary conditions
    // ============================================================
    m_param->periodic.is_valid = true;
    m_param->periodic.range_min[0] = x_min;
    m_param->periodic.range_max[0] = x_max;
#if DIM >= 2
    m_param->periodic.range_min[1] = -1.0;
    m_param->periodic.range_max[1] = 1.0;
#endif
#if DIM == 3
    m_param->periodic.range_min[2] = -1.0;
    m_param->periodic.range_max[2] = 1.0;
#endif

    // ============================================================
    // Compute expected diffusion timescale
    // ============================================================
    // For GSPH without grad-h: D_eff ≈ ε*c_s*h with ε ≈ 0.35
    const real epsilon = 0.35;
    const real h_avg = h_init;  // Approximate smoothing length
    const real D_eff = epsilon * c_s * h_avg;
    const real Gamma_theory = D_eff * k * k;
    const real tau_decay = 1.0 / Gamma_theory;
    const real e_folding_time = tau_decay;  // Time for amplitude to decay by factor of e

    std::cout << "\n=== Diffusion Theory Prediction ===" << std::endl;
    std::cout << "  Effective diffusivity D_eff = ε·c_s·h = " << D_eff << std::endl;
    std::cout << "    (assuming ε = " << epsilon << ", h = " << h_avg << ")" << std::endl;
    std::cout << "  Growth rate Γ = D·k² = " << Gamma_theory << std::endl;
    std::cout << "  e-folding time τ = 1/Γ = " << e_folding_time << std::endl;
    std::cout << "  Expected: A(t) = A₀·exp(-Γt)" << std::endl;
    std::cout << "\n  To verify: track amplitude decay and check:" << std::endl;
    std::cout << "    ln(A(t)/A₀) vs t should have slope -Γ" << std::endl;
    std::cout << "================================================\n" << std::endl;

    // Store useful parameters for post-processing
    m_sample_parameters["L"] = L;
    m_sample_parameters["k"] = k;
    m_sample_parameters["c_s"] = c_s;
    m_sample_parameters["D_eff_theory"] = D_eff;
    m_sample_parameters["Gamma_theory"] = Gamma_theory;
    m_sample_parameters["tau_decay"] = tau_decay;
}

}  // namespace sph
