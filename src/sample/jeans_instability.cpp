/**
 * @file jeans_instability.cpp
 * @brief 1D Jeans instability test for grad-h correction validation
 *
 * Theory (Jeans Instability):
 *   - Uniform gas with small density perturbation
 *   - Perturbations with wavelength > λ_J grow exponentially
 *   - λ_J = c_s * sqrt(π / (G * ρ_0))  (Jeans length)
 *   - Growth rate: ω² = c_s² * k² - 4πGρ₀
 *     - ω² > 0: stable oscillation
 *     - ω² < 0: exponential growth (instability)
 *
 * Test Strategy:
 *   - Create periodic box with ρ(x) = ρ₀(1 + A*sin(kx))
 *   - Choose k such that λ > λ_J (unstable)
 *   - Monitor perturbation amplitude growth
 *   - With grad-h: should match theory
 *   - Without grad-h: may have different growth rate
 */

#include "solver.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "particle.hpp"
#include "exception.hpp"

#include <cmath>
#include <vector>
#include <iostream>

namespace sph {

void Solver::make_jeans_instability()
{
    // ============================================================
    // Parameters from JSON config
    // ============================================================
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real wavelength = boost::any_cast<real>(m_sample_parameters["wavelength"]);
    const real amplitude = boost::any_cast<real>(m_sample_parameters["amplitude"]);
    const real rho_0 = boost::any_cast<real>(m_sample_parameters["rho_0"]);
    const real c_s = boost::any_cast<real>(m_sample_parameters["c_s"]);
    const real gamma = m_param->physics.gamma;
    const real G = m_param->gravity.constant;

    std::cout << "\n=== Jeans Instability Test ===" << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "  N particles     = " << N << std::endl;
    std::cout << "  Wavelength λ    = " << wavelength << std::endl;
    std::cout << "  Amplitude A     = " << amplitude << std::endl;
    std::cout << "  Mean density ρ₀ = " << rho_0 << std::endl;
    std::cout << "  Sound speed c_s = " << c_s << std::endl;
    std::cout << "  G               = " << G << std::endl;
    std::cout << "  γ               = " << gamma << std::endl;

    // Compute Jeans length and growth rate
    const real k = 2.0 * M_PI / wavelength;
    const real lambda_J = c_s * std::sqrt(M_PI / (G * rho_0));
    const real omega_sq = c_s * c_s * k * k - 4.0 * M_PI * G * rho_0;
    
    std::cout << "\nJeans Analysis:" << std::endl;
    std::cout << "  Jeans length λ_J = " << lambda_J << std::endl;
    std::cout << "  λ/λ_J = " << wavelength / lambda_J << std::endl;
    std::cout << "  ω² = c_s²k² - 4πGρ₀ = " << omega_sq << std::endl;
    
    if (omega_sq < 0) {
        real growth_rate = std::sqrt(-omega_sq);
        real e_fold_time = 1.0 / growth_rate;
        std::cout << "  MODE: UNSTABLE (λ > λ_J)" << std::endl;
        std::cout << "  Growth rate Γ = " << growth_rate << std::endl;
        std::cout << "  e-folding time = " << e_fold_time << std::endl;
    } else {
        real oscillation_freq = std::sqrt(omega_sq);
        real oscillation_period = 2.0 * M_PI / oscillation_freq;
        std::cout << "  MODE: STABLE (λ < λ_J)" << std::endl;
        std::cout << "  Oscillation frequency ω = " << oscillation_freq << std::endl;
        std::cout << "  Oscillation period T = " << oscillation_period << std::endl;
    }

    // Domain: one wavelength with periodic boundaries
    const real L = wavelength;
    const real x_min = 0.0;
    const real x_max = L;

    // Compute pressure from sound speed and density
    // c_s² = γP/ρ → P = ρ * c_s² / γ
    const real P_0 = rho_0 * c_s * c_s / gamma;
    const real u_0 = P_0 / ((gamma - 1.0) * rho_0);  // Internal energy

    std::cout << "\nThermodynamics:" << std::endl;
    std::cout << "  Background pressure P₀ = " << P_0 << std::endl;
    std::cout << "  Internal energy u₀ = " << u_0 << std::endl;

    // ============================================================
    // Equal-mass particle placement using CDF inversion
    // ============================================================
    const real total_mass = rho_0 * L;
    const real m_particle = total_mass / N;

    std::vector<real> x_positions(N);
    
    // Newton-Raphson for CDF inversion
    const real a_coeff = amplitude / (2.0 * M_PI);
    
    auto compute_cdf = [&](real x) -> real {
        return x / L - a_coeff * (std::cos(2.0 * M_PI * x / L) - 1.0);
    };
    
    auto compute_pdf = [&](real x) -> real {
        return (1.0 + amplitude * std::sin(2.0 * M_PI * x / L)) / L;
    };

    for (int i = 0; i < N; ++i) {
        real target = (static_cast<real>(i) + 0.5) / N;
        real x = target * L;
        
        for (int iter = 0; iter < 50; ++iter) {
            real F_x = compute_cdf(x);
            real f_x = compute_pdf(x);
            real dx = (target - F_x) / f_x;
            x += dx;
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
        rho_local[i] = rho_0 * (1.0 + amplitude * std::sin(2.0 * M_PI * x_positions[i] / L));
    }

    // ============================================================
    // Create particles
    // ============================================================
    auto& particles = m_sim->get_particles();
    particles.resize(N);
    m_sim->set_particle_num(N);

    const int n_neighbor = m_param->physics.neighbor_number;
    const real h_init = L / N * static_cast<real>(n_neighbor) / 2.0;

    for (int i = 0; i < N; ++i) {
        particles[i].id = i;
        particles[i].mass = m_particle;
        
        particles[i].pos = vec_t(0.0);
        particles[i].pos[0] = x_positions[i];
        
        // Zero initial velocity
        particles[i].vel = vec_t(0.0);
        
        // Uniform pressure, varying density
        particles[i].dens = rho_local[i];
        particles[i].pres = P_0;  // Uniform pressure initially
        particles[i].ene = P_0 / ((gamma - 1.0) * rho_local[i]);
        particles[i].sound = c_s;
        
        particles[i].sml = h_init;
        particles[i].acc = vec_t(0.0);
        particles[i].grav_acc = vec_t(0.0);
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

    // Store useful parameters for post-processing
    m_sample_parameters["L"] = L;
    m_sample_parameters["k"] = k;
    m_sample_parameters["lambda_J"] = lambda_J;
    m_sample_parameters["omega_sq"] = omega_sq;
    if (omega_sq < 0) {
        m_sample_parameters["growth_rate"] = std::sqrt(-omega_sq);
    }

    std::cout << "\n=== Jeans Instability Setup Complete ===" << std::endl;
}

}  // namespace sph
