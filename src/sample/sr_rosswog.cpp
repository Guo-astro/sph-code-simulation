/**
 * @file sr_rosswog.cpp
 * @brief Special Relativistic SPH initial conditions from Rosswog (2010)
 * 
 * Reference: Rosswog, S. (2010) "Conservative, special-relativistic smoothed
 * particle hydrodynamics" arXiv:0907.4890v3
 * 
 * This file implements the following Rosswog tests:
 *   - Test 3: Perturbed shock tube (sinusoidal density perturbation)
 *   - Test 5: Sine wave advection (periodic domain)
 *   - Test 6: Square wave advection (periodic domain)
 *   - Test 7: Simple wave (isentropic flow)
 * 
 * Tests 1, 2, and 4 are implemented in sr_sod.cpp as:
 *   - Test 1: "rosswog_test1" (same as "blast_wave")
 *   - Test 2: "rosswog_test2" (same as "strong_blast")
 *   - Test 4: "rosswog_test4" (wall shock, γ≈50,000)
 */

#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"
#include "logger.hpp"
#include "srgsph/sr_primitive_recovery.hpp"

#include <cmath>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace sph {

void Solver::make_sr_rosswog()
{
#if DIM != 1
    THROW_ERROR("SR Rosswog tests require DIM == 1");
#else

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real gamma = m_param->physics.gamma;
    const real c_speed = m_param->srgsph.c_speed;
    const real c2 = c_speed * c_speed;

    // Get test type from configuration
    std::string test_type = "perturbed_shock_tube";
    if (m_sample_parameters.count("testType")) {
        test_type = boost::any_cast<std::string>(m_sample_parameters["testType"]);
    }

    WRITE_LOG << "=== Rosswog (2010) SR-SPH Test Generator ===";
    WRITE_LOG << "Test type: " << test_type;
    WRITE_LOG << "Total particles requested: " << N;
    WRITE_LOG << "Adiabatic index γ = " << gamma;

    // Domain boundaries
    const real x_min = m_param->periodic.range_min[0];
    const real x_max = m_param->periodic.range_max[0];
    const real domain_length = x_max - x_min;

    // Allocate particles
    auto& particles = m_sim->get_particles();

    if (test_type == "perturbed_shock_tube" || test_type == "test3") {
        // =================================================================
        // Test 3: Perturbed Shock Tube
        // Left:  (N, v, P) = (5, 0, 50)
        // Right: (N, v, P) = (2 + 0.3*sin(50x), 0, 5)
        // Domain: x ∈ [-0.5, 0.5], discontinuity at x = 0
        // =================================================================
        WRITE_LOG << "Generating: PERTURBED SHOCK TUBE (Rosswog Test 3)";
        WRITE_LOG << "Left state:  (n, v, P) = (5, 0, 50)";
        WRITE_LOG << "Right state: (n, v, P) = (2 + 0.3*sin(50x), 0, 5)";

        const real n_left = 5.0;
        const real v_left = 0.0;
        const real P_left = 50.0;

        const real n_right_base = 2.0;
        const real n_right_amp = 0.3;
        const real n_right_freq = 50.0;
        const real v_right = 0.0;
        const real P_right = 5.0;

        const real x_disc = 0.0;  // Discontinuity at x=0

        // Calculate particle distribution based on total mass
        // Left: n_left * (x_disc - x_min) = 5 * 0.5 = 2.5
        // Right: approx 1.0 (integral of perturbed density)
        const real mass_left = n_left * (x_disc - x_min);
        const real mass_right = n_right_base * (x_max - x_disc);
        const real total_mass = mass_left + mass_right;

        const int N_left = static_cast<int>(N * mass_left / total_mass);
        const int N_right = N - N_left;

        WRITE_LOG << "Particle distribution: Left=" << N_left << ", Right=" << N_right;

        // Left side: uniform spacing
        const real dx_left = (x_disc - x_min) / N_left;
        const real particle_mass = n_left * dx_left;

        particles.resize(N);
        int idx = 0;

        // Create left-side particles (uniform density)
        for (int i = 0; i < N_left; ++i) {
            particles[idx].id = idx;
            particles[idx].pos[0] = x_min + (i + 0.5) * dx_left;
            particles[idx].vel[0] = v_left;
            particles[idx].dens = n_left;
            particles[idx].pres = P_left;
            particles[idx].mass = particle_mass;
            particles[idx].ene = P_left / ((gamma - 1.0) * n_left);
            idx++;
        }

        // Right side: non-uniform spacing to match perturbed density
        real x = x_disc;
        while (idx < N && x < x_max) {
            real n_local = n_right_base + n_right_amp * std::sin(n_right_freq * x);
            real dx_local = particle_mass / n_local;
            
            x += dx_local * 0.5;
            if (x >= x_max) break;
            
            particles[idx].id = idx;
            particles[idx].pos[0] = x;
            particles[idx].vel[0] = v_right;
            particles[idx].dens = n_local;
            particles[idx].pres = P_right;
            particles[idx].mass = particle_mass;
            particles[idx].ene = P_right / ((gamma - 1.0) * n_local);
            idx++;
            
            x += dx_local * 0.5;
        }

        particles.resize(idx);
        m_sim->set_particle_num(idx);
        WRITE_LOG << "Created " << idx << " particles for perturbed shock tube";

    } else if (test_type == "sine_advection" || test_type == "test5") {
        // =================================================================
        // Test 5: Sine Wave Advection
        // Periodic domain with smooth sinusoidal density profile
        // v_0 = 0.997 (γ ≈ 12.92)
        // n(x) = 1 + 0.5*sin(2π*x)
        // P = 1 (uniform, cold limit)
        // =================================================================
        WRITE_LOG << "Generating: SINE WAVE ADVECTION (Rosswog Test 5)";
        
        real v_0 = 0.997;
        if (m_sample_parameters.count("advectionVelocity")) {
            v_0 = boost::any_cast<real>(m_sample_parameters["advectionVelocity"]);
        }
        
        const real gamma_lorentz = 1.0 / std::sqrt(1.0 - v_0 * v_0);
        WRITE_LOG << "Advection velocity v_0 = " << v_0 << ", Lorentz factor γ = " << gamma_lorentz;
        
        const real n_base = 1.0;
        const real amplitude = 0.5;
        const real P_0 = 1.0;
        const real total_mass = domain_length;  // For A=0.5, integral is exactly L
        const real particle_mass = total_mass / N;
        
        particles.resize(N);
        
        real x = x_min;
        for (int i = 0; i < N; ++i) {
            real n_local = n_base + amplitude * std::sin(2.0 * M_PI * (x - x_min) / domain_length);
            real dx_local = particle_mass / n_local;
            
            real x_particle = x + dx_local * 0.5;
            if (x_particle > x_max) x_particle = x_max - 1e-10;
            
            n_local = n_base + amplitude * std::sin(2.0 * M_PI * (x_particle - x_min) / domain_length);
            
            particles[i].id = i;
            particles[i].pos[0] = x_particle;
            particles[i].vel[0] = v_0;
            particles[i].dens = n_local;
            particles[i].pres = P_0;
            particles[i].mass = particle_mass;
            particles[i].ene = P_0 / ((gamma - 1.0) * n_local);
            
            x += dx_local;
            if (x >= x_max) break;
        }
        
        m_sim->set_particle_num(N);
        WRITE_LOG << "Created " << N << " particles for sine advection";
        WRITE_LOG << "NOTE: Periodic boundary conditions required!";

    } else if (test_type == "square_advection" || test_type == "test6") {
        // =================================================================
        // Test 6: Square Wave Advection
        // Periodic domain with discontinuous square wave density
        // v_0 = 0.997 (γ ≈ 12.92)
        // n(x) = 4 for x ∈ [0.25, 0.75], n(x) = 1 otherwise
        // P = 1 (uniform)
        // =================================================================
        WRITE_LOG << "Generating: SQUARE WAVE ADVECTION (Rosswog Test 6)";
        
        real v_0 = 0.997;
        if (m_sample_parameters.count("advectionVelocity")) {
            v_0 = boost::any_cast<real>(m_sample_parameters["advectionVelocity"]);
        }
        
        const real gamma_lorentz = 1.0 / std::sqrt(1.0 - v_0 * v_0);
        WRITE_LOG << "Advection velocity v_0 = " << v_0 << ", Lorentz factor γ = " << gamma_lorentz;
        
        const real n_high = 4.0;
        const real n_low = 1.0;
        const real x_low1 = x_min + 0.25 * domain_length;
        const real x_high = x_min + 0.75 * domain_length;
        const real P_0 = 1.0;
        
        // Total mass: 0.25*1 + 0.5*4 + 0.25*1 = 2.5 (for unit domain)
        const real total_mass = 0.25 * n_low * domain_length + 
                                0.5 * n_high * domain_length + 
                                0.25 * n_low * domain_length;
        
        const int N_low1 = static_cast<int>(N * 0.25 * n_low * domain_length / total_mass);
        const int N_high = static_cast<int>(N * 0.5 * n_high * domain_length / total_mass);
        const int N_low2 = N - N_low1 - N_high;
        
        WRITE_LOG << "Particle distribution: Low1=" << N_low1 << ", High=" << N_high << ", Low2=" << N_low2;
        
        particles.resize(N);
        int idx = 0;
        
        // Region 1: [x_min, x_low1), n = 1
        real dx1 = (x_low1 - x_min) / N_low1;
        for (int i = 0; i < N_low1; ++i) {
            particles[idx].id = idx;
            particles[idx].pos[0] = x_min + (i + 0.5) * dx1;
            particles[idx].vel[0] = v_0;
            particles[idx].dens = n_low;
            particles[idx].pres = P_0;
            particles[idx].mass = n_low * dx1;
            particles[idx].ene = P_0 / ((gamma - 1.0) * n_low);
            idx++;
        }
        
        // Region 2: [x_low1, x_high], n = 4
        real dx2 = (x_high - x_low1) / N_high;
        for (int i = 0; i < N_high; ++i) {
            particles[idx].id = idx;
            particles[idx].pos[0] = x_low1 + (i + 0.5) * dx2;
            particles[idx].vel[0] = v_0;
            particles[idx].dens = n_high;
            particles[idx].pres = P_0;
            particles[idx].mass = n_high * dx2;
            particles[idx].ene = P_0 / ((gamma - 1.0) * n_high);
            idx++;
        }
        
        // Region 3: (x_high, x_max], n = 1
        real dx3 = (x_max - x_high) / N_low2;
        for (int i = 0; i < N_low2; ++i) {
            particles[idx].id = idx;
            particles[idx].pos[0] = x_high + (i + 0.5) * dx3;
            particles[idx].vel[0] = v_0;
            particles[idx].dens = n_low;
            particles[idx].pres = P_0;
            particles[idx].mass = n_low * dx3;
            particles[idx].ene = P_0 / ((gamma - 1.0) * n_low);
            idx++;
        }
        
        m_sim->set_particle_num(N);
        WRITE_LOG << "Created " << N << " particles for square advection";
        WRITE_LOG << "NOTE: Periodic boundary conditions required!";

    } else if (test_type == "simple_wave" || test_type == "test7") {
        // =================================================================
        // Test 7: Simple Wave (Isentropic Flow)
        // Smooth, isentropic expansion/compression wave
        // Radiation-dominated EOS recommended (γ = 4/3)
        // v_max = 0.7 (maximum velocity amplitude)
        // =================================================================
        WRITE_LOG << "Generating: SIMPLE WAVE (Rosswog Test 7)";
        WRITE_LOG << "Radiation-dominated EOS recommended (γ = 4/3)";
        
        if (std::abs(gamma - 4.0/3.0) > 0.01) {
            WRITE_LOG << "WARNING: Simple wave test designed for γ = 4/3, got γ = " << gamma;
        }
        
        real v_max = 0.7;
        if (m_sample_parameters.count("vMax")) {
            v_max = boost::any_cast<real>(m_sample_parameters["vMax"]);
        }
        
        WRITE_LOG << "Maximum velocity v_max = " << v_max;
        
        const real n_0 = 1.0;
        const real P_0 = 100.0;
        const real dx = domain_length / N;
        
        particles.resize(N);
        
        for (int i = 0; i < N; ++i) {
            real x = x_min + (i + 0.5) * dx;
            
            // Simple wave velocity profile: v(x) = v_max * sin(π*(x-x_min)/L)
            real v = v_max * std::sin(M_PI * (x - x_min) / domain_length);
            
            // Relativistic density correction for isentropic flow
            real gamma_v_sq = 1.0 / (1.0 - v * v);
            real n_local = n_0 / std::pow(gamma_v_sq, 1.5);
            if (n_local < 0.1 * n_0) n_local = 0.1 * n_0;
            
            // Pressure from isentropic relation
            const real K = P_0 / std::pow(n_0, gamma);
            real P_local = K * std::pow(n_local, gamma);
            
            particles[i].id = i;
            particles[i].pos[0] = x;
            particles[i].vel[0] = v;
            particles[i].dens = n_local;
            particles[i].pres = P_local;
            particles[i].mass = n_0 * dx;
            particles[i].ene = P_local / ((gamma - 1.0) * n_local);
        }
        
        m_sim->set_particle_num(N);
        WRITE_LOG << "Created " << N << " particles for simple wave";

    } else {
        THROW_ERROR("Unknown Rosswog test type: " + test_type + 
                    ". Valid options: perturbed_shock_tube (test3), sine_advection (test5), "
                    "square_advection (test6), simple_wave (test7)");
    }

    // =================================================================
    // Initialize SR-GSPH canonical variables for all particles
    // Following sr_sod.cpp pattern from Kitajima et al. (2025)
    // =================================================================
    auto& p = m_sim->get_particles();
    const int num = m_sim->get_particle_num();

    WRITE_LOG << "Initializing SR-GSPH canonical variables...";
    
    for (int i = 0; i < num; ++i) {
        // Lorentz factor
        real v2 = 0.0;
        for (int d = 0; d < DIM; ++d) {
            v2 += p[i].vel[d] * p[i].vel[d];
        }
        p[i].gamma_lor = 1.0 / std::sqrt(1.0 - v2 / c2);
        
        // Specific internal energy (from pressure and density)
        real u_init = p[i].pres / ((gamma - 1.0) * p[i].dens);
        p[i].ene = u_init;
        
        // Specific enthalpy: H = 1 + u/c² + P/(n*c²)
        real H = 1.0 + u_init / c2 + p[i].pres / (p[i].dens * c2);
        p[i].enthalpy = H;
        
        // Sound speed
        p[i].sound = std::sqrt((gamma - 1.0) * (H - 1.0) / H) * c_speed;
        
        // Baryon number per particle
        // nu = N * V_lab where V_lab is the lab-frame volume element
        // For equal-mass particles: nu ≈ mass / (gamma * n)
        p[i].nu = p[i].mass;  // Will be recalculated during simulation
        
        // Conserved variables (Kitajima et al. 2025 formulation)
        // N = γ * n (comoving number density boosted to lab frame)
        real N_conserved = p[i].gamma_lor * p[i].dens;
        p[i].N = N_conserved;
        p[i].dens = N_conserved;  // SPH uses N as density
        
        // S = γ * H * v (relativistic momentum density)
        vec_t S_conserved;
        for (int d = 0; d < DIM; ++d) {
            S_conserved[d] = p[i].gamma_lor * H * p[i].vel[d];
        }
        p[i].S = S_conserved;
        
        // e = (H * (X * γ² - 1) + 1) / (X * γ)
        // where X = γ_eos / (γ_eos - 1)
        const real X = gamma / (gamma - 1.0);
        real e_conserved = (H * (X * p[i].gamma_lor * p[i].gamma_lor - 1.0) + 1.0) / (X * p[i].gamma_lor);
        p[i].e = e_conserved;
        
        // Initialize derivatives to zero
        p[i].dS = vec_t(0.0);
        p[i].de = 0.0;
        p[i].dS_old = vec_t(0.0);
        p[i].de_old = 0.0;
        
        // Smoothing length (will be recalculated)
        p[i].sml = 2.0 * std::pow(p[i].mass / p[i].N, 1.0 / DIM);
    }
    
    WRITE_LOG << "SR-GSPH initialization complete";

#endif  // DIM != 1
}

}  // namespace sph
