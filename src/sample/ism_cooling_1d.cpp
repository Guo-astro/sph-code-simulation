/**
 * @file ism_cooling_1d.cpp
 * @brief 1D ISM cooling test - CNM thermal relaxation benchmark
 * 
 * Based on Koyama & Inutsuka (2000) Model C10:
 *   - Uniform CNM at n=10 cm^-3, T=107 K
 *   - Small temperature perturbations to trigger relaxation
 *   - Test thermal evolution to equilibrium state
 * 
 * Equal mass particle initialization (like Sod shock tube)
 */

#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"
#include <cmath>
#include <cstdlib>

namespace sph
{

void Solver::make_ism_cooling_1d()
{
#if DIM != 1
    THROW_ERROR("ISM cooling 1D test requires DIM == 1");
#else

    std::cout << "\n=== ISM Cooling 1D Test - CNM Thermal Relaxation ===" << std::endl;
    std::cout << "Based on Koyama & Inutsuka (2000) Model C10" << std::endl;
    std::cout << "Initial: Uniform CNM n=10 cm^-3, T=107 K" << std::endl;
    
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real domain_size = 1.0;  // Domain [0, 1] in code units
    
    // === Physical Parameters (K&I Model C10) ===
    const real n_H_cgs = 10.0;        // Hydrogen number density [cm^-3]
    const real T_init_K = 107.0;      // Initial temperature [K]
    const real mu = 1.27;             // Mean molecular weight (mostly HI)
    
    // In code units: density ~ number density, temperature ~ energy
    // We set density_to_n_H = 1.0 in config, so rho_code = n_H_cgs directly
    const real rho_code = n_H_cgs;    // Code units (number density)
    const real T_code = T_init_K;     // Code units (temperature in K)
    
    // Equal mass particle initialization (like Sod shock tube)
    const real dx = domain_size / N;           // Uniform spacing
    const real mass = rho_code * dx;           // Equal mass per particle
    
    std::cout << "\nPhysical Parameters:" << std::endl;
    std::cout << "  n_H = " << n_H_cgs << " cm^-3" << std::endl;
    std::cout << "  T_init = " << T_init_K << " K" << std::endl;
    std::cout << "  mu = " << mu << std::endl;
    std::cout << "  rho_code = " << rho_code << std::endl;
    std::cout << "  T_code = " << T_code << std::endl;
    
    std::cout << "\nParticle Setup:" << std::endl;
    std::cout << "  N = " << N << " particles" << std::endl;
    std::cout << "  dx = " << dx << " (uniform spacing)" << std::endl;
    std::cout << "  mass = " << mass << " (equal mass)" << std::endl;
    std::cout << "  Total mass = " << N * mass << std::endl;
    
    std::vector<SPHParticle> p(N);
    const real gamma = m_param->physics.gamma;
    
    // Seed random number generator
    srand(42);  // Fixed seed for reproducibility
    
    for(int i = 0; i < N; ++i) {
        auto & p_i = p[i];
        
        // Position: uniform spacing in [0, 1]
        p_i.pos[0] = (i + 0.5) * dx;
        
        // Velocity: zero (hydrostatic equilibrium)
        p_i.vel[0] = 0.0;
        
        // Equal mass initialization
        p_i.mass = mass;
        
        // Uniform density (will be recalculated by initial_smoothing)
        p_i.dens = rho_code;
        
        // Uniform temperature (NO perturbations for thermal relaxation test)
        const real T_particle = T_code;
        
        // Pressure from ideal gas law: P = n * k_B * T / (mu * m_H)
        // In simplified code units: P = n_H * T
        p_i.pres = p_i.dens * T_particle;
        
        // Internal energy: e = P / ((gamma - 1) * rho)
        p_i.ene = p_i.pres / ((gamma - 1.0) * p_i.dens);
        
        p_i.id = i;
        
        // Debug output for first, middle, and last particles
        if(i == 0 || i == N-1 || i == N/2) {
            std::cout << "  Particle " << i << ": x=" << p_i.pos[0] 
                      << ", dens=" << p_i.dens 
                      << ", mass=" << p_i.mass
                      << ", T=" << T_particle << " K"
                      << ", pres=" << p_i.pres
                      << ", ene=" << p_i.ene << std::endl;
        }
    }
    
    std::cout << "\n✓ Particle initialization complete" << std::endl;
    std::cout << "Note: density will be recalculated by initial_smoothing() from neighbors\n" << std::endl;
    
    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
    
    WRITE_LOG << "ISM Cooling 1D Test - CNM Thermal Relaxation";
    WRITE_LOG << "* Based on K&I (2000) Model C10";
    WRITE_LOG << "* Number of particles: " << N;
    WRITE_LOG << "* Initial density: " << n_H_cgs << " cm^-3 (uniform CNM)";
    WRITE_LOG << "* Initial temperature: " << T_init_K << " K (with ±10% perturbations)";
    WRITE_LOG << "* Domain: [0, " << domain_size << "]";
    WRITE_LOG << "* Equal mass initialization: mass = " << mass;

#endif
}

}
