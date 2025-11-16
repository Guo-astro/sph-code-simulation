// Special Relativistic Sod Shock Tube Test
// Based on Kitajima et al. (2025) arXiv:2510.18251v1 Section 3.1.1

#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"
#include "logger.hpp"
#include "srgsph/sr_primitive_recovery.hpp"
#include <cmath>

namespace sph
{

void Solver::make_sr_sod()
{
#if DIM != 1
    THROW_ERROR("SR Sod test requires DIM == 1");
#else

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real gamma = m_param->physics.gamma;
    const real c_speed = m_param->srgsph.c_speed;
    
    // Different baryon number test (default: false)
    bool different_nu = false;
    if (m_sample_parameters.count("different_nu")) {
        different_nu = boost::any_cast<bool>(m_sample_parameters["different_nu"]);
    }
    
    // Initial conditions from paper Section 3.1.1
    // Left state:  P=1.0, n=1.0, v=0
    // Right state: P=0.1, n=0.125, v=0
    // Particles: 3200 left, 400 right (8:1 ratio)
    
    const int N_left = N * 8;
    const int N_right = N;
    const int num = N_left + N_right;
    
    // Domain: x ∈ [-0.5, 0.5], interface at x=0
    const real x_left_start = -0.5;
    const real x_left_end = 0.0;
    const real x_right_start = 0.0;
    const real x_right_end = 0.5;
    
    const real dx_left = (x_left_end - x_left_start) / N_left;
    const real dx_right = (x_right_end - x_right_start) / N_right;
    
    // Left state primitive variables
    const real P_left = 1.0;
    const real n_left = 1.0;  // baryon number density
    const real v_left = 0.0;
    
    // Right state primitive variables
    const real P_right = 0.1;
    const real n_right = 0.125;
    const real v_right = 0.0;
    
    // Compute rest mass density from baryon number density
    // ρ = n × m_baryon (set m_baryon = 1 for dimensionless units)
    const real rho_left = n_left;
    const real rho_right = n_right;
    
    // Baryon number per particle
    real nu_left, nu_right;
    if (different_nu) {
        // Different baryon numbers: ν_left ≠ ν_right
        // Total baryon in left half = total baryon in right half
        // N_left × ν_left = N_right × ν_right
        nu_left = 1.0;
        nu_right = static_cast<real>(N_left) / static_cast<real>(N_right) * nu_left;
    } else {
        // Same baryon number: ν_left = ν_right
        nu_left = 0.5 / n_left / N_left;  // from total mass in left half
        nu_right = nu_left;
    }
    
    std::vector<SPHParticle> p(num);
    
    // Initialize left state particles
    for (int i = 0; i < N_left; ++i) {
        auto& p_i = p[i];
        p_i.id = i;
        p_i.pos[0] = x_left_start + (i + 0.5) * dx_left;
        
        // Compute SPH mass from baryon number and density
        // mass = ν × n / ρ = ν (since we set ρ = n)
        p_i.mass = nu_left;
        p_i.nu = nu_left;
        
        // Compute conserved variables from primitives
        vec_t vel;
        vel[0] = v_left;
        #if DIM >= 2
        vel[1] = 0.0;
        #endif
        #if DIM == 3
        vel[2] = 0.0;
        #endif
        
        // Create primitive variables struct
        sph::srgsph::PrimitiveVariables prim;
        prim.vel = vel;
        prim.density = rho_left;
        prim.pressure = P_left;
        prim.gamma_lor = 1.0;  // v=0
        
        // Compute enthalpy: H = 1 + u/c² + P/(nc²)
        const real c2 = c_speed * c_speed;
        const real u = P_left / ((gamma - 1.0) * rho_left);  // Specific internal energy
        prim.enthalpy = 1.0 + u / c2 + P_left / (rho_left * c2);
        prim.sound_speed = std::sqrt(gamma * P_left / (rho_left * prim.enthalpy));
        
        // Baryon number density N = γn
        p_i.N = prim.gamma_lor * rho_left;
        
        // Convert to conserved variables
        sph::srgsph::PrimitiveRecovery::primitive_to_conserved(
            prim, p_i.N, c_speed,
            p_i.S, p_i.e
        );
        
        // Store primitive variables (will be recovered in pre-interaction)
        p_i.vel = vel;
        p_i.dens = rho_left;
        p_i.pres = P_left;
        p_i.ene = P_left / ((gamma - 1.0) * rho_left);  // internal energy
        
        // Initial Lorentz factor (v=0 → γ=1)
        p_i.gamma_lor = 1.0;
    }
    
    // Initialize right state particles
    for (int i = 0; i < N_right; ++i) {
        auto& p_i = p[N_left + i];
        p_i.id = N_left + i;
        p_i.pos[0] = x_right_start + (i + 0.5) * dx_right;
        
        p_i.mass = nu_right;
        p_i.nu = nu_right;
        
        vec_t vel;
        vel[0] = v_right;
        #if DIM >= 2
        vel[1] = 0.0;
        #endif
        #if DIM == 3
        vel[2] = 0.0;
        #endif
        
        // Create primitive variables struct
        sph::srgsph::PrimitiveVariables prim;
        prim.vel = vel;
        prim.density = rho_right;
        prim.pressure = P_right;
        prim.gamma_lor = 1.0;  // v=0
        
        // Compute enthalpy
        const real c2 = c_speed * c_speed;
        const real u = P_right / ((gamma - 1.0) * rho_right);
        prim.enthalpy = 1.0 + u / c2 + P_right / (rho_right * c2);
        prim.sound_speed = std::sqrt(gamma * P_right / (rho_right * prim.enthalpy));
        
        // Baryon number density
        p_i.N = prim.gamma_lor * rho_right;
        
        // Convert to conserved variables
        sph::srgsph::PrimitiveRecovery::primitive_to_conserved(
            prim, p_i.N, c_speed,
            p_i.S, p_i.e
        );
        
        p_i.vel = vel;
        p_i.dens = rho_right;
        p_i.pres = P_right;
        p_i.ene = P_right / ((gamma - 1.0) * rho_right);
        
        p_i.gamma_lor = 1.0;
    }
    
    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
    
    WRITE_LOG << "SR Sod shock tube initialized:";
    WRITE_LOG << "  Left:  " << N_left << " particles, P=" << P_left 
             << ", n=" << n_left << ", ν=" << nu_left;
    WRITE_LOG << "  Right: " << N_right << " particles, P=" << P_right 
             << ", n=" << n_right << ", ν=" << nu_right;
    WRITE_LOG << "  Different ν: " << (different_nu ? "YES" : "NO");
#endif
}

}
