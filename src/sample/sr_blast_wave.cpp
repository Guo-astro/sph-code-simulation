// Special Relativistic Blast Wave Tests
// Based on Kitajima et al. (2025) arXiv:2510.18251v1 Sections 3.1.2 and 3.1.3

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

// Standard blast wave (Section 3.1.2)
void Solver::make_sr_blast_wave()
{
#if DIM != 1
    THROW_ERROR("SR blast wave test requires DIM == 1");
#else

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real gamma = m_param->physics.gamma;
    const real c_speed = m_param->srgsph.c_speed;
    
    // Initial conditions from paper Section 3.1.2
    // Left state:  P=40/3, n=10.0, v=0
    // Right state: P=10^-6, n=1.0, v=0
    // End time: t=0.4
    
    const int N_left = N;
    const int N_right = N;
    const int num = N_left + N_right;
    
    // Domain: x ∈ [-0.5, 0.5], interface at x=0
    const real x_left_start = -0.5;
    const real x_left_end = 0.0;
    const real x_right_start = 0.0;
    const real x_right_end = 0.5;
    
    const real dx_left = (x_left_end - x_left_start) / N_left;
    const real dx_right = (x_right_end - x_right_start) / N_right;
    
    // Left state: high pressure blast
    const real P_left = 40.0 / 3.0;
    const real n_left = 10.0;
    const real v_left = 0.0;
    const real rho_left = n_left;
    
    // Right state: low pressure ambient
    const real P_right = 1.0e-6;
    const real n_right = 1.0;
    const real v_right = 0.0;
    const real rho_right = n_right;
    
    // Same baryon number per particle
    const real nu = 0.5 / n_left / N_left;
    
    std::vector<SPHParticle> p(num);
    
    // Initialize left state (blast region)
    for (int i = 0; i < N_left; ++i) {
        auto& p_i = p[i];
        p_i.id = i;
        p_i.pos[0] = x_left_start + (i + 0.5) * dx_left;
        p_i.mass = nu;
        p_i.nu = nu;
        
        vec_t vel;
        vel[0] = v_left;
        #if DIM >= 2
        vel[1] = 0.0;
        #endif
        #if DIM == 3
        vel[2] = 0.0;
        #endif
        
        sph::srgsph::primitive_to_conserved(
            vel, rho_left, P_left, gamma, c_speed,
            p_i.S, p_i.e, p_i.N
        );
        
        p_i.vel = vel;
        p_i.dens = rho_left;
        p_i.pres = P_left;
        p_i.ene = P_left / ((gamma - 1.0) * rho_left);
        p_i.gamma_lor = 1.0;
    }
    
    // Initialize right state (ambient medium)
    for (int i = 0; i < N_right; ++i) {
        auto& p_i = p[N_left + i];
        p_i.id = N_left + i;
        p_i.pos[0] = x_right_start + (i + 0.5) * dx_right;
        p_i.mass = nu;
        p_i.nu = nu;
        
        vec_t vel;
        vel[0] = v_right;
        #if DIM >= 2
        vel[1] = 0.0;
        #endif
        #if DIM == 3
        vel[2] = 0.0;
        #endif
        
        sph::srgsph::primitive_to_conserved(
            vel, rho_right, P_right, gamma, c_speed,
            p_i.S, p_i.e, p_i.N
        );
        
        p_i.vel = vel;
        p_i.dens = rho_right;
        p_i.pres = P_right;
        p_i.ene = P_right / ((gamma - 1.0) * rho_right);
        p_i.gamma_lor = 1.0;
    }
    
    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
    
    WRITE_LOG << "SR standard blast wave initialized:";
    WRITE_LOG << "  Left:  P=" << P_left << ", n=" << n_left;
    WRITE_LOG << "  Right: P=" << P_right << ", n=" << n_right;
#endif
}

// Strong blast wave (Section 3.1.3)
void Solver::make_sr_strong_blast()
{
#if DIM != 1
    THROW_ERROR("SR strong blast test requires DIM == 1");
#else

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real gamma = m_param->physics.gamma;
    const real c_speed = m_param->srgsph.c_speed;
    
    // Initial conditions from paper Section 3.1.3
    // Left state:  P=1000, n=1.0, v=0
    // Right state: P=0.01, n=1.0, v=0
    // End time: t=0.16
    // Tests C_smooth sensitivity
    
    const int N_left = N;
    const int N_right = N;
    const int num = N_left + N_right;
    
    const real x_left_start = -0.5;
    const real x_left_end = 0.0;
    const real x_right_start = 0.0;
    const real x_right_end = 0.5;
    
    const real dx_left = (x_left_end - x_left_start) / N_left;
    const real dx_right = (x_right_end - x_right_start) / N_right;
    
    // Left state: very high pressure
    const real P_left = 1000.0;
    const real n_left = 1.0;
    const real v_left = 0.0;
    const real rho_left = n_left;
    
    // Right state: very low pressure
    const real P_right = 0.01;
    const real n_right = 1.0;
    const real v_right = 0.0;
    const real rho_right = n_right;
    
    // Same baryon number and density → same mass
    const real nu = 0.5 / n_left / N_left;
    
    std::vector<SPHParticle> p(num);
    
    // Initialize left state
    for (int i = 0; i < N_left; ++i) {
        auto& p_i = p[i];
        p_i.id = i;
        p_i.pos[0] = x_left_start + (i + 0.5) * dx_left;
        p_i.mass = nu;
        p_i.nu = nu;
        
        vec_t vel;
        vel[0] = v_left;
        #if DIM >= 2
        vel[1] = 0.0;
        #endif
        #if DIM == 3
        vel[2] = 0.0;
        #endif
        
        sph::srgsph::primitive_to_conserved(
            vel, rho_left, P_left, gamma, c_speed,
            p_i.S, p_i.e, p_i.N
        );
        
        p_i.vel = vel;
        p_i.dens = rho_left;
        p_i.pres = P_left;
        p_i.ene = P_left / ((gamma - 1.0) * rho_left);
        p_i.gamma_lor = 1.0;
    }
    
    // Initialize right state
    for (int i = 0; i < N_right; ++i) {
        auto& p_i = p[N_left + i];
        p_i.id = N_left + i;
        p_i.pos[0] = x_right_start + (i + 0.5) * dx_right;
        p_i.mass = nu;
        p_i.nu = nu;
        
        vec_t vel;
        vel[0] = v_right;
        #if DIM >= 2
        vel[1] = 0.0;
        #endif
        #if DIM == 3
        vel[2] = 0.0;
        #endif
        
        sph::srgsph::primitive_to_conserved(
            vel, rho_right, P_right, gamma, c_speed,
            p_i.S, p_i.e, p_i.N
        );
        
        p_i.vel = vel;
        p_i.dens = rho_right;
        p_i.pres = P_right;
        p_i.ene = P_right / ((gamma - 1.0) * rho_right);
        p_i.gamma_lor = 1.0;
    }
    
    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
    
    WRITE_LOG << "SR strong blast wave initialized:";
    WRITE_LOG << "  Left:  P=" << P_left << ", n=" << n_left;
    WRITE_LOG << "  Right: P=" << P_right << ", n=" << n_right;
    WRITE_LOG << "  Pressure ratio: " << P_left / P_right;
#endif
}

}
