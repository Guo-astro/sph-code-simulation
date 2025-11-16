// Ultra-Relativistic Shock Tests
// Based on Kitajima et al. (2025) arXiv:2510.18251v1 Section 3.2

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

void Solver::make_sr_ultra_relat()
{
#if DIM != 1
    THROW_ERROR("Ultra-relativistic shock test requires DIM == 1");
#else

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real gamma = m_param->physics.gamma;
    const real c_speed = m_param->srgsph.c_speed;
    
    // Velocity of left state (default: 0.9c)
    real v_left_input = 0.9;
    if (m_sample_parameters.count("v_left")) {
        v_left_input = boost::any_cast<real>(m_sample_parameters["v_left"]);
    }
    
    // Initial conditions from paper Section 3.2
    // Left state:  v_x = v_left (0.9 to 0.999999999), P=10, n=1.0
    // Right state: v_x = 0, P=10, n=1.0
    // Tests extreme Lorentz factors
    
    const int N_left = N;
    const int N_right = N;
    const int num = N_left + N_right;
    
    const real x_left_start = -0.5;
    const real x_left_end = 0.0;
    const real x_right_start = 0.0;
    const real x_right_end = 0.5;
    
    const real dx_left = (x_left_end - x_left_start) / N_left;
    const real dx_right = (x_right_end - x_right_start) / N_right;
    
    // Left state: moving at relativistic speed
    const real P_left = 10.0;
    const real n_left = 1.0;
    const real v_left = v_left_input * c_speed;
    const real rho_left = n_left;
    
    // Right state: at rest
    const real P_right = 10.0;
    const real n_right = 1.0;
    const real v_right = 0.0;
    const real rho_right = n_right;
    
    // Compute Lorentz factor for left state
    const real v2_left = v_left * v_left;
    const real gamma_lor_left = 1.0 / std::sqrt(1.0 - v2_left / (c_speed * c_speed));
    
    // Same baryon number
    const real nu = 0.5 / n_left / N_left;
    
    std::vector<SPHParticle> p(num);
    
    // Initialize left state (relativistic flow)
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
        p_i.gamma_lor = gamma_lor_left;
    }
    
    // Initialize right state (at rest)
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
    
    WRITE_LOG << "Ultra-relativistic shock initialized:";
    WRITE_LOG << "  Left:  v=" << v_left / c_speed << "c, γ=" << gamma_lor_left 
             << ", P=" << P_left << ", n=" << n_left;
    WRITE_LOG << "  Right: v=0, γ=1, P=" << P_right << ", n=" << n_right;
#endif
}

}
