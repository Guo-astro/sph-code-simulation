// 2D Shock Tube (Sod shock tube extended to 2D)
// Classic Riemann problem in 2D: left and right states separated at x=0.5
// Domain: [0,1] x [0,0.2], 200x40 particles
// Left state: rho=1.0, P=1.0, v=0
// Right state: rho=0.125, P=0.1, v=0

#include "../../include/solver.hpp"
#include "../../include/simulation.hpp"
#include "../../include/particle.hpp"
#include "../../include/exception.hpp"
#include "../../include/parameters.hpp"
#include <cmath>

namespace sph {

void Solver::make_shock_tube_2d()
{
#if DIM != 2
    THROW_ERROR("DIM != 2");
#else
    const int Nx = boost::any_cast<int>(m_sample_parameters["Nx"]);
    const int Ny = boost::any_cast<int>(m_sample_parameters["Ny"]);
    const real Lx = 1.0;  // Domain: [0, 1] in x
    const real Ly = 0.2;  // Domain: [0, 0.2] in y
    const real dx = Lx / Nx;
    const real dy = Ly / Ny;
    
    const int num = Nx * Ny;
    std::vector<SPHParticle> p(num);
    
    // Physical parameters
    const real gamma = m_param->physics.gamma;
    const real x_discontinuity = 0.5;  // Discontinuity at x=0.5
    
    // Left state (x < 0.5)
    const real rho_left = 1.0;
    const real pres_left = 1.0;
    
    // Right state (x > 0.5)
    const real rho_right = 0.125;
    const real pres_right = 0.1;
    
    // Initialize particles on uniform grid
    real x = 0.0 + dx * 0.5;
    real y = 0.0 + dy * 0.5;
    
    for(int i = 0; i < num; ++i) {
        auto & p_i = p[i];
        
        p_i.pos[0] = x;
        p_i.pos[1] = y;
        p_i.vel[0] = 0.0;
        p_i.vel[1] = 0.0;
        p_i.id = i;
        
        // Set properties based on position
        if(x < x_discontinuity) {
            // Left state
            p_i.dens = rho_left;
            p_i.pres = pres_left;
            p_i.mass = rho_left * dx * dy;
        } else {
            // Right state
            p_i.dens = rho_right;
            p_i.pres = pres_right;
            p_i.mass = rho_right * dx * dy;
        }
        
        p_i.ene = p_i.pres / ((gamma - 1.0) * p_i.dens);
        
        // Move to next grid position
        x += dx;
        if(x >= Lx) {
            x = 0.0 + dx * 0.5;
            y += dy;
        }
    }
    
    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
#endif
}

} // namespace sph
