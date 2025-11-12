// 2D Pressure Equilibrium Test
// From Hopkins 2015 DISPH paper Section 4.3
// Domain: [0,1] x [0,1]
// High density square region: rho=4.0 in [0.25, 0.75] x [0.25, 0.75]
// Low density background: rho=1.0 outside square
// Uniform pressure: P=2.5 everywhere
// Tests maintenance of pressure equilibrium across density discontinuity

#include "../../include/solver.hpp"
#include "../../include/simulation.hpp"
#include "../../include/particle.hpp"
#include "../../include/exception.hpp"
#include "../../include/parameters.hpp"
#include <cmath>

namespace sph {

void Solver::make_pressure_equilibrium()
{
#if DIM != 2
    THROW_ERROR("DIM != 2");
#else
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real domain_size = 1.0;  // Domain: [0, 1] x [0, 1]
    const real dx = domain_size / N;
    
    const int num = N * N;
    std::vector<SPHParticle> p(num);
    
    // Physical parameters
    const real gamma = m_param->physics.gamma;
    const real pressure = 2.5;  // Uniform pressure everywhere
    
    // Density values
    const real rho_high = 4.0;  // High density in square
    const real rho_low = 1.0;   // Low density outside
    
    // Square region boundaries [0.25, 0.75] x [0.25, 0.75]
    const real square_min = 0.25;
    const real square_max = 0.75;
    
    // Initialize particles on uniform grid
    real x = 0.0 + dx * 0.5;
    real y = 0.0 + dx * 0.5;
    
    for(int i = 0; i < num; ++i) {
        auto & p_i = p[i];
        
        p_i.pos[0] = x;
        p_i.pos[1] = y;
        p_i.vel[0] = 0.0;
        p_i.vel[1] = 0.0;
        p_i.pres = pressure;  // Uniform pressure
        p_i.id = i;
        
        // Check if particle is inside the high-density square
        bool inside_square = (x >= square_min && x <= square_max &&
                            y >= square_min && y <= square_max);
        
        if(inside_square) {
            // High density region
            p_i.dens = rho_high;
            p_i.mass = rho_high * dx * dx;
        } else {
            // Low density background
            p_i.dens = rho_low;
            p_i.mass = rho_low * dx * dx;
        }
        
        p_i.ene = pressure / ((gamma - 1.0) * p_i.dens);
        
        // Move to next grid position
        x += dx;
        if(x >= domain_size) {
            x = 0.0 + dx * 0.5;
            y += dx;
        }
    }
    
    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
#endif
}

} // namespace sph
