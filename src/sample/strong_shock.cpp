// 1D Strong Shock Test
// From Hopkins 2015 DISPH paper Section 4.2
// Extreme pressure ratio: P_left/P_right = 10,000
// Left state: rho=1.0, P=1000.0, v=0
// Right state: rho=1.0, P=0.1, v=0
// Tests ability to handle strong shocks without oscillations

#include "../../include/solver.hpp"
#include "../../include/simulation.hpp"
#include "../../include/particle.hpp"
#include "../../include/exception.hpp"
#include "../../include/parameters.hpp"
#include <cmath>

namespace sph {

void Solver::make_strong_shock()
{
#if DIM != 1
    THROW_ERROR("DIM != 1");
#else
    // ============================================================================
    // SSOT: Strong Shock Domain Configuration
    // Domain: [-0.5, 0.5] (size = 1.0)
    // This is the Single Source of Truth for the strong shock test domain.
    // All configuration files and visualization scripts must match this domain.
    // ============================================================================
    
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real domain_size = 1.0;  // Domain: [-0.5, 0.5]
    const real dx = domain_size / N;
    
    std::vector<SPHParticle> p(N);
    
    // Physical parameters
    const real gamma = m_param->physics.gamma;
    const real x_discontinuity = 0.0;  // Discontinuity at x=0
    
    // Left state (high pressure)
    const real rho_left = 1.0;
    const real pres_left = 1000.0;
    
    // Right state (low pressure)
    const real rho_right = 1.0;
    const real pres_right = 0.1;
    
    const real mass = domain_size / N;  // Uniform mass (density will vary)
    
    real x = -0.5 + dx * 0.5;
    
    for(int i = 0; i < N; ++i) {
        auto & p_i = p[i];
        
        p_i.pos[0] = x;
        p_i.vel[0] = 0.0;
        p_i.mass = mass;
        p_i.id = i;
        
        // Set properties based on position
        if(x < x_discontinuity) {
            // Left state: high pressure
            p_i.dens = rho_left;
            p_i.pres = pres_left;
        } else {
            // Right state: low pressure
            p_i.dens = rho_right;
            p_i.pres = pres_right;
        }
        
        p_i.ene = p_i.pres / ((gamma - 1.0) * p_i.dens);
        
        x += dx;
    }
    
    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
#endif
}

} // namespace sph
