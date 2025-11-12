// 1D Vacuum Formation Test
// From Hopkins 2015 DISPH paper Section 4.1
// Two gas clouds moving apart with velocities v = ±2.0
// Initial state: rho=1.0, P=0.4, v=-2.0 (left) and v=+2.0 (right)
// Creates vacuum at center as clouds separate

#include "../../include/solver.hpp"
#include "../../include/simulation.hpp"
#include "../../include/particle.hpp"
#include "../../include/exception.hpp"
#include "../../include/parameters.hpp"
#include <cmath>

namespace sph {

void Solver::make_vacuum()
{
#if DIM != 1
    THROW_ERROR("DIM != 1");
#else
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    
    // Place N/2 particles in left region and N/2 in right region
    const int N_left = N / 2;   // 400 particles in left region
    const int N_right = N / 2;  // 400 particles in right region
    const int N_total = N_left + N_right;
    
    // Each region occupies half of the domain [-0.5, 0.5]
    const real left_region_min = -0.5;
    const real left_region_max = 0.0;
    const real right_region_min = 0.0;
    const real right_region_max = 0.5;
    
    const real left_region_size = left_region_max - left_region_min;  // 0.5
    const real right_region_size = right_region_max - right_region_min;  // 0.5
    
    const real dx_left = left_region_size / N_left;
    const real dx_right = right_region_size / N_right;
    
    std::vector<SPHParticle> p(N_total);
    
    // Physical parameters
    const real gamma = m_param->physics.gamma;
    const real density = 1.0;
    const real pressure = 0.4;
    const real velocity_magnitude = 2.0;
    
    // Mass for each particle (total mass = density * total_domain_size)
    const real total_domain_size = 1.0;
    const real mass = density * total_domain_size / N_total;
    
    int particle_id = 0;
    
    // Create left region particles
    real x = left_region_min + dx_left * 0.5;
    for(int i = 0; i < N_left; ++i) {
        auto & p_i = p[particle_id];
        
        p_i.pos[0] = x;
        p_i.dens = density;
        p_i.pres = pressure;
        p_i.mass = mass;
        p_i.ene = pressure / ((gamma - 1.0) * density);
        p_i.vel[0] = -velocity_magnitude;  // Moving left
        p_i.id = particle_id;
        
        x += dx_left;
        particle_id++;
    }
    
    // Create right region particles
    x = right_region_min + dx_right * 0.5;
    for(int i = 0; i < N_right; ++i) {
        auto & p_i = p[particle_id];
        
        p_i.pos[0] = x;
        p_i.dens = density;
        p_i.pres = pressure;
        p_i.mass = mass;
        p_i.ene = pressure / ((gamma - 1.0) * density);
        p_i.vel[0] = velocity_magnitude;   // Moving right
        p_i.id = particle_id;
        
        x += dx_right;
        particle_id++;
    }
    
    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
#endif
}

} // namespace sph
