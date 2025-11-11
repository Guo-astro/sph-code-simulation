// Sedov blast wave test (2D)
// Self-similar expansion of a strong shock wave
// See Sedov (1959) "Similarity and Dimensional Methods in Mechanics"

// Initial condition generator for 2D Sedov blast wave
// References: Sedov (1959), Kamm & Timmes (2007)

#include "../../include/solver.hpp"
#include <cmath>

namespace sph {

void Solver::make_sedov()
{
#if DIM != 2
    THROW_ERROR("DIM != 2");
#else
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real domain_size = 1.0;  // Domain: [-0.5, 0.5] x [-0.5, 0.5]
    const real dx = domain_size / N;
    
    const int num = N * N;
    std::vector<SPHParticle> p(num);
    
    // Physical parameters
    const real gamma = m_param->physics.gamma;
    const real background_density = 1.0;
    const real background_pressure = 1.0e-6;  // Very low pressure
    const real background_energy = background_pressure / ((gamma - 1.0) * background_density);
    
    // Blast energy setup
    const real total_energy = 1.0;  // Total energy in the blast
    const real blast_radius = 2.0 * dx;  // Radius of energy deposition
    
    // Initialize uniform grid
    real x = -0.5 + dx * 0.5;
    real y = -0.5 + dx * 0.5;
    const real mass = background_density * domain_size * domain_size / num;
    
    // Count particles in blast region for energy distribution
    int blast_particle_count = 0;
    for(int i = 0; i < num; ++i) {
        real pos_x = x;
        real pos_y = y;
        real r = std::sqrt(pos_x * pos_x + pos_y * pos_y);
        
        if(r < blast_radius) {
            blast_particle_count++;
        }
        
        x += dx;
        if(x > 0.5) {
            x = -0.5 + dx * 0.5;
            y += dx;
        }
    }
    
    // Energy per particle in blast region
    // Ensure at least one particle gets the energy
    if(blast_particle_count == 0) {
        blast_particle_count = 1;
    }
    const real blast_energy_per_particle = total_energy / (blast_particle_count * mass);
    
    // Reset counters
    x = -0.5 + dx * 0.5;
    y = -0.5 + dx * 0.5;
    
    // Set particle properties
    for(int i = 0; i < num; ++i) {
        auto & p_i = p[i];
        
        p_i.pos[0] = x;
        p_i.pos[1] = y;
        p_i.vel[0] = 0.0;
        p_i.vel[1] = 0.0;
        p_i.dens = background_density;
        p_i.mass = mass;
        p_i.id = i;
        
        // Check if particle is in blast region
        real r = std::sqrt(x * x + y * y);
        
        if(r < blast_radius) {
            // High energy in blast region
            p_i.ene = blast_energy_per_particle;
            p_i.pres = (gamma - 1.0) * p_i.dens * p_i.ene;
        } else {
            // Background energy
            p_i.ene = background_energy;
            p_i.pres = background_pressure;
        }
        
        x += dx;
        if(x > 0.5) {
            x = -0.5 + dx * 0.5;
            y += dx;
        }
    }
    
    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
#endif
}

} // namespace sph
