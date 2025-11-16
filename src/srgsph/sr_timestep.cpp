#include "parameters.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "openmp.hpp"
#include "srgsph/sr_timestep.hpp"

#include <iostream>

namespace sph
{
namespace srgsph
{

void TimeStep::calculation(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    const int num = sim->get_particle_num();

    omp_real dt_min(std::numeric_limits<real>::max());
    
    real min_sound = std::numeric_limits<real>::max();
    real max_sound = 0.0;
    real min_h = std::numeric_limits<real>::max();
    real max_h = 0.0;
    int zero_sound_count = 0;
    
    for(int i = 0; i < num; ++i) {
        min_sound = std::min(min_sound, particles[i].sound);
        max_sound = std::max(max_sound, particles[i].sound);
        min_h = std::min(min_h, particles[i].sml);
        max_h = std::max(max_h, particles[i].sml);
        if(particles[i].sound < 1.0e-10) zero_sound_count++;
    }
    
    static int call_count = 0;
    if(call_count < 5 || zero_sound_count > 0) {
        std::cerr << "Timestep calculation #" << call_count << ":" << std::endl;
        std::cerr << "  Sound speed range: [" << min_sound << ", " << max_sound << "]" << std::endl;
        std::cerr << "  Smoothing length range: [" << min_h << ", " << max_h << "]" << std::endl;
        std::cerr << "  Particles with sound ~ 0: " << zero_sound_count << "/" << num << std::endl;
        call_count++;
    }
    
#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        // Eq. 73 in SR-GSPH paper: dt = C_CFL * min(h_i / c_s_i)
        const real c_s = std::max(particles[i].sound, real(1.0e-10));
        const real dt_sound_i = c_sound * particles[i].sml / c_s;
        if(dt_sound_i < dt_min.get()) {
            dt_min.get() = dt_sound_i;
        }
    }
    
    const real dt_final = dt_min.min();
    if(dt_final < 1e-10) {
        std::cerr << "WARNING: Extremely small timestep: dt = " << dt_final << std::endl;
    }
    
    sim->set_dt(dt_final);
}

}
}
