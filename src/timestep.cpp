#include "parameters.hpp"
#include "timestep.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "openmp.hpp"

#include <algorithm>
#include <cmath>

namespace sph
{

void TimeStep::initialize(std::shared_ptr<SPHParameters> param)
{
    c_sound = param->cfl.sound;
    c_force = param->cfl.force;
}

void TimeStep::calculation(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    const int num = sim->get_particle_num();

    omp_real dt_min(std::numeric_limits<real>::max());
#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        // Total acceleration CFL (includes fluid + gravity from previous step)
        const real acc_abs = std::abs(particles[i].acc);
        if(acc_abs > 0.0) {
            const real dt_force_i = c_force * std::sqrt(particles[i].sml / acc_abs);
            if(dt_force_i < dt_min.get()) {
                dt_min.get() = dt_force_i;
            }
        }
        
        // Separate gravity CFL using grav_acc (ensures gravity dynamics are resolved)
        // This is important for self-gravitating systems
        const real grav_abs = std::abs(particles[i].grav_acc);
        if(grav_abs > 0.0) {
            const real dt_grav_i = c_force * std::sqrt(particles[i].sml / grav_abs);
            if(dt_grav_i < dt_min.get()) {
                dt_min.get() = dt_grav_i;
            }
        }
    }

    const real dt_sound_i = c_sound * sim->get_h_per_v_sig();
    
    sim->set_dt(std::min(dt_sound_i, dt_min.min()));
}

}