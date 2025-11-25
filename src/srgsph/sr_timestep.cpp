#include "srgsph/sr_timestep.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "openmp.hpp"
#include <limits>
#include <algorithm>

namespace sph
{
namespace srgsph
{

void TimeStep::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::TimeStep::initialize(param);
}

/**
 * Calculate timestep for Special Relativistic GSPH
 * Based on Eq. 73 in Kitajima et al. (2025):
 * dt = C_CFL * min_i(h_i / c_{s,i})
 *
 * Unlike non-relativistic SPH which uses signal velocity between particles,
 * SR-GSPH uses a simpler criterion based only on sound speed at each particle.
 */
void TimeStep::calculation(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    const int num = sim->get_particle_num();

    // Thread-safe minimum using omp_real
    omp_real min_h_over_cs(std::numeric_limits<real>::max());

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        const auto & p_i = particles[i];

        // Compute h / c_s for this particle
        if (p_i.sound > 0.0) {
            const real h_over_cs = p_i.sml / p_i.sound;

            // Update thread-local minimum
            if (min_h_over_cs.get() > h_over_cs) {
                min_h_over_cs.get() = h_over_cs;
            }
        }
    }

    // Apply CFL condition
    const real dt_sound = c_sound * min_h_over_cs.min();

    // Store minimum h/v_sig for reference (used by base class)
    sim->set_h_per_v_sig(min_h_over_cs.min());

    // Set timestep (CFL limited)
    sim->set_dt(dt_sound);
}

} // namespace srgsph
} // namespace sph
