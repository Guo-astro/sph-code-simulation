#include "srgsph/sr_timestep.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "openmp.hpp"
#include "logger.hpp"
#include <limits>
#include <algorithm>
#include <cmath>

namespace sph
{
namespace srgsph
{

void TimeStep::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::TimeStep::initialize(param);
    m_c_speed = param->srgsph.c_speed;
}

/**
 * Calculate timestep for Special Relativistic GSPH
 *
 * Following Kitajima, Inutsuka & Seno (2025) SRGSPH paper:
 * "we find that a simpler time step criterion, similar to that used in
 * non-relativistic cases, works well without the more elaborate methods
 * employed in Chow & Monaghan (1997) and Rosswog (2010)."
 *
 * The CFL condition used is simply:
 *   Δt = C_CFL * min_i(h_i / c_{s,i})
 *
 * where c_s is the sound speed.
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

        // Skip ghost particles
        if (p_i.is_ghost) continue;

        // Skip if sound speed is zero or negative
        if (p_i.sound <= 0.0) continue;

        // Simple CFL: h / c_s (as in Kitajima et al. 2025)
        const real h_over_cs = p_i.sml / p_i.sound;

        // Update thread-local minimum
        if (min_h_over_cs.get() > h_over_cs) {
            min_h_over_cs.get() = h_over_cs;
        }
    }

    // Apply CFL condition: dt = C_CFL * min(h/c_s)
    const real dt_cfl = c_sound * min_h_over_cs.min();

    // Store minimum h/v_sig for reference (used by base class)
    sim->set_h_per_v_sig(min_h_over_cs.min());

    // Set timestep (CFL limited)
    sim->set_dt(dt_cfl);

    // Debug: report timestep info on first few calls
    static int call_count = 0;
    if (call_count < 5) {
        WRITE_LOG << "[SR TIMESTEP] dt=" << dt_cfl
                  << " min(h/c_s)=" << min_h_over_cs.min()
                  << " CFL=" << c_sound;
        ++call_count;
    }
}

} // namespace srgsph
} // namespace sph
