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
    m_fixed_dt = param->grgsph.geodesic_mode ? 0.01 : 0.0;  // Fixed dt for geodesic tests
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
    omp_real min_dt_force(std::numeric_limits<real>::max());

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        const auto & p_i = particles[i];

        // Skip ghost particles
        if (p_i.is_ghost) continue;

        // === Sound CFL: h / c_s (as in Kitajima et al. 2025) ===
        if (p_i.sound > 0.0) {
            const real h_over_cs = p_i.sml / p_i.sound;
            if (min_h_over_cs.get() > h_over_cs) {
                min_h_over_cs.get() = h_over_cs;
            }
        }

        // === Force/acceleration CFL: sqrt(h / |a|) ===
        // Important for GR-GSPH with metric source terms (gravity)
        const real acc_abs = std::abs(p_i.acc);
        if (acc_abs > 0.0) {
            const real dt_force_i = c_force * std::sqrt(p_i.sml / acc_abs);
            if (min_dt_force.get() > dt_force_i) {
                min_dt_force.get() = dt_force_i;
            }
        }
    }

    // Apply CFL condition: dt = min(dt_sound, dt_force)
    const real dt_sound = c_sound * min_h_over_cs.min();
    const real dt_force = min_dt_force.min();
    real dt = std::min(dt_sound, dt_force);

    // For geodesic tests, use fixed timestep if specified
    if (m_fixed_dt > 0.0) {
        dt = m_fixed_dt;
    }

    // Store minimum h/v_sig for reference (used by base class)
    sim->set_h_per_v_sig(min_h_over_cs.min());

    // Set timestep (CFL limited)
    sim->set_dt(dt);

    // Debug: report timestep info on first few calls
    static int call_count = 0;
    if (call_count < 5) {
        WRITE_LOG << "[SR TIMESTEP] dt=" << dt
                  << " dt_sound=" << dt_sound
                  << " dt_force=" << dt_force;
        ++call_count;
    }
}

} // namespace srgsph
} // namespace sph
