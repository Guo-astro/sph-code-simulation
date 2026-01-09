/**
 * SRMHD Timestep Calculation
 *
 * CFL condition based on relativistic fast magnetosonic speed:
 *   dt = C_CFL * min_i(h_i / c_{f,i})
 *
 * where c_f is the fast magnetosonic speed combining sound and Alfven speeds.
 *
 * Following Kitajima et al. (2025): "a simpler time step criterion,
 * similar to that used in non-relativistic cases, works well"
 */

#include "srmhd/srmhd_timestep.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "openmp.hpp"
#include "logger.hpp"
#include <limits>
#include <algorithm>
#include <cmath>

namespace sph
{
namespace srmhd
{

void TimeStep::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::TimeStep::initialize(param);
    m_c_speed = param->srgsph.c_speed;

    std::cout << "[SRMHD TimeStep] Initialized with c=" << m_c_speed << std::endl;
}

void TimeStep::calculation(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    const int num = sim->get_particle_num();

    // Thread-safe minimum using omp_real
    omp_real min_h_over_cf(std::numeric_limits<real>::max());

#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        const auto & p_i = particles[i];

        // Skip ghost particles
        if (p_i.is_ghost) continue;

        // Fast magnetosonic CFL: h / c_f
        // p_i.sound stores the fast magnetosonic speed in SRMHD
        if (p_i.sound > 0.0) {
            const real h_over_cf = p_i.sml / p_i.sound;
            if (min_h_over_cf.get() > h_over_cf) {
                min_h_over_cf.get() = h_over_cf;
            }
        }

        // NOTE: For relativistic MHD, we do NOT use the non-relativistic
        // acceleration-based CFL formula sqrt(h/|a|) because:
        // 1. acc = dS is the rate of change of canonical momentum S = γHv,
        //    NOT the physical acceleration dv/dt
        // 2. |dS| can be arbitrarily large even for small velocity changes
        // 3. The correct CFL is based on signal speeds (fast magnetosonic)
    }

    // Apply CFL condition: dt = C_CFL * min(h / c_f)
    const real dt_fast = c_sound * min_h_over_cf.min();
    real dt = dt_fast;

    // Store minimum h/v_sig for reference
    sim->set_h_per_v_sig(min_h_over_cf.min());

    // Set timestep
    sim->set_dt(dt);

    // Debug: report timestep info and identify culprit particle
    static int call_count = 0;
    if (call_count < 5 || dt < 1e-6) {
        // Find the particle that limits dt
        int culprit = -1;
        real min_ratio = std::numeric_limits<real>::max();
        for (int i = 0; i < num; ++i) {
            const auto & p_i = particles[i];
            if (p_i.is_ghost) continue;
            if (p_i.sound > 0.0) {
                const real ratio = p_i.sml / p_i.sound;
                if (ratio < min_ratio) {
                    min_ratio = ratio;
                    culprit = i;
                }
            }
        }
        if (culprit >= 0 && (call_count < 5 || dt < 1e-6)) {
            const auto& p = particles[culprit];
            std::cout << "[TIMESTEP DEBUG] dt=" << dt << " particle=" << culprit
                      << " x=" << p.pos[0] << " h=" << p.sml << " c_f=" << p.sound
                      << " rho=" << p.dens << " P=" << p.pres
                      << " gamma_lor=" << p.gamma_lor << std::endl;
        }
        ++call_count;
    }
}

} // namespace srmhd
} // namespace sph
