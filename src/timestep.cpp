/**
 * @file timestep.cpp
 * @brief Timestep calculation with CFL, force, and gravity constraints
 *
 * Implements the timestep constraints:
 *   1. Sound CFL: dt_sound = c_sound * h / c_s  (for GSPH: use sound speed, not v_sig)
 *   2. Force CFL: dt_force = 0.3 * sqrt(h / |a_total|)
 *   3. Gravity CFL: dt_grav = 0.3 * sqrt(h / |a_grav|)
 *
 * The final timestep is the minimum of all constraints.
 *
 * For external BH simulations:
 *   - BH acceleration is added to p[i].acc (total acceleration)
 *   - dt_force constraint captures BH effects automatically
 *   - grav_acc contains self-gravity only (if enabled)
 */

#include "parameters.hpp"
#include "timestep.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "openmp.hpp"

#include <algorithm>
#include <cmath>

namespace sph
{

// Force CFL coefficient (fixed at 0.3 for stability)
constexpr real C_FORCE_CFL = 0.3;

void TimeStep::initialize(std::shared_ptr<SPHParameters> param)
{
    c_sound = param->cfl.sound;
    c_force = param->cfl.force;  // kept for compatibility but not used for force CFL
}

void TimeStep::calculation(std::shared_ptr<Simulation> sim)
{
    auto& particles = sim->get_particles();
    const int num = sim->get_particle_num();

    omp_real dt_min(std::numeric_limits<real>::max());
    omp_real h_per_cs_min(std::numeric_limits<real>::max());

#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        // Skip ghost particles for timestep calculation
        if (particles[i].is_ghost) continue;

        // === Sound CFL (GSPH: use h/c_s, not h/v_sig) ===
        // dt_sound = c_sound * h / c_s
        const real cs = particles[i].sound;
        if (cs > 0.0) {
            const real h_per_cs_i = particles[i].sml / cs;
            if (h_per_cs_i < h_per_cs_min.get()) {
                h_per_cs_min.get() = h_per_cs_i;
            }
        }

        // === Force/acceleration CFL ===
        // dt_force = 0.3 * sqrt(h / |a|)
        // where a = total acceleration (fluid + self-gravity + external BH)
        const real acc_abs = std::abs(particles[i].acc);
        if (acc_abs > 0.0) {
            const real dt_force_i = C_FORCE_CFL * std::sqrt(particles[i].sml / acc_abs);
            if (dt_force_i < dt_min.get()) {
                dt_min.get() = dt_force_i;
            }
        }

        // === Gravity CFL (self-gravity only) ===
        // dt_grav = 0.3 * sqrt(h / |a_grav|)
        const real grav_abs = std::abs(particles[i].grav_acc);
        if (grav_abs > 0.0) {
            const real dt_grav_i = C_FORCE_CFL * std::sqrt(particles[i].sml / grav_abs);
            if (dt_grav_i < dt_min.get()) {
                dt_min.get() = dt_grav_i;
            }
        }
    }

    // Sound CFL: dt_sound = c_sound * min(h / c_s)
    const real dt_sound = c_sound * h_per_cs_min.min();

    // Final timestep: minimum of all constraints
    sim->set_dt(std::min(dt_sound, dt_min.min()));
}

}