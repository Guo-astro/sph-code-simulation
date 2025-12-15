/**
 * @file timestep.cpp
 * @brief Timestep calculation with CFL, force, and gravity constraints
 *
 * Implements the timestep constraints described in timestep.hpp:
 *   1. Sound CFL: dt_sound = c_sound * h / v_sig
 *   2. Force CFL: dt_force = c_force * sqrt(h / |a_total|)
 *   3. Gravity CFL: dt_grav = c_force * sqrt(h / |a_grav|)
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

void TimeStep::initialize(std::shared_ptr<SPHParameters> param)
{
    c_sound = param->cfl.sound;
    c_force = param->cfl.force;
}

void TimeStep::calculation(std::shared_ptr<Simulation> sim)
{
    auto& particles = sim->get_particles();
    const int num = sim->get_particle_num();

    omp_real dt_min(std::numeric_limits<real>::max());

#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        // Skip ghost particles for timestep calculation
        if (particles[i].is_ghost) continue;

        // === Force/acceleration CFL ===
        // dt_force = c_force * sqrt(h / |a|)
        // where a = total acceleration (fluid + self-gravity + external BH)
        // This is the PRIMARY constraint for systems with external BH gravity.
        const real acc_abs = std::abs(particles[i].acc);
        if (acc_abs > 0.0) {
            const real dt_force_i = c_force * std::sqrt(particles[i].sml / acc_abs);
            if (dt_force_i < dt_min.get()) {
                dt_min.get() = dt_force_i;
            }
        }

        // === Gravity CFL (self-gravity only) ===
        // dt_grav = c_force * sqrt(h / |a_grav|)
        // where a_grav = self-gravity acceleration from BH-tree
        // Important for self-gravitating clouds to resolve Jeans instability.
        // NOTE: External BH acceleration is in p[i].acc, not grav_acc.
        const real grav_abs = std::abs(particles[i].grav_acc);
        if (grav_abs > 0.0) {
            const real dt_grav_i = c_force * std::sqrt(particles[i].sml / grav_abs);
            if (dt_grav_i < dt_min.get()) {
                dt_min.get() = dt_grav_i;
            }
        }
    }

    // === Sound CFL ===
    // dt_sound = c_sound * min(h / v_sig)
    // where v_sig is signal speed from Riemann solver (stored in simulation)
    const real dt_sound_i = c_sound * sim->get_h_per_v_sig();

    // Final timestep: minimum of all constraints
    sim->set_dt(std::min(dt_sound_i, dt_min.min()));
}

}