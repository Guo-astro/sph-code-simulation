#pragma once

/**
 * @file timestep.hpp
 * @brief Timestep calculation with CFL and gravity constraints
 *
 * Computes the global timestep as the minimum of multiple constraints:
 *
 * 1. Sound speed CFL (hydro stability):
 *      dt_sound = C_CFL * min_i(h_i / v_sig_i)
 *    where v_sig is the signal speed from the Riemann solver.
 *
 * 2. Force/acceleration limiter (dynamics stability):
 *      dt_force = C_force * min_i(sqrt(h_i / |a_i|))
 *    where a_i is the total acceleration (fluid + gravity + external).
 *
 * 3. Gravity limiter (important when external BH is present):
 *      dt_grav = eta_t * min_i(sqrt(h_i / |a_grav_i|))
 *    where a_grav includes self-gravity and external BH acceleration.
 *
 * For simulations with external black holes:
 *   - The BH acceleration is added to p[i].acc in PointMassBlackHole::calculation()
 *   - dt_force automatically includes BH effects via total acceleration
 *   - Typical safe values: C_CFL ~ 0.2-0.4, C_force ~ 0.1-0.2
 *   - Near BH (r ~ epsilon), use smaller eta_t ~ 0.1 for stability
 *
 * References:
 *   - Monaghan (1992): SPH CFL condition
 *   - Springel (2005): Gravity timestep limiter
 */

#include "defines.hpp"
#include "module.hpp"

namespace sph
{

class TimeStep : public Module {
protected:
    real c_sound;  ///< CFL coefficient for sound speed: dt_s = c_sound * h / v_sig
    real c_force;  ///< CFL coefficient for acceleration: dt_f = c_force * sqrt(h / |a|)

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    virtual void calculation(std::shared_ptr<Simulation> sim) override;
};

}