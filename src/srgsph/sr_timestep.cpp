/**
 * SR-GSPH Adaptive Timestep Calculation
 * 
 * Implements CFL condition for relativistic SPH with variable smoothing length.
 * 
 * Key Equation (Kitajima et al. 2025, Eq. 73):
 *   Δt = C_CFL · min_i [h_i / c_s,i]
 * 
 * where:
 * - C_CFL: Courant-Friedrichs-Lewy safety factor (typically 0.3-0.5)
 * - h_i: Variable smoothing length (adapts to local resolution)
 * - c_s,i: Sound speed (signal propagation speed)
 * 
 * Physical Interpretation:
 * - h_i/c_s,i: Time for sound wave to cross kernel support
 * - Adaptive h → Timestep tracks local resolution automatically
 * - Smaller h near shocks/discontinuities → smaller Δt (stability)
 * - Larger h in smooth regions → larger Δt (efficiency)
 * 
 * CFL Condition:
 * - Ensures numerical stability: information propagates < 1 kernel radius per step
 * - Prevents particles from moving faster than numerical domain of dependence
 * - Critical for capturing shocks and contact discontinuities
 * 
 * Variable-h Benefits:
 * - Fixed-h: dt = C_CFL · h_global / max(c_s) → overly conservative in low-res regions
 * - Variable-h: dt adapts to local h_i/c_s,i → optimal timestep everywhere
 * 
 * Parameters:
 * - c_sound: CFL coefficient (set in parameters, typically 0.3-0.5)
 * - Lower c_sound: more stable, slower simulation
 * - Higher c_sound: faster but may violate CFL in extreme flows
 */
#include "defines.hpp"
#include "srgsph/sr_timestep.hpp"
#include "simulation.hpp"
#include "openmp.hpp"

#include <iostream>

namespace sph
{
namespace srgsph
{

/**
 * Compute adaptive timestep using CFL condition
 * 
 * Algorithm:
 * 1. Loop over all particles in parallel
 * 2. Compute local dt_i = h_i / c_s,i for each particle
 * 3. Find global minimum: dt_min = min_i(dt_i)
 * 4. Apply CFL safety factor: Δt = C_CFL · dt_min
 * 5. Set simulation timestep
 * 
 * Edge Cases:
 * - If c_s,i = 0 (vacuum): skip particle (keep dt_min from others)
 * - If all c_s = 0: dt_min = infinity → handled by base class
 * 
 * Performance:
 * - OpenMP reduction on dt_min (thread-safe minimum)
 * - Linear complexity O(N) where N = number of particles
 */
void TimeStep::calculation(std::shared_ptr<Simulation> sim)
{
    // SR-GSPH timestep: dt = C_CFL * min(h_i / c_s_i)
    // Based on Eq. 73 in Kitajima et al. (2025)
    
    auto &particles = sim->get_particles();
    const int num = sim->get_particle_num();
    
    omp_real dt_min(std::numeric_limits<real>::max());
    
#pragma omp parallel for
    for (int i = 0; i < num; ++i)
    {
        const auto &p_i = particles[i];
        
        // Compute timestep based on sound speed
        const real h_i = p_i.sml;
        const real cs_i = p_i.sound;
        
        if (cs_i > 0.0)
        {
            const real dt_i = h_i / cs_i;
            if (dt_min.get() > dt_i)
            {
                dt_min.get() = dt_i;
            }
        }
    }
    
    // Store h/c_s, which will be multiplied by c_sound CFL coefficient in base class
    sim->set_h_per_v_sig(dt_min.min());
    
    // Apply CFL coefficient and set the timestep
    const real dt_sound = c_sound * dt_min.min();
    sim->set_dt(dt_sound);
}

} // namespace srgsph
} // namespace sph
