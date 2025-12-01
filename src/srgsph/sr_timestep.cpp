#include "srgsph/sr_timestep.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "openmp.hpp"
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
 * Compute relativistic characteristic speed from Pons et al. (2000) Eq. 6
 * 
 * λ± = [v^x(1-cs²) ± cs√((1-v²)(1 - v²cs² - (v^x)²(1-cs²)))] / (1 - v²cs²)
 * 
 * With high tangent velocity v^t, the total v² = vx² + vt² is large,
 * which significantly increases λ compared to the non-relativistic case.
 * 
 * @param vx Normal velocity component
 * @param vt Tangent velocity magnitude  
 * @param cs Sound speed (in units of c)
 * @param c  Speed of light (default 1.0)
 * @return Maximum absolute characteristic speed
 */
real TimeStep::compute_characteristic_speed(real vx, real vt, real cs, real c)
{
    // Normalize to c=1
    const real vx_norm = vx / c;
    const real vt_norm = vt / c;
    const real cs_norm = cs / c;
    
    const real v2 = vx_norm * vx_norm + vt_norm * vt_norm;
    const real cs2 = cs_norm * cs_norm;
    
    // Check for superluminal - should not happen but be safe
    if (v2 >= 1.0) {
        return c;  // Return speed of light
    }
    
    // Denominator: 1 - v²cs²
    const real denom = 1.0 - v2 * cs2;
    if (std::abs(denom) < 1e-15) {
        return c;  // Degenerate case
    }
    
    // Discriminant under the square root:
    // (1-v²)(1 - v²cs² - vx²(1-cs²))
    // = (1-v²)(1 - v²cs² - vx² + vx²cs²)
    // = (1-v²)(1 - vx² - v²cs² + vx²cs²)
    // = (1-v²)(1 - vx² - cs²(v² - vx²))
    // = (1-v²)(1 - vx² - cs²·vt²)
    const real one_minus_v2 = 1.0 - v2;
    const real inner = 1.0 - v2 * cs2 - vx_norm * vx_norm * (1.0 - cs2);
    const real discriminant = one_minus_v2 * inner;
    
    if (discriminant < 0.0) {
        // This can happen near the light cone - use c
        return c;
    }
    
    const real sqrt_disc = std::sqrt(discriminant);
    const real vx_term = vx_norm * (1.0 - cs2);
    
    // λ+ and λ- (characteristic speeds for right/left going waves)
    const real lambda_plus = (vx_term + cs_norm * sqrt_disc) / denom;
    const real lambda_minus = (vx_term - cs_norm * sqrt_disc) / denom;
    
    // Return maximum absolute value, scaled back by c
    return c * std::max(std::abs(lambda_plus), std::abs(lambda_minus));
}

/**
 * Calculate timestep for Special Relativistic GSPH using proper characteristic speeds
 * 
 * Unlike the simple h/cs criterion, this uses:
 * dt = C_CFL * min_i(h_i / λ_max_i)
 * 
 * where λ_max is the relativistic characteristic speed that accounts for
 * both normal and tangent velocities. This is crucial for problems with
 * high tangent velocity where the simple h/cs criterion gives too large dt.
 */
void TimeStep::calculation(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    const int num = sim->get_particle_num();

    // Thread-safe minimum using omp_real
    omp_real min_h_over_lambda(std::numeric_limits<real>::max());

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        const auto & p_i = particles[i];

        // Skip if sound speed is zero or negative
        if (p_i.sound <= 0.0) continue;
        
        // Get velocity components
        const real vx = std::abs(p_i.vel);  // Normal velocity magnitude
        const real vt = std::abs(p_i.vel_t);  // Tangent velocity
        
        // Compute relativistic characteristic speed
        const real lambda_max = compute_characteristic_speed(vx, vt, p_i.sound, m_c_speed);
        
        // Ensure lambda_max is positive and reasonable
        const real lambda_safe = std::max(lambda_max, p_i.sound);  // At minimum, use sound speed
        
        // Compute h / λ_max for this particle
        const real h_over_lambda = p_i.sml / lambda_safe;

        // Update thread-local minimum
        if (min_h_over_lambda.get() > h_over_lambda) {
            min_h_over_lambda.get() = h_over_lambda;
        }
    }

    // Apply CFL condition
    const real dt_cfl = c_sound * min_h_over_lambda.min();

    // Store minimum h/v_sig for reference (used by base class)
    sim->set_h_per_v_sig(min_h_over_lambda.min());

    // Set timestep (CFL limited)
    sim->set_dt(dt_cfl);
}

} // namespace srgsph
} // namespace sph
