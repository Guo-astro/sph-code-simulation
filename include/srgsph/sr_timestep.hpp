#pragma once

#include "timestep.hpp"
#include "defines.hpp"

namespace sph
{
namespace srgsph
{

/**
 * Special Relativistic timestep calculation with proper characteristic speeds
 * 
 * Uses the relativistic characteristic velocity from Pons et al. (2000) Eq. 6:
 * λ± = [v^x(1-cs²) ± cs√((1-v²)(1 - v²cs² - (v^x)²(1-cs²)))] / (1 - v²cs²)
 * 
 * This properly accounts for tangent velocity v^t which reduces the available
 * "velocity budget" for the normal direction and affects wave propagation speeds.
 * 
 * The CFL condition becomes: dt = C_CFL * min_i(h_i / λ_max_i)
 * where λ_max is the maximum characteristic speed at each particle.
 */
class TimeStep : public sph::TimeStep {
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
    
    /**
     * Compute relativistic characteristic speed (maximum of λ+ and λ-)
     * From Pons et al. (2000) Eq. 6
     * 
     * @param vx Normal velocity component
     * @param vt Tangent velocity magnitude (can be 0)
     * @param cs Sound speed
     * @param c  Speed of light
     * @return Maximum characteristic speed |λ|
     */
    static real compute_characteristic_speed(real vx, real vt, real cs, real c = 1.0);

private:
    real m_c_speed;  // Speed of light
};

}
}
