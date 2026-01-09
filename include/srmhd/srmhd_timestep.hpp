#pragma once

#include "timestep.hpp"
#include "defines.hpp"

namespace sph
{
namespace srmhd
{

/**
 * Special Relativistic MHD Timestep Calculation
 *
 * Uses CFL condition based on relativistic fast magnetosonic speed:
 *
 *   dt = C_CFL * min_i(h_i / c_{f,i})
 *
 * where c_f is the relativistic fast magnetosonic speed:
 *   c_f^2 = (c_s^2 + c_A^2 + sqrt((c_s^2 + c_A^2)^2 - 4*c_s^2*c_A^2))/ 2
 *
 * For typical cases where c_f < c, this is more restrictive than the
 * light-crossing time h/c.
 *
 * Following Kitajima et al. (2025): "a simpler time step criterion,
 * similar to that used in non-relativistic cases, works well"
 */
class TimeStep : public sph::TimeStep {
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;

private:
    real m_c_speed;   // Speed of light
};

}
}
