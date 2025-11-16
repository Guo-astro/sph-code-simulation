#pragma once

#include "timestep.hpp"

namespace sph
{
namespace srgsph
{

/**
 * Special Relativistic timestep calculation
 * Based on Eq. 73 in Kitajima et al. (2025): dt = C_CFL * min(h_i / c_s_i)
 * 
 * Unlike non-relativistic SPH which uses signal velocity v_sig between particles,
 * SR-GSPH uses a simpler criterion based only on the sound speed at each particle.
 */
class TimeStep : public sph::TimeStep {
public:
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}
