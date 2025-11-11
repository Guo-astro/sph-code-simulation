#pragma once

#include "fluid_force.hpp"
#include <functional>

namespace sph
{
namespace gdisph
{

class FluidForce : public sph::FluidForce {
    bool m_is_2nd_order;
    real m_gamma;

    // (velocity, density, pressure, sound speed)
    std::function<void(const real[], const real[], real & pstar, real & vstar)> m_solver;

    void hll_solver();
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}
