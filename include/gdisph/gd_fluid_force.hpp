#pragma once

#include "fluid_force.hpp"
#include <functional>
#include <memory>

namespace sph {
namespace thermal {
    class InoueInutsukaCooling;  // Forward declaration
}
}

namespace sph
{
namespace gdisph
{

class FluidForce : public sph::FluidForce {
    bool m_is_2nd_order;
    real m_gamma;
    
    // Thermal cooling (optional)
    bool m_enable_cooling;
    std::shared_ptr<thermal::InoueInutsukaCooling> m_cooling;
    real m_thermal_relax_time;
    real m_density_to_n_H;

    // (velocity, density, pressure, sound speed)
    std::function<void(const real[], const real[], real & pstar, real & vstar)> m_solver;

    void hll_solver();
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}

