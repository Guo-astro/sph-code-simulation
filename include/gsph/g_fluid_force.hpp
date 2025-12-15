#pragma once

#include <functional>
#include <memory>
#include "fluid_force.hpp"

namespace sph {
namespace thermal {
    class InoueInutsukaCooling;  // Forward declaration
}
}

namespace sph
{
namespace gsph
{

class FluidForce : public sph::FluidForce {
    bool m_is_2nd_order;
    real m_gamma;

    // Thermal cooling (optional) - Inoue & Inutsuka (2008)
    // Analytic fit to Koyama & Inutsuka (2000) cooling
    bool m_enable_cooling;
    std::shared_ptr<thermal::InoueInutsukaCooling> m_cooling;
    real m_density_to_n_H;    // Code density -> n_H [cm^-3]
    real m_u_to_cgs;          // Code energy -> erg/g
    real m_t_to_cgs;          // Code time -> seconds

    // (velocity, density, pressure, sound speed)
    std::function<void(const real[], const real[], real & pstar, real & vstar)> m_solver;

    void hll_solver();
    void iterative_solver();  // van Leer (1997) iterative Riemann solver
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}