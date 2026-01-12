#pragma once

#include <functional>
#include <memory>
#include "fluid_force.hpp"

namespace sph {
namespace thermal {
    class KoyamaInutsukaCooling;  // Forward declaration
}
}

namespace sph
{
namespace gsph
{

class FluidForce : public sph::FluidForce {
    bool m_is_2nd_order;
    real m_gamma;

    // Volume-based approach (Kitajima et al.)
    bool m_use_volume_based;  // Use V² instead of 1/ρ² in force formula

    // Thermal cooling (optional) - Koyama & Inutsuka (2000)
    // Full tabulated cooling with column density dependence
    bool m_enable_cooling;
    std::shared_ptr<thermal::KoyamaInutsukaCooling> m_cooling;
    real m_density_to_n_H;    // Code density -> n_H [cm^-3]
    real m_u_to_cgs;          // Code energy -> erg/g
    real m_t_to_cgs;          // Code time -> seconds
    real m_N_H_column;        // Column density [cm^-2]

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