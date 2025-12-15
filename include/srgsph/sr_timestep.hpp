#pragma once

#include "timestep.hpp"
#include "defines.hpp"

namespace sph
{
namespace srgsph
{

/**
 * Special Relativistic timestep calculation
 * 
 * Following Kitajima, Inutsuka & Seno (2025) SRGSPH paper, we use a simple
 * CFL condition based on sound speed rather than elaborate relativistic
 * characteristic speeds:
 * 
 *   Δt = C_CFL * min_i(h_i / c_{s,i})
 * 
 * The paper notes: "we find that a simpler time step criterion, similar to
 * that used in non-relativistic cases, works well without the more elaborate
 * methods employed in Chow & Monaghan (1997) and Rosswog (2010)."
 * 
 * This is more stable for problems with high tangent velocity where the
 * Pons et al. (2000) characteristic speed formula causes severe timestep
 * collapse as v² → 1.
 */
class TimeStep : public sph::TimeStep {
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;

private:
    real m_c_speed;  // Speed of light (for reference)
};

}
}
