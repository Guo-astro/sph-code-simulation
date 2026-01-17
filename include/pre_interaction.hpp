#pragma once

#include <vector>

#include "module.hpp"
#include "particle.hpp"

namespace sph
{
class Periodic;
class KernelFunction;

class PreInteraction : public Module {
protected:
    bool m_use_balsara_switch;
    bool m_use_time_dependent_av;
    bool m_use_gradh;  // Enable grad-h correction (Ω factor)
    real m_alpha_max;
    real m_alpha_min;
    real m_epsilon; // tau = h / (epsilon * c)
    real m_gamma;
    int  m_neighbor_number;
    real m_kernel_ratio;
    real m_c_smooth;  // C_smooth: smoothing length expansion factor for h-adaptation
    bool m_iteration;
    bool m_first;
    bool m_preserve_initial_density;  // Skip density recalculation (for shock tubes)

    // Jeans length resolution check (Truelove et al. 1997 criterion)
    bool m_gravity_enabled;
    real m_G;  // Gravitational constant
    bool m_jeans_warning_issued;  // Avoid spamming warnings

    virtual real newton_raphson(
        const SPHParticle & p_i,
        const std::vector<SPHParticle> & particles,
        const std::vector<int> & neighbor_list,
        const int n_neighbor,
        const Periodic * periodic,
        const KernelFunction * kernel
    );
    void initial_smoothing(std::shared_ptr<Simulation> sim);

public:
    virtual void initialize(std::shared_ptr<SPHParameters> param) override;
    virtual void calculation(std::shared_ptr<Simulation> sim) override;
};
}
