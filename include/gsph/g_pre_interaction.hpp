#pragma once

#include "pre_interaction.hpp"

namespace sph
{
namespace gsph
{

class PreInteraction : public sph::PreInteraction {
    bool m_is_2nd_order;
    bool m_use_gradh;  // Enable grad-h correction factor (default: true)

    // Volume-based approach parameters (Kitajima et al.)
    bool m_use_volume_based;  // Use volume-based density instead of mass-based
    real m_eta;       // Smoothing length parameter: h = η * V_p^(1/d)
    real m_c_smooth;  // Smoothing expansion factor: V_p* uses h* = c_smooth * h
    bool m_first_calculation;  // First timestep flag

    // Helper functions for volume-based approach
    real compute_volume(
        const SPHParticle & p_i,
        const std::vector<SPHParticle> & particles,
        const std::vector<int> & neighbor_list,
        int n_neighbor,
        const class Periodic * periodic,
        const class KernelFunction * kernel,
        real h,
        real & gradh_out);

    real compute_smoothing_length_volume(
        const SPHParticle & p_i,
        const std::vector<SPHParticle> & particles,
        const std::vector<int> & neighbor_list,
        int n_neighbor,
        const class Periodic * periodic,
        const class KernelFunction * kernel);

protected:
    // Override initial_smoothing to support volume-based density
    void initial_smoothing(std::shared_ptr<Simulation> sim);

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}
