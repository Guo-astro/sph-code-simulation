#pragma once

#include "pre_interaction.hpp"
#include "particle.hpp"
#include "defines.hpp"

namespace sph
{
namespace srmhd
{

/**
 * Special Relativistic MHD Pre-Interaction Module
 *
 * Performs per-timestep setup before force calculation:
 * 1. Compute smoothing length h (volume-based iteration)
 * 2. Compute particle volume V_p and lab-frame density N
 * 3. Recover primitive variables from conserved (including B)
 * 4. Compute divergence of B for Powell correction
 * 5. Compute gradients for MUSCL reconstruction (if enabled)
 * 6. Compute relativistic fast magnetosonic speed
 *
 * Based on:
 * - Kitajima et al. (2025) for relativistic density/primitive recovery
 * - Iwasaki & Inutsuka (2011) for MHD-specific quantities
 */
class PreInteraction : public sph::PreInteraction {
    real m_gamma;           // Ratio of specific heats
    real m_c_speed;         // Speed of light
    bool m_is_2nd_order;    // Enable MUSCL reconstruction
    bool m_iteration;       // Enable iterative smoothing length
    real m_eta;             // Smoothing length parameter
    real m_c_smooth;        // Volume sampling factor
    bool m_use_powell;      // Enable Powell divergence cleaning
    bool m_first_calculation;  // First timestep flag

    /**
     * Compute particle volume V_p and grad-h correction factor Omega
     * V_p = [sum_j W(r_ij, h)]^(-1)
     * Omega = 1 / (1 + h * sum_j dW/dh / (d * sum_j W))
     */
    real compute_volume(
        const SPHParticle & p_i,
        const std::vector<SPHParticle> & particles,
        const std::vector<int> & neighbor_list,
        const int n_neighbor,
        const Periodic * periodic,
        const KernelFunction * kernel,
        const real h,
        real & gradh_out
    );

    /**
     * Compute variable smoothing length using Kitajima volume-based approach
     * h = eta * [V_p^*]^(1/d) where V_p^* uses C_smooth*h kernel
     */
    real compute_smoothing_length(
        const SPHParticle & p_i,
        const std::vector<SPHParticle> & particles,
        const std::vector<int> & neighbor_list,
        const int n_neighbor,
        const Periodic * periodic,
        const KernelFunction * kernel,
        const real search_radius
    );

    /**
     * Compute divergence of B using SPH formulation
     * div(B)_i = sum_j m_j * B_parallel^* * F_ij
     */
    real compute_div_B(
        const SPHParticle & p_i,
        const std::vector<SPHParticle> & particles,
        const std::vector<int> & neighbor_list,
        const int n_neighbor,
        const Periodic * periodic,
        const KernelFunction * kernel
    );

    /**
     * Initial smoothing length calculation for first timestep
     */
    void initial_smoothing(std::shared_ptr<Simulation> sim);

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}
