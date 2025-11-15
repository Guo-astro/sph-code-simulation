#pragma once

#include "pre_interaction.hpp"

namespace sph
{
namespace srgsph
{

/**
 * Special Relativistic GSPH Pre-Interaction Phase
 * Based on Kitajima, Inutsuka, Seno (2025) arXiv:2510.18251v1
 * 
 * Computes:
 * - Particle volume Vp (volume-based approach)
 * - Variable smoothing length h with gather method
 * - Baryon number density N = ν/Vp
 * - Density, pressure, and gradients for MUSCL reconstruction
 */
class PreInteraction : public sph::PreInteraction {
    real m_eta;        // Smoothing length parameter (default: 1.0)
    real m_c_smooth;   // Gradient smoothing parameter (default: 2.0)
    real m_c_speed;    // Speed of light
    real m_gamma;      // EOS gamma
    bool m_first;      // First call flag

    /**
     * Initial smoothing length estimation
     */
    void initial_smoothing(std::shared_ptr<Simulation> sim);

    /**
     * Compute particle volume Vp(x_i) = [Σ_j W(x_i-x_j, h_i)]^(-1)
     * Eq. 33 in paper
     */
    real compute_volume(
        const SPHParticle & p_i,
        const std::vector<SPHParticle> & particles,
        const std::vector<int> & neighbor_list,
        const int n_neighbor,
        const Periodic * periodic,
        const KernelFunction * kernel,
        const real h
    );

    /**
     * Iterative computation of smoothing length
     * Eqs. 35-36: h = η Vp*^(1/d), Vp* uses C_smooth*h
     */
    real compute_smoothing_length(
        const SPHParticle & p_i,
        const std::vector<SPHParticle> & particles,
        const std::vector<int> & neighbor_list,
        const int n_neighbor,
        const Periodic * periodic,
        const KernelFunction * kernel
    );

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    
    /**
     * Main calculation loop
     * Computes volume, smoothing length, baryon number density,
     * and gradients for MUSCL reconstruction
     */
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}
