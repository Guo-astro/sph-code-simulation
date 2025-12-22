#pragma once

#include "pre_interaction.hpp"
#include <vector>

namespace sph
{
namespace srgsph
{

/**
 * Special Relativistic GSPH Pre-Interaction Phase
 * Based on Kitajima, Inutsuka, Seno (2025) arXiv:2510.18251v1
 * 
 * Volume-Based Variable Smoothing Length Formulation (§2.3):
 * - Smoothing length h(x) varies in space: h(x) = η * [V_p*(x)]^(1/d)
 * - Particle volume: V_p(x) = [Σ_j W(x-x_j, h(x))]^(-1)  (Eq. 220)
 * - Volume for h calculation: V_p*(x) = [Σ_j W(x-x_j, C_smooth*h(x))]^(-1)  (Eq. 233)
 * - Number density: N(x) = ν(x) / V_p(x)  (Eq. 239)
 * - Parameters: η=1.0, C_smooth=2.0
 * 
 * Computes:
 * - Variable smoothing length h via iterative volume-based approach
 * - Baryon number density N from particle volume (not kernel sum)
 * - Primitive variable recovery (velocity, pressure, sound speed)
 * - Gradients for MUSCL reconstruction
 */
class PreInteraction : public sph::PreInteraction {
protected:
    real m_eta;           // Smoothing length parameter η (default: 1.0)
    real m_c_smooth;      // Kernel expansion for h iteration C_smooth (default: 2.0)
    real m_gamma;         // EOS gamma
    real m_c_speed;       // Speed of light
    bool m_is_2nd_order;  // Whether to use 2nd order MUSCL reconstruction
    bool m_iteration;     // Whether to iterate smoothing length (from iterativeSmoothingLength)
    bool m_first;         // First timestep flag
    bool m_first_calculation; // First calculation flag (for initializing conserved from primitives)

    /**
     * Compute particle volume V_p(x_i) = [Σ_j W(x_i - x_j, h)]^(-1)
     * and grad-h correction factor Ω_i = 1 / (1 + h * Σ dW/dh / (D * Σ W))
     * Used for number density: N_i = ν_i / V_p(x_i)
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
     * Iteratively solves: h = η * [V_p*]^(1/d) where V_p* = [Σ_j W(x-x_j, C_smooth*h)]^(-1)
     * (Eqs. 231-233 in Kitajima et al.)
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
     * Initial smoothing length calculation for first timestep
     * Uses h = η * (ν/N)^(1/d) as initial guess
     */
    void initial_smoothing(std::shared_ptr<Simulation> sim);

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    
    /**
     * Main calculation loop
     * - First step: iterate Eqs. (231-233) until convergence for each particle
     * - Later steps: use previous h as initial guess
     * - Compute N = ν/V_p (volume-based, not kernel sum)
     * - Recover primitives and compute gradients for MUSCL
     */
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}
