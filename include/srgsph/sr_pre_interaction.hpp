#pragma once

#include "pre_interaction.hpp"

namespace sph
{
namespace srgsph
{

/**
 * Special Relativistic GSPH Pre-Interaction Phase
 * Based on Kitajima, Inutsuka, Seno (2025) arXiv:2510.18251v1 §2.2
 * 
 * Fixed Smoothing Length Formulation:
 * - Single constant h for all particles and all times
 * - Direct kernel sum for number density: N(x_i) = Σ_j ν_j W(x_i - x_j; h)
 * - Exact momentum and energy conservation
 * 
 * Computes:
 * - Baryon number density N via direct kernel sum (Eq. 10)
 * - Primitive variable recovery (velocity, pressure, sound speed)
 * - Gradients for MUSCL reconstruction
 */
class PreInteraction : public sph::PreInteraction {
    real m_h;             // Fixed smoothing length (constant for all particles)
    real m_gamma;         // EOS gamma
    real m_c_speed;       // Speed of light
    bool m_is_2nd_order;  // Whether to use 2nd order MUSCL reconstruction

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    
    /**
     * Main calculation loop
     * Computes number density via direct kernel sum (Eq. 10),
     * recovers primitives, and computes gradients for MUSCL
     */
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}
