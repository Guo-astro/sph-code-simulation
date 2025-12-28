#pragma once

#include "pre_interaction.hpp"

namespace sph
{
namespace gspmhd
{

/**
 * @brief GSPMHD Pre-interaction module
 *
 * Computes density, smoothing length, and divergence of B for MHD simulations.
 * Based on Iwasaki & Inutsuka (2011).
 */
class PreInteraction : public sph::PreInteraction {
    real m_gamma;
    bool m_use_powell;

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}
