#pragma once

#include "pre_interaction.hpp"

namespace sph
{
namespace gsph
{

class PreInteraction : public sph::PreInteraction {
    bool m_is_2nd_order;
    bool m_use_gradh;  // Enable grad-h correction factor (default: true)
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}
