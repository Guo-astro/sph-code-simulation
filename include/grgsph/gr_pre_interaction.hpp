#pragma once

#include "srgsph/sr_pre_interaction.hpp"
#include "grgsph/gr_metric.hpp"
#include <memory>

namespace sph {
namespace grgsph {

/**
 * GR-GSPH Pre-Interaction
 *
 * Extends SR-GSPH pre-interaction with metric-aware computations:
 * - Lorentz factor uses proper metric: Γ = 1/√(1 - γ_ij V^i V^j)
 * - Coordinate to Eulerian velocity transformation: V^i = (v^i + β^i)/α
 * - Conserved variables computed with metric-aware Lorentz factor
 *
 * This ensures consistency between the pre-interaction (primitive <-> conserved)
 * and the force calculation (which uses metric for source terms).
 */
class GRPreInteraction : public srgsph::PreInteraction {
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;

    // Get/set metric
    void set_metric(std::unique_ptr<MetricBase> metric) { m_metric = std::move(metric); }
    void set_metric_ptr(const MetricBase* metric) { m_metric_ptr = metric; }
    const MetricBase* get_metric() const {
        return m_metric ? m_metric.get() : m_metric_ptr;
    }

private:
    std::unique_ptr<MetricBase> m_metric;  // Owned metric (if created here)
    const MetricBase* m_metric_ptr = nullptr;  // Non-owning pointer (if shared with force)
};

} // namespace grgsph
} // namespace sph
