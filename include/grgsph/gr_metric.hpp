#pragma once

#include "defines.hpp"
#include "vector_type.hpp"
#include <cmath>
#include <memory>

namespace sph {
namespace grgsph {

/**
 * 3+1 Decomposition of Spacetime Metric
 *
 * The line element in 3+1 form:
 *   ds² = -α² dt² + γ_ij (dx^i + β^i dt)(dx^j + β^j dt)
 *
 * where:
 *   α     = lapse function
 *   β^i   = shift vector
 *   γ_ij  = spatial 3-metric
 *
 * References:
 *   - Liptai & Price (2019) "General relativistic SPH"
 *   - Appendix A & B for Schwarzschild/Kerr in Cartesian coordinates
 */
struct Metric31 {
    // Primary metric components
    real alpha;           // Lapse function α
    real beta[3];         // Shift vector β^i
    real gamma_ij[3][3];  // Spatial 3-metric (covariant)
    real gamma_inv[3][3]; // Inverse spatial 3-metric (contravariant)
    real sqrt_gamma;      // √(det γ_ij)
    real sqrt_neg_g;      // √(-g) = α√γ

    // 4-metric components (for source terms)
    real g_tt;            // g_00
    real g_ti[3];         // g_0i (= -α² β_i for standard 3+1)

    // Constructor initializes to Minkowski
    Metric31() {
        alpha = 1.0;
        sqrt_gamma = 1.0;
        sqrt_neg_g = 1.0;
        g_tt = -1.0;
        for (int i = 0; i < 3; ++i) {
            beta[i] = 0.0;
            g_ti[i] = 0.0;
            for (int j = 0; j < 3; ++j) {
                gamma_ij[i][j] = (i == j) ? 1.0 : 0.0;
                gamma_inv[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
    }

    // Compute Lorentz factor from Eulerian velocity V^i
    // Γ = 1/√(1 - γ_ij V^i V^j)
    real lorentz_factor(const real V[3]) const {
        real V_sq = 0.0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                V_sq += gamma_ij[i][j] * V[i] * V[j];
            }
        }
        V_sq = std::min(V_sq, 0.9999);  // Clamp to subluminal
        return 1.0 / std::sqrt(1.0 - V_sq);
    }

    // Convert coordinate velocity v^i to Eulerian velocity V^i
    // V^i = (v^i + β^i) / α
    void coord_to_eulerian(const real v[3], real V[3]) const {
        for (int i = 0; i < 3; ++i) {
            V[i] = (v[i] + beta[i]) / alpha;
        }
    }

    // Convert Eulerian velocity V^i to coordinate velocity v^i
    // v^i = α V^i - β^i
    void eulerian_to_coord(const real V[3], real v[3]) const {
        for (int i = 0; i < 3; ++i) {
            v[i] = alpha * V[i] - beta[i];
        }
    }

    // Lower index: V_i = γ_ij V^j
    void lower_index(const real V_up[3], real V_down[3]) const {
        for (int i = 0; i < 3; ++i) {
            V_down[i] = 0.0;
            for (int j = 0; j < 3; ++j) {
                V_down[i] += gamma_ij[i][j] * V_up[j];
            }
        }
    }

    // Raise index: V^i = γ^ij V_j
    void raise_index(const real V_down[3], real V_up[3]) const {
        for (int i = 0; i < 3; ++i) {
            V_up[i] = 0.0;
            for (int j = 0; j < 3; ++j) {
                V_up[i] += gamma_inv[i][j] * V_down[j];
            }
        }
    }
};

/**
 * Metric derivatives for gravitational source terms
 *
 * The source term f_i = (√-g / 2ρ*) T^μν ∂g_μν/∂x^i
 * requires derivatives of the metric components.
 */
struct MetricDerivatives {
    real dg_tt[3];           // ∂g_tt/∂x^i
    real dg_ti[3][3];        // ∂g_ti/∂x^j
    real dgamma_ij[3][3][3]; // ∂γ_ij/∂x^k
};

/**
 * Abstract base class for spacetime metrics
 */
class MetricBase {
public:
    virtual ~MetricBase() = default;

    // Compute metric at position
    virtual void compute(const vec_t& pos, Metric31& metric) const = 0;

    // Compute metric derivatives at position
    virtual void compute_derivatives(const vec_t& pos, MetricDerivatives& derivs) const = 0;

    // Convenience: compute both metric and derivatives
    void compute_all(const vec_t& pos, Metric31& metric, MetricDerivatives& derivs) const {
        compute(pos, metric);
        compute_derivatives(pos, derivs);
    }
};

/**
 * Minkowski metric (flat spacetime)
 * Used for SR tests and as a reference
 */
class MinkowskiMetric : public MetricBase {
public:
    void compute(const vec_t& pos, Metric31& metric) const override {
        metric = Metric31();  // Default constructor gives Minkowski
    }

    void compute_derivatives(const vec_t& pos, MetricDerivatives& derivs) const override {
        // All derivatives are zero for flat spacetime
        for (int i = 0; i < 3; ++i) {
            derivs.dg_tt[i] = 0.0;
            for (int j = 0; j < 3; ++j) {
                derivs.dg_ti[i][j] = 0.0;
                for (int k = 0; k < 3; ++k) {
                    derivs.dgamma_ij[i][j][k] = 0.0;
                }
            }
        }
    }
};

/**
 * Schwarzschild metric in Cartesian-like coordinates
 *
 * Line element (spherical):
 *   ds² = -(1-2M/r) dt² + dr²/(1-2M/r) + r²(dθ² + sin²θ dφ²)
 *
 * Transformed to Cartesian (x,y,z) following GRSPH Appendix A.
 * This avoids coordinate singularities except at r=0.
 *
 * Key properties:
 *   - α² = 1 - 2M/r (lapse squared)
 *   - β^i = 0 (no shift)
 *   - √(-g) = 1 in Cartesian coordinates
 */
class SchwarzschildMetric : public MetricBase {
    real M;  // Black hole mass

public:
    explicit SchwarzschildMetric(real mass) : M(mass) {}

    void compute(const vec_t& pos, Metric31& metric) const override;
    void compute_derivatives(const vec_t& pos, MetricDerivatives& derivs) const override;

    real get_mass() const { return M; }
    real schwarzschild_radius() const { return 2.0 * M; }
};

/**
 * Kerr metric in Cartesian-like Boyer-Lindquist coordinates
 *
 * Includes black hole spin (angular momentum parameter a).
 * Transformed to Cartesian following GRSPH Appendix B.
 *
 * Key properties:
 *   - a = spin parameter, |a| ≤ M
 *   - Frame dragging: g_tφ ≠ 0 (non-zero shift)
 *   - Ergosphere: α² < 0 outside horizon but inside ergosphere
 */
class KerrMetric : public MetricBase {
    real M;  // Black hole mass
    real a;  // Spin parameter

public:
    KerrMetric(real mass, real spin) : M(mass), a(spin) {}

    void compute(const vec_t& pos, Metric31& metric) const override;
    void compute_derivatives(const vec_t& pos, MetricDerivatives& derivs) const override;

    real get_mass() const { return M; }
    real get_spin() const { return a; }
    real outer_horizon() const { return M + std::sqrt(M*M - a*a); }
    real inner_horizon() const { return M - std::sqrt(M*M - a*a); }
};

/**
 * Factory function to create metric from string name
 */
std::unique_ptr<MetricBase> create_metric(const std::string& name, real M = 1.0, real a = 0.0);

} // namespace grgsph
} // namespace sph
