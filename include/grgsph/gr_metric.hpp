#pragma once

#include "defines.hpp"
#include "vector_type.hpp"
#include <cmath>
#include <algorithm>
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

    // Inner product with spatial metric: <u,v> = γ_ij u^i v^j
    real metric_inner_product(const real u[3], const real v[3]) const {
        real result = 0.0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result += gamma_ij[i][j] * u[i] * v[j];
            }
        }
        return result;
    }

    // Norm with spatial metric: |v| = sqrt(γ_ij v^i v^j)
    real metric_norm(const real v[3]) const {
        return std::sqrt(std::max(1e-30, metric_inner_product(v, v)));
    }
};

/**
 * Orthonormal Tetrad for GR-GSPH Riemann Solver
 *
 * Constructs a local orthonormal frame {e^(n), e^(t1), e^(t2)} at each
 * interaction point, using the spatial metric γ_ij for orthonormalization.
 *
 * The tetrad satisfies: γ_ij e^i_(a) e^j_(b) = δ_ab
 *
 * References:
 *   - Liptai & Price (2019) Section 2.4: "The Riemann problem is solved
 *     in a local orthonormal frame aligned with the line of sight"
 */
struct Tetrad {
    real n[3];     // Unit normal (line-of-sight direction) e^i_(n)
    real t1[3];    // First tangent direction e^i_(t1)
    real t2[3];    // Second tangent direction e^i_(t2)

    /**
     * Construct orthonormal tetrad from line-of-sight vector
     *
     * Uses Gram-Schmidt orthonormalization with the spatial metric γ_ij:
     * 1. n^i = los^i / |los|_γ  (normalize line-of-sight)
     * 2. t1^i = (ref - <ref,n>_γ n) / |...|_γ  (Gram-Schmidt)
     * 3. t2^i = ε^ijk n_j t1_k / √γ  (cross product for 3D)
     *
     * @param los Line-of-sight vector (separation r_ij)
     * @param metric Spatial metric for orthonormalization
     */
    void construct(const real los[3], const Metric31& metric);

    /**
     * Project Eulerian velocity onto tetrad frame
     *
     * V^(a) = γ_ij V^i e^j_(a)  (tetrad-frame components)
     *
     * @param V Eulerian velocity V^i in coordinate frame
     * @param V_n Normal component (along line-of-sight)
     * @param V_t1 First tangent component
     * @param V_t2 Second tangent component
     */
    void project_velocity(const real V[3], const Metric31& metric,
                          real& V_n, real& V_t1, real& V_t2) const;

    /**
     * Reconstruct coordinate velocity from tetrad components
     *
     * V^i = V^(n) n^i + V^(t1) t1^i + V^(t2) t2^i
     *
     * @param V_n Normal component
     * @param V_t1 First tangent component
     * @param V_t2 Second tangent component
     * @param V Output: Eulerian velocity in coordinate frame
     */
    void reconstruct_velocity(real V_n, real V_t1, real V_t2, real V[3]) const;
};

/**
 * Metric spatial derivatives for gravitational momentum source terms
 *
 * The momentum source f_i = (√-g / 2N*) T^μν ∂g_μν/∂x^i
 * requires spatial derivatives of the metric components.
 */
struct MetricDerivatives {
    real dg_tt[3];           // ∂g_tt/∂x^i
    real dg_ti[3][3];        // ∂g_ti/∂x^j
    real dgamma_ij[3][3][3]; // ∂γ_ij/∂x^k
};

/**
 * Metric time derivatives for gravitational energy source (work) term
 *
 * The gravitational work Λ = -(√-g / 2N*) T^μν ∂g_μν/∂t (Formulation §5.3)
 * requires time derivatives of the metric components.
 *
 * For stationary spacetimes (Schwarzschild, Kerr, Minkowski), all derivatives are zero.
 * Non-zero for:
 *   - Cosmological backgrounds (FLRW expansion)
 *   - Gravitational wave spacetimes
 *   - Dynamical black hole spacetimes (accretion, mergers)
 */
struct MetricTimeDerivatives {
    real dg_tt_dt;           // ∂g_tt/∂t
    real dg_ti_dt[3];        // ∂g_ti/∂t
    real dgamma_ij_dt[3][3]; // ∂γ_ij/∂t

    // Constructor initializes to zero (stationary spacetime)
    MetricTimeDerivatives() : dg_tt_dt(0.0) {
        for (int i = 0; i < 3; ++i) {
            dg_ti_dt[i] = 0.0;
            for (int j = 0; j < 3; ++j) {
                dgamma_ij_dt[i][j] = 0.0;
            }
        }
    }

    // Check if metric is time-dependent (any non-zero derivative)
    bool is_time_dependent() const {
        constexpr real eps = 1e-30;
        if (std::abs(dg_tt_dt) > eps) return true;
        for (int i = 0; i < 3; ++i) {
            if (std::abs(dg_ti_dt[i]) > eps) return true;
            for (int j = 0; j < 3; ++j) {
                if (std::abs(dgamma_ij_dt[i][j]) > eps) return true;
            }
        }
        return false;
    }
};

/**
 * Abstract base class for spacetime metrics
 *
 * Provides interface for computing:
 *   1. Metric components (α, β^i, γ_ij) at a position
 *   2. Spatial derivatives ∂g_μν/∂x^i for momentum source
 *   3. Time derivatives ∂g_μν/∂t for energy source (gravitational work)
 */
class MetricBase {
public:
    virtual ~MetricBase() = default;

    // Compute metric at position
    virtual void compute(const vec_t& pos, Metric31& metric) const = 0;

    // Compute spatial metric derivatives at position (for momentum source)
    virtual void compute_derivatives(const vec_t& pos, MetricDerivatives& derivs) const = 0;

    // Compute time derivatives of metric at position (for energy source)
    // Default implementation: stationary metric (all time derivatives zero)
    virtual void compute_time_derivatives(const vec_t& pos, real time,
                                          MetricTimeDerivatives& time_derivs) const {
        time_derivs = MetricTimeDerivatives();  // All zeros for stationary metrics
    }

    // Check if metric is stationary (no time dependence)
    // Override in derived classes for time-dependent metrics
    virtual bool is_stationary() const { return true; }

    // Convenience: compute both metric and spatial derivatives
    void compute_all(const vec_t& pos, Metric31& metric, MetricDerivatives& derivs) const {
        compute(pos, metric);
        compute_derivatives(pos, derivs);
    }

    // Full computation including time derivatives (for dynamical spacetimes)
    void compute_full(const vec_t& pos, real time, Metric31& metric,
                      MetricDerivatives& derivs, MetricTimeDerivatives& time_derivs) const {
        compute(pos, metric);
        compute_derivatives(pos, derivs);
        compute_time_derivatives(pos, time, time_derivs);
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
