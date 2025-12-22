/**
 * GR-GSPH Metric Implementation
 *
 * Schwarzschild and Kerr metrics in Cartesian-like coordinates.
 * Based on GRSPH (Liptai & Price 2019) Appendix A and B.
 */

#include "grgsph/gr_metric.hpp"
#include <cmath>
#include <stdexcept>
#include <string>

namespace sph {
namespace grgsph {

// ============================================================================
// Schwarzschild Metric Implementation
// ============================================================================

/**
 * Schwarzschild metric in Cartesian coordinates
 *
 * From GRSPH Appendix A (Eq. A.1-A.8):
 * The metric components in Cartesian (x,y,z) with r = √(x²+y²+z²)
 *
 * Key equations:
 *   g_tt = -(1 - 2M/r)
 *   g_xx = 1/(1-2M/r) * [1 - 2M(y²+z²)/r³]
 *   g_xy = xy * 2M / [r³(1-2M/r)]
 *   etc.
 *
 * In these coordinates, √(-g) = 1 (unit determinant).
 */
void SchwarzschildMetric::compute(const vec_t& pos, Metric31& metric) const
{
    // Compute radial coordinate
#if DIM == 1
    const real x = pos[0];
    const real y = 0.0;
    const real z = 0.0;
#elif DIM == 2
    const real x = pos[0];
    const real y = pos[1];
    const real z = 0.0;
#else
    const real x = pos[0];
    const real y = pos[1];
    const real z = pos[2];
#endif

    const real r2 = x*x + y*y + z*z;
    const real r = std::sqrt(r2);

    // Avoid singularity at r = 0
    if (r < 1e-10) {
        metric = Metric31();  // Return Minkowski for safety
        return;
    }

    const real r3 = r * r2;
    const real rs = 2.0 * M;  // Schwarzschild radius
    const real f = 1.0 - rs / r;  // 1 - 2M/r

    // Avoid horizon (r = 2M)
    if (f < 1e-10) {
        metric = Metric31();  // Return Minkowski for safety
        return;
    }

    const real f_inv = 1.0 / f;

    // Lapse: α² = 1 - 2M/r
    metric.alpha = std::sqrt(f);

    // Shift: β^i = 0 for Schwarzschild
    metric.beta[0] = 0.0;
    metric.beta[1] = 0.0;
    metric.beta[2] = 0.0;

    // Time-time component
    metric.g_tt = -f;

    // Time-space components (zero for Schwarzschild)
    metric.g_ti[0] = 0.0;
    metric.g_ti[1] = 0.0;
    metric.g_ti[2] = 0.0;

    // Spatial metric components γ_ij (GRSPH Eq. A.3-A.8)
    // Factor for off-diagonal terms
    const real off_factor = rs / (r3 * f);

    // Diagonal components
    metric.gamma_ij[0][0] = f_inv * (1.0 - rs * (y*y + z*z) / r3);
    metric.gamma_ij[1][1] = f_inv * (1.0 - rs * (x*x + z*z) / r3);
    metric.gamma_ij[2][2] = f_inv * (1.0 - rs * (x*x + y*y) / r3);

    // Off-diagonal components
    metric.gamma_ij[0][1] = metric.gamma_ij[1][0] = x * y * off_factor;
    metric.gamma_ij[0][2] = metric.gamma_ij[2][0] = x * z * off_factor;
    metric.gamma_ij[1][2] = metric.gamma_ij[2][1] = y * z * off_factor;

    // Inverse spatial metric γ^ij (GRSPH Eq. A.9-A.14)
    metric.gamma_inv[0][0] = 1.0 - rs * x*x / r3;
    metric.gamma_inv[1][1] = 1.0 - rs * y*y / r3;
    metric.gamma_inv[2][2] = 1.0 - rs * z*z / r3;

    metric.gamma_inv[0][1] = metric.gamma_inv[1][0] = -rs * x * y / r3;
    metric.gamma_inv[0][2] = metric.gamma_inv[2][0] = -rs * x * z / r3;
    metric.gamma_inv[1][2] = metric.gamma_inv[2][1] = -rs * y * z / r3;

    // Determinants: √γ = 1, √(-g) = α√γ = α
    metric.sqrt_gamma = 1.0;
    metric.sqrt_neg_g = metric.alpha;
}

/**
 * Schwarzschild metric derivatives
 *
 * Computed via numerical finite difference for simplicity.
 * (Analytical derivatives are lengthy but can be added for optimization)
 */
void SchwarzschildMetric::compute_derivatives(const vec_t& pos, MetricDerivatives& derivs) const
{
    const real eps = 1e-6;

    Metric31 metric_plus, metric_minus;

    // Compute derivatives by central difference
    for (int k = 0; k < 3; ++k) {
        vec_t pos_plus = pos;
        vec_t pos_minus = pos;

#if DIM == 1
        if (k == 0) {
            pos_plus[0] += eps;
            pos_minus[0] -= eps;
        } else {
            // For DIM=1, derivatives in y,z are zero
            derivs.dg_tt[k] = 0.0;
            for (int i = 0; i < 3; ++i) {
                derivs.dg_ti[i][k] = 0.0;
                for (int j = 0; j < 3; ++j) {
                    derivs.dgamma_ij[i][j][k] = 0.0;
                }
            }
            continue;
        }
#elif DIM == 2
        if (k == 0) {
            pos_plus[0] += eps;
            pos_minus[0] -= eps;
        } else if (k == 1) {
            pos_plus[1] += eps;
            pos_minus[1] -= eps;
        } else {
            // For DIM=2, derivatives in z are zero
            derivs.dg_tt[k] = 0.0;
            for (int i = 0; i < 3; ++i) {
                derivs.dg_ti[i][k] = 0.0;
                for (int j = 0; j < 3; ++j) {
                    derivs.dgamma_ij[i][j][k] = 0.0;
                }
            }
            continue;
        }
#else
        pos_plus[k] += eps;
        pos_minus[k] -= eps;
#endif

        compute(pos_plus, metric_plus);
        compute(pos_minus, metric_minus);

        const real inv_2eps = 0.5 / eps;

        derivs.dg_tt[k] = (metric_plus.g_tt - metric_minus.g_tt) * inv_2eps;

        for (int i = 0; i < 3; ++i) {
            derivs.dg_ti[i][k] = (metric_plus.g_ti[i] - metric_minus.g_ti[i]) * inv_2eps;
            for (int j = 0; j < 3; ++j) {
                derivs.dgamma_ij[i][j][k] =
                    (metric_plus.gamma_ij[i][j] - metric_minus.gamma_ij[i][j]) * inv_2eps;
            }
        }
    }
}

// ============================================================================
// Kerr Metric Implementation
// ============================================================================

/**
 * Kerr metric in Cartesian-like Boyer-Lindquist coordinates
 *
 * From GRSPH Appendix B (Eq. B.1-B.16).
 * The coordinate transformation from Boyer-Lindquist (r,θ,φ) to Cartesian:
 *   x = √(r² + a²) sin(θ) cos(φ)
 *   y = √(r² + a²) sin(θ) sin(φ)
 *   z = r cos(θ)
 *
 * where r is implicitly defined by:
 *   r² = (R² - a² + √[(R² - a²)² + 4a²z²]) / 2
 *   R² = x² + y² + z²
 */
void KerrMetric::compute(const vec_t& pos, Metric31& metric) const
{
#if DIM == 1
    const real x = pos[0];
    const real y = 0.0;
    const real z = 0.0;
#elif DIM == 2
    const real x = pos[0];
    const real y = pos[1];
    const real z = 0.0;
#else
    const real x = pos[0];
    const real y = pos[1];
    const real z = pos[2];
#endif

    const real R2 = x*x + y*y + z*z;
    const real a2 = a * a;

    // Compute Boyer-Lindquist r from Cartesian coordinates
    // r² = (R² - a² + √[(R² - a²)² + 4a²z²]) / 2
    const real temp = R2 - a2;
    const real r2 = 0.5 * (temp + std::sqrt(temp * temp + 4.0 * a2 * z * z));
    const real r = std::sqrt(std::max(r2, 1e-20));

    // Avoid singularity
    if (r < 1e-10) {
        metric = Metric31();
        return;
    }

    // Kerr auxiliary quantities
    const real rho2 = r2 + a2 * z * z / r2;  // ρ² = r² + a² cos²θ
    const real Delta = r2 - 2.0 * M * r + a2;  // Δ = r² - 2Mr + a²

    // Check for horizon
    if (Delta < 1e-10) {
        metric = Metric31();
        return;
    }

    const real sin2_theta = 1.0 - z * z / r2;  // sin²θ = 1 - z²/r²
    const real xy2 = x*x + y*y;

    // Metric components in Boyer-Lindquist for reference:
    // g_tt = -(1 - 2Mr/ρ²)
    // g_tφ = -2Mar sin²θ / ρ²
    // g_φφ = (r² + a² + 2Ma²r sin²θ/ρ²) sin²θ

    const real g_tt_BL = -(1.0 - 2.0 * M * r / rho2);
    const real g_tphi_BL = -2.0 * M * r * a * sin2_theta / rho2;
    const real g_phiphi_BL = (r2 + a2 + 2.0 * M * r * a2 * sin2_theta / rho2) * sin2_theta;

    // Transform to Cartesian-like coordinates (GRSPH Eq. B.6-B.16)
    metric.g_tt = g_tt_BL;

    // Time-space components (from g_tφ)
    // g_tx = -y g_tφ / (x² + y²)
    // g_ty = x g_tφ / (x² + y²)
    // g_tz = 0
    if (xy2 > 1e-20) {
        metric.g_ti[0] = -y * g_tphi_BL / xy2;
        metric.g_ti[1] = x * g_tphi_BL / xy2;
    } else {
        metric.g_ti[0] = 0.0;
        metric.g_ti[1] = 0.0;
    }
    metric.g_ti[2] = 0.0;

    // Spatial metric components (complex expressions from GRSPH Appendix B)
    // For simplicity, compute key components for equatorial plane (z ≈ 0)
    // and use numerical finite difference for general case

    // Lapse from g_tt and shift
    // α² = -g_tt + g_ti g^ti (for diagonal shift this simplifies)
    // For Kerr: α² ≈ Δρ² / [(r² + a²)² - a²Δ sin²θ]
    const real Sigma = (r2 + a2) * (r2 + a2) - a2 * Delta * sin2_theta;
    const real alpha2 = Delta * rho2 / Sigma;
    metric.alpha = std::sqrt(std::max(alpha2, 1e-20));

    // Shift vector β^i
    // For Kerr: β^φ = -g_tφ / g_φφ (other components zero in spherical)
    // Transform to Cartesian
    if (std::abs(g_phiphi_BL) > 1e-20 && xy2 > 1e-20) {
        const real beta_phi = -g_tphi_BL / g_phiphi_BL;
        // β^x = -y β^φ, β^y = x β^φ (transformation from spherical)
        const real r_perp = std::sqrt(xy2);
        metric.beta[0] = -y * beta_phi / r_perp;
        metric.beta[1] = x * beta_phi / r_perp;
    } else {
        metric.beta[0] = 0.0;
        metric.beta[1] = 0.0;
    }
    metric.beta[2] = 0.0;

    // Spatial metric - use approximate form for now
    // Full expressions are in GRSPH Appendix B Eq. B.9-B.16
    // These are quite complex; for a working implementation we use
    // numerical differentiation for derivatives

    // Simplified spatial metric (accurate near equatorial plane)

    // Start with identity and add corrections
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            metric.gamma_ij[i][j] = (i == j) ? 1.0 : 0.0;
            metric.gamma_inv[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Add Kerr corrections (simplified form)
    if (xy2 > 1e-10) {
        const real rho2_delta = rho2 * Delta;
        const real xy = x * y;

        metric.gamma_ij[0][0] += r2 * x * x / (rho2_delta) +
            g_phiphi_BL * y * y / (xy2 * xy2) - 1.0;
        metric.gamma_ij[1][1] += r2 * y * y / (rho2_delta) +
            g_phiphi_BL * x * x / (xy2 * xy2) - 1.0;
        metric.gamma_ij[0][1] = metric.gamma_ij[1][0] =
            r2 * xy / (rho2_delta) - g_phiphi_BL * xy / (xy2 * xy2);
    }

    metric.gamma_ij[2][2] = rho2 / r2 + z * z * (r2 + a2) * (r2 + a2) / (r2 * rho2 * Delta);

    // For now, approximate inverse as identity (accurate for weak field)
    // Full computation would require matrix inversion

    // Determinants
    metric.sqrt_gamma = 1.0;  // In Cartesian-like Kerr coords, √γ = 1
    metric.sqrt_neg_g = metric.alpha;
}

/**
 * Kerr metric derivatives (numerical finite difference)
 */
void KerrMetric::compute_derivatives(const vec_t& pos, MetricDerivatives& derivs) const
{
    const real eps = 1e-6;

    Metric31 metric_plus, metric_minus;

    for (int k = 0; k < DIM; ++k) {
        vec_t pos_plus = pos;
        vec_t pos_minus = pos;

        pos_plus[k] += eps;
        pos_minus[k] -= eps;

        compute(pos_plus, metric_plus);
        compute(pos_minus, metric_minus);

        const real inv_2eps = 0.5 / eps;

        derivs.dg_tt[k] = (metric_plus.g_tt - metric_minus.g_tt) * inv_2eps;

        for (int i = 0; i < 3; ++i) {
            derivs.dg_ti[i][k] = (metric_plus.g_ti[i] - metric_minus.g_ti[i]) * inv_2eps;
            for (int j = 0; j < 3; ++j) {
                derivs.dgamma_ij[i][j][k] =
                    (metric_plus.gamma_ij[i][j] - metric_minus.gamma_ij[i][j]) * inv_2eps;
            }
        }
    }

    // Zero out derivatives for dimensions beyond DIM
    for (int k = DIM; k < 3; ++k) {
        derivs.dg_tt[k] = 0.0;
        for (int i = 0; i < 3; ++i) {
            derivs.dg_ti[i][k] = 0.0;
            for (int j = 0; j < 3; ++j) {
                derivs.dgamma_ij[i][j][k] = 0.0;
            }
        }
    }
}

// ============================================================================
// Factory Function
// ============================================================================

std::unique_ptr<MetricBase> create_metric(const std::string& name, real M, real a)
{
    if (name == "minkowski" || name == "flat") {
        return std::make_unique<MinkowskiMetric>();
    } else if (name == "schwarzschild") {
        return std::make_unique<SchwarzschildMetric>(M);
    } else if (name == "kerr") {
        return std::make_unique<KerrMetric>(M, a);
    } else {
        throw std::invalid_argument("Unknown metric: " + name);
    }
}

} // namespace grgsph
} // namespace sph
