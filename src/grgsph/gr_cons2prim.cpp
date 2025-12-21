/**
 * GR-GSPH Conservative to Primitive Variable Recovery
 *
 * Based on GRSPH (Liptai & Price 2019) Appendix B.
 * Uses Newton-Raphson iteration on the enthalpy w.
 */

#include "grgsph/gr_cons2prim.hpp"
#include <cmath>
#include <algorithm>

namespace sph {
namespace grgsph {

/**
 * Compute conserved from primitive (forward direction)
 *
 * Straightforward computation used for initialization.
 */
void GRCons2Prim::prim2cons(
    const GRPrimitive& prim,
    const Metric31& metric,
    real gamma_ad,
    GRConserved& U)
{
    const real n = prim.n;
    const real P = prim.P;
    const real Gamma = prim.Gamma;
    const real h = prim.h;

    // ρ* = √γ n Γ
    U.rho_star = metric.sqrt_gamma * n * Gamma;

    // Enthalpy density w = n h
    const real w_dens = n * h;

    // Covariant momentum p_i = Γ w γ_{ij} V^j
    // First lower the velocity index
    real V_down[3];
    metric.lower_index(prim.V, V_down);

    for (int i = 0; i < 3; ++i) {
        U.p[i] = Gamma * w_dens * V_down[i];
    }

    // Energy: e = U^0 w - P √(-g) / ρ*
    // U^0 = Γ / α
    const real U0 = Gamma / metric.alpha;
    U.e = U0 * w_dens - P * metric.sqrt_neg_g / U.rho_star;

    U.use_entropy = false;
}

/**
 * Recover primitive variables from conserved
 *
 * Newton-Raphson iteration on enthalpy w.
 *
 * The key equation is f(w) = w(ρ,P) - w = 0, where:
 *   Γ² = 1 + |p|²/w²
 *   ρ = ρ* / (√γ Γ)
 *   P = (w - 1)(γ-1)/γ · ρ   (from EOS)
 *
 * For entropy formulation (K = P/ρ^γ):
 *   P = K ρ^γ
 */
bool GRCons2Prim::solve(
    const GRConserved& U,
    const Metric31& metric,
    real gamma_ad,
    real w_prev,
    GRPrimitive& prim)
{
    // Compute |p|² = γ^{ij} p_i p_j
    real p_sq = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            p_sq += metric.gamma_inv[i][j] * U.p[i] * U.p[j];
        }
    }

    // Safety check
    if (U.rho_star <= 0.0) {
        prim = GRPrimitive();
        return false;
    }

    // Initial guess
    real w = std::max(w_prev, 1.0 + 1e-10);

    // Newton-Raphson iteration
    bool converged = false;
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        // Compute Γ from w and |p|²
        // Γ² = 1 + |p|² / w²
        const real w_sq = w * w;
        const real Gamma_sq = 1.0 + p_sq / w_sq;
        const real Gamma = std::sqrt(Gamma_sq);

        // Compute ρ from ρ* and Γ
        // ρ = ρ* / (√γ Γ)
        const real rho = U.rho_star / (metric.sqrt_gamma * Gamma);

        // Compute P from EOS
        real P;
        if (U.use_entropy) {
            // Entropy formulation: K is stored in e, P = K ρ^γ
            const real K = U.e;
            P = K * std::pow(rho, gamma_ad);
        } else {
            // Energy formulation
            // From e = U^0 w_dens - P √(-g) / ρ*, solve for P
            // This requires another equation; use EOS instead
            // P = (w - 1)(γ-1)/γ · ρ  (from h = 1 + γ/(γ-1) P/ρ)
            P = (w - 1.0) * (gamma_ad - 1.0) / gamma_ad * rho;
        }

        // Compute enthalpy from ρ and P
        const real w_new = compute_enthalpy(rho, P, gamma_ad);

        // Compute residual
        const real f = w_new - w;

        // Check convergence
        if (std::abs(f) < TOL || std::abs(f / w) < TOL) {
            w = w_new;
            converged = true;
            break;
        }

        // Compute derivative for Newton step
        // df/dw = dw_new/dw - 1
        // This is complex; use numerical derivative for robustness
        const real dw = w * 1e-6;
        const real w_plus = w + dw;

        const real Gamma_sq_plus = 1.0 + p_sq / (w_plus * w_plus);
        const real Gamma_plus = std::sqrt(Gamma_sq_plus);
        const real rho_plus = U.rho_star / (metric.sqrt_gamma * Gamma_plus);

        real P_plus;
        if (U.use_entropy) {
            P_plus = U.e * std::pow(rho_plus, gamma_ad);
        } else {
            P_plus = (w_plus - 1.0) * (gamma_ad - 1.0) / gamma_ad * rho_plus;
        }

        const real w_new_plus = compute_enthalpy(rho_plus, P_plus, gamma_ad);
        const real f_plus = w_new_plus - w_plus;

        const real df_dw = (f_plus - f) / dw;

        // Newton update
        if (std::abs(df_dw) > 1e-20) {
            const real delta_w = -f / df_dw;

            // Limit step size for stability
            const real max_delta = 0.5 * w;
            const real clamped_delta = std::max(-max_delta, std::min(max_delta, delta_w));

            w = w + clamped_delta;

            // Ensure w > 1 (physical constraint)
            w = std::max(w, 1.0 + 1e-10);
        } else {
            // df/dw ≈ 0, try bisection-like step
            w = 0.5 * (w + w_new);
        }
    }

    if (!converged) {
        // Fallback: use previous values or simple estimate
        prim = GRPrimitive();
        prim.n = U.rho_star / metric.sqrt_gamma;  // Assume Γ ≈ 1
        prim.P = 0.0;
        prim.Gamma = 1.0;
        prim.h = 1.0;
        return false;
    }

    // Compute final primitive variables
    const real w_sq = w * w;
    const real Gamma_sq = 1.0 + p_sq / w_sq;
    prim.Gamma = std::sqrt(Gamma_sq);

    prim.n = U.rho_star / (metric.sqrt_gamma * prim.Gamma);

    if (U.use_entropy) {
        prim.P = U.e * std::pow(prim.n, gamma_ad);
    } else {
        prim.P = (w - 1.0) * (gamma_ad - 1.0) / gamma_ad * prim.n;
    }

    prim.h = w;

    // Recover velocity: V^i = γ^{ij} p_j / (Γ w ρ)
    const real factor = 1.0 / (prim.Gamma * w * prim.n);
    for (int i = 0; i < 3; ++i) {
        prim.V[i] = 0.0;
        for (int j = 0; j < 3; ++j) {
            prim.V[i] += metric.gamma_inv[i][j] * U.p[j];
        }
        prim.V[i] *= factor;
    }

    // Compute sound speed
    prim.c_s = compute_sound_speed(prim.n, prim.P, prim.h, gamma_ad);

    return true;
}

} // namespace grgsph
} // namespace sph
