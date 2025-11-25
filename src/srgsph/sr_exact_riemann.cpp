#include "srgsph/sr_exact_riemann.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

namespace sph {
namespace srgsph {
namespace riemann {

/**
 * Calculate velocity behind shock wave
 * Based on Python solve_shock() (lines 31-154 in riemann_solver.py)
 */
static bool solve_shock(
    real p_target, const RiemannState& state, real vta, real gamma, real c_speed,
    bool is_left_wave, real& vxb, real& rhob, real& hb, real& vtb)
{
    if (p_target <= 0.0) return false;

    const real pa = std::max(state.P, 1e-10);
    const real pb = p_target;
    const real rhoa = state.n;
    const real va = state.v;

    // Compute enthalpy
    const real ha = 1.0 + (gamma / (gamma - 1.0)) * (pa / rhoa);

    // Quadratic for h_b (Eq 27 in Pons et al.)
    const real A = (gamma - 1.0) * (pa - pb) / (gamma * pb);
    const real B = ha * (pa - pb) / rhoa;

    const real qa = 1.0 + A;
    const real qb = -A;
    const real qc = B - ha * ha;

    const real delta = qb * qb - 4.0 * qa * qc;
    if (delta < 0.0) return false;

    hb = (-qb + std::sqrt(delta)) / (2.0 * qa);

    // Recover rho_b
    rhob = (gamma / (gamma - 1.0)) * pb / (hb - 1.0);

    // Mass flux j^2
    const real denom = (hb / rhob - ha / rhoa);
    if (std::abs(denom) < 1e-15) {
        vxb = va;
        return true;
    }

    const real j2 = -(pb - pa) / denom;
    if (j2 <= 0.0) {
        vxb = va;
        return true;
    }

    const real j = std::sqrt(j2);

    // Lorentz factor for state a - MUST include tangential velocity!
    // Pons et al. (2000): Use FULL Lorentz factor from total velocity
    const real v2a_total = va * va + vta * vta;  // Total velocity squared
    const real Wa = 1.0 / std::sqrt(1.0 - v2a_total);
    const real Da = rhoa * Wa;  // Lab-frame density D = ρ*W

    // Shock speed (Pons et al. 2000, Eq. 26) - Python line 108
    // NOTE: This is DIFFERENT from MM94! Uses (1-vx^2) factor for tangential velocity coupling
    const real term = j * j + Da * Da * (1.0 - va * va);
    const real sqrt_term = std::sqrt(term);
    const real denom_vs = Da * Da + j * j;

    real Vs;
    if (is_left_wave) {
        Vs = (Da * Da * va - j * sqrt_term) / denom_vs;
    } else {
        Vs = (Da * Da * va + j * sqrt_term) / denom_vs;
    }

    const real Ws = 1.0 / std::sqrt(1.0 - Vs * Vs);
    const real j_signed = is_left_wave ? -j : j;

    if (std::abs(j_signed) < 1e-15) {
        vxb = va;
        return true;
    }

    // Velocity behind shock
    const real num = ha * Wa * va + Ws * (pb - pa) / j_signed;
    const real den = ha * Wa + (pb - pa) * (Ws * va / j_signed + 1.0 / Da);

    if (std::abs(den) < 1e-15) return false;

    vxb = num / den;

    // Tangential velocity (Eq 25 in Pons et al.)
    // Python lines 148-149
    if (vta > 1e-10) {
        const real factor_num = 1.0 - vxb * vxb;
        const real factor_den = hb * hb + (ha * Wa * vta) * (ha * Wa * vta);
        if (factor_num > 0.0 && factor_den > 0.0) {
            vtb = vta * ha * Wa * std::sqrt(factor_num / factor_den);
        } else {
            vtb = 0.0;
        }
    } else {
        vtb = 0.0;
    }

    return true;
}

/**
 * Calculate velocity behind rarefaction wave
 * Based on Python solve_rarefaction_analytical() (lines 156-201 in riemann_solver.py)
 */
static bool solve_rarefaction_zero_tangent(
    real p_target, const RiemannState& state, real gamma,
    bool is_left_wave, real& vxb, real& rhob, real& hb, real& vtb)
{
    const real pa = std::max(state.P, 1e-10);
    const real rhoa = state.n;

    // Isentropic relation
    const real const_entropy = pa / std::pow(rhoa, gamma);
    rhob = std::pow(p_target / const_entropy, 1.0 / gamma);

    // Sound speeds
    const real ua = pa / ((gamma - 1.0) * rhoa);
    const real ha_calc = 1.0 + ua + pa / rhoa;
    const real csa = std::sqrt(gamma * pa / (rhoa * ha_calc));

    const real ub = p_target / ((gamma - 1.0) * rhob);
    hb = 1.0 + ub + p_target / rhob;
    const real csb = std::sqrt(gamma * p_target / (rhob * hb));

    const real sign = is_left_wave ? 1.0 : -1.0;
    const real sqrt_g_minus_1 = std::sqrt(gamma - 1.0);

    const real term_v = (1.0 + state.v) / (1.0 - state.v);
    const real term_ca = (sqrt_g_minus_1 + csa) / (sqrt_g_minus_1 - csa);
    const real term_cb = (sqrt_g_minus_1 + csb) / (sqrt_g_minus_1 - csb);

    const real exponent = sign * 2.0 / sqrt_g_minus_1;
    const real base = term_ca / term_cb;

    const real A = term_v * std::pow(base, exponent);
    vxb = (A - 1.0) / (A + 1.0);

    vtb = 0.0;

    return true;
}

static bool solve_rarefaction(
    real p_target, const RiemannState& state, real vta, real gamma, real c_speed,
    bool is_left_wave, real& vxb, real& rhob, real& hb, real& vtb)
{
    const real pa = std::max(state.P, 1e-10);
    if (vta < 1e-6) {
        return solve_rarefaction_zero_tangent(p_target, state, gamma, is_left_wave, vxb, rhob, hb, vtb);
    }

    if (p_target > pa) {
        return false; // Should be handled by shock solver
    }

    const real rhoa = state.n;
    const real const_entropy = pa / std::pow(rhoa, gamma);
    const real ha = 1.0 + (gamma / (gamma - 1.0)) * (pa / rhoa);

    const real v2_total_a = state.v * state.v + vta * vta;
    const real v2_clamped = std::min(v2_total_a, 1.0 - 1e-12);
    const real Wa = 1.0 / std::sqrt(1.0 - v2_clamped);
    const real Kt = ha * Wa * vta;

    const real sign = is_left_wave ? -1.0 : 1.0;        // Eq. 17 in Pons et al.
    const real xi_sign = is_left_wave ? -1.0 : 1.0;      // Eq. 14 sign convention

    const real pressure_span = std::abs(pa - p_target);
    if (pressure_span < 1e-12) {
        vxb = state.v;
        rhob = rhoa;
        hb = ha;
        vtb = vta;
        return true;
    }

    const int raw_steps = static_cast<int>(pressure_span / std::max(pa * 1e-3, 1e-6)) + 64;
    const int steps = std::max(64, std::min(2048, raw_steps));

    auto deriv = [&](real p, real vx) -> real {
        const real p_safe = std::max(p, 1e-12);
        const real rho = std::pow(p_safe / const_entropy, 1.0 / gamma);
        const real h = 1.0 + (gamma / (gamma - 1.0)) * (p_safe / rho);

        const real K_over_h = Kt / h;
        const real denom_v = std::max(1e-12, 1.0 - vx * vx);
        const real W2 = (1.0 + K_over_h * K_over_h) / denom_v;
        if (W2 <= 0.0) return std::numeric_limits<real>::quiet_NaN();
        const real W = std::sqrt(W2);
        const real vt = Kt / (h * W);

        const real cs2 = std::max(0.0, gamma * p_safe / (rho * h));
        const real cs = std::sqrt(cs2);

        const real v2 = vx * vx + vt * vt;
        const real term1 = vx * (1.0 - cs2);
        real term2_inner = (1.0 - v2) * (1.0 - v2 * cs2 - vx * vx * (1.0 - cs2));
        term2_inner = std::max(term2_inner, 0.0);
        const real term2 = cs * std::sqrt(term2_inner);
        const real denom = std::max(1e-12, 1.0 - v2 * cs2);

        const real xi = (term1 + xi_sign * term2) / denom;
        const real one_minus_xiv = std::max(1e-8, 1.0 - xi * vx);
        const real g = (vt * vt * (xi * xi - 1.0)) / (one_minus_xiv * one_minus_xiv);

        const real denom_deriv = rho * h * W2 * cs * std::sqrt(1.0 + g);
        if (denom_deriv <= 0.0 || !std::isfinite(denom_deriv)) {
            return std::numeric_limits<real>::quiet_NaN();
        }
        return sign / denom_deriv;
    };

    real p_curr = pa;
    real vx_curr = state.v;

    for (int step = 0; step < steps; ++step) {
        const real p_next = pa + (p_target - pa) * static_cast<real>(step + 1) / steps;
        const real dp = p_next - p_curr;

        const real k1 = deriv(p_curr, vx_curr);
        if (!std::isfinite(k1)) return false;
        const real k2 = deriv(p_curr + 0.5 * dp, vx_curr + 0.5 * dp * k1);
        if (!std::isfinite(k2)) return false;
        const real k3 = deriv(p_curr + 0.5 * dp, vx_curr + 0.5 * dp * k2);
        if (!std::isfinite(k3)) return false;
        const real k4 = deriv(p_curr + dp, vx_curr + dp * k3);
        if (!std::isfinite(k4)) return false;

        vx_curr += (dp / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        p_curr = p_next;
    }

    vxb = vx_curr;

    const real pb = p_target;
    rhob = std::pow(pb / const_entropy, 1.0 / gamma);
    hb = 1.0 + (gamma / (gamma - 1.0)) * (pb / rhob);

    const real Kt_hb = Kt / hb;
    const real W2_b = (1.0 + Kt_hb * Kt_hb) / std::max(1e-12, 1.0 - vxb * vxb);
    if (W2_b <= 0.0) return false;
    const real Wb = std::sqrt(W2_b);
    vtb = Kt / (hb * Wb);

    return true;
}

/**
 * Wave curve function: velocity as function of pressure
 */
static real wave_curve(real p, const RiemannState& state, real vta, real gamma, real c_speed, bool is_left_wave)
{
    real vxb, rhob, hb, vtb;
    bool success;

    if (p > state.P) {
        success = solve_shock(p, state, vta, gamma, c_speed, is_left_wave, vxb, rhob, hb, vtb);
    } else {
        success = solve_rarefaction(p, state, vta, gamma, c_speed, is_left_wave, vxb, rhob, hb, vtb);
    }

    return success ? vxb : std::numeric_limits<real>::quiet_NaN();
}

/**
 * Exact Riemann solver using Newton-Raphson iteration
 * Based on Python solve_riemann() (lines 346-443)
 */
bool exact_riemann_solver(
    const RiemannState& left,
    const RiemannState& right,
    real vt_left,
    real vt_right,
    real gamma_c,
    real c_speed,
    real& P_star,
    real& v_star,
    real& vt_star,
    int max_iter,
    real tol)
{
    // Check for identical states
    if (std::abs(left.P - right.P) < 1e-6 &&
        std::abs(left.n - right.n) < 1e-6 &&
        std::abs(left.v - right.v) < 1e-6) {
        P_star = left.P;
        v_star = left.v;
        vt_star = vt_left;
        return true;
    }

    // Function to minimize: f(p) = v_left(p) - v_right(p)
    auto func = [&](real p) -> real {
        const real vl = wave_curve(p, left, vt_left, gamma_c, c_speed, true);
        const real vr = wave_curve(p, right, vt_right, gamma_c, c_speed, false);
        if (!std::isfinite(vl) || !std::isfinite(vr)) {
            return std::numeric_limits<real>::quiet_NaN();
        }
        return vl - vr;
    };

    // Bracket the root
    real p_min = std::min(left.P, right.P) * 0.001;
    real p_max = std::max(left.P, right.P) * 1000.0;

    real f_min = func(p_min);
    real f_max = func(p_max);

    auto is_valid = [](real value) {
        return std::isfinite(value);
    };

    if (!is_valid(f_min) || !is_valid(f_max)) {
        return false;
    }

    // Check if we can bracket
    if (f_min * f_max > 0.0) {
        // Try expanding bracket
        p_min *= 0.01;
        p_max *= 100.0;
        f_min = func(p_min);
        f_max = func(p_max);

        if (!is_valid(f_min) || !is_valid(f_max) || f_min * f_max > 0.0) {
            return false; // Cannot bracket
        }
    }

    // Brent's method (simplified bisection for robustness)
    real p_lower = p_min;
    real p_upper = p_max;
    P_star = 0.5 * (p_lower + p_upper);

    for (int iter = 0; iter < max_iter; ++iter) {
        const real f = func(P_star);
        if (!is_valid(f)) {
            return false;
        }

        if (std::abs(f) < tol || std::abs(p_upper - p_lower) < tol * P_star) {
            break;
        }

        if (f * f_min < 0.0) {
            p_upper = P_star;
            f_max = f;
        } else {
            p_lower = P_star;
            f_min = f;
        }

        P_star = 0.5 * (p_lower + p_upper);
    }

    // Compute v_star and vt_star (Python lines 420-443)
    real vl_star, rhol, hl, vtl_star;
    real vr_star, rhor, hr, vtr_star;

    // Get full star states from both sides
    if (P_star > left.P) {
        solve_shock(P_star, left, vt_left, gamma_c, c_speed, true, vl_star, rhol, hl, vtl_star);
    } else {
        solve_rarefaction(P_star, left, vt_left, gamma_c, c_speed, true, vl_star, rhol, hl, vtl_star);
    }

    if (P_star > right.P) {
        solve_shock(P_star, right, vt_right, gamma_c, c_speed, false, vr_star, rhor, hr, vtr_star);
    } else {
        solve_rarefaction(P_star, right, vt_right, gamma_c, c_speed, false, vr_star, rhor, hr, vtr_star);
    }

    // Average normal velocity
    v_star = 0.5 * (vl_star + vr_star);

    // Select tangential velocity based on direction (upwinding)
    // Python lines 426-443: if v_star > 0, use left state; else use right state
    if (v_star > 0.0) {
        vt_star = vtl_star;
    } else {
        vt_star = vtr_star;
    }

    return true;
}

/**
 * HLLC Riemann solver (approximate but robust)
 * Based on standard HLLC formulation for relativistic hydro
 */
void hllc_riemann_solver(
    const RiemannState& left,
    const RiemannState& right,
    real vt_left,
    real vt_right,
    real gamma_c,
    real c_speed,
    real& P_star,
    real& v_star,
    real& vt_star)
{
    // Compute wave speeds using simple estimates
    const real cs_l = left.c_s;
    const real cs_r = right.c_s;

    // Left and right wave speeds (simplified)
    const real S_L = std::min(left.v - cs_l, right.v - cs_r);
    const real S_R = std::max(left.v + cs_l, right.v + cs_r);

    // Contact wave speed (star region velocity)
    v_star = (right.P - left.P + left.n * left.v * (S_L - left.v) -
              right.n * right.v * (S_R - right.v)) /
             (left.n * (S_L - left.v) - right.n * (S_R - right.v));

    // Clamp velocity
    if (v_star >= 1.0) v_star = 0.99;
    if (v_star <= -1.0) v_star = -0.99;

    // Star region pressure (HLLC estimate)
    P_star = left.P + left.n * (left.v - S_L) * (left.v - v_star);

    // Safety bounds
    if (P_star < 0.0) P_star = 0.5 * (left.P + right.P);

    // Tangential velocity (upwinding based on contact wave)
    if (v_star > 0.0) {
        vt_star = vt_left;
    } else {
        vt_star = vt_right;
    }
}

} // namespace riemann
} // namespace srgsph
} // namespace sph
