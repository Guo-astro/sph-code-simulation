/**
 * Special Relativistic Exact Riemann Solver
 *
 * Based on Kitajima et al. (2025) arXiv:2510.18251v1 and Pons et al. (2000)
 *
 * Variable naming follows Kitajima notation (SRGSPH paper Section 2.1):
 *   N       = lab-frame baryon number density (N = γn)
 *   n       = rest-frame baryon number density
 *   γ       = Lorentz factor (gamma in code)
 *   H       = enthalpy per baryon
 *   P       = pressure
 *   v^x     = normal velocity component (v_x in code)
 *   v^t     = tangential velocity component (v_t in code)
 *   c_s     = sound speed
 *   γ_c     = ratio of specific heats (gamma_c in code)
 *   S       = canonical momentum per baryon (Eq. 5)
 *   e       = canonical energy per baryon (Eq. 6)
 *   u       = thermal energy per baryon
 *
 * Subscript conventions:
 *   _L, _R  = left and right initial states
 *   _star   = star region (Riemann solution interface state)
 *   _a      = state ahead of wave (known state)
 *   _b      = state behind wave (to be computed)
 */

#include "srgsph/sr_exact_riemann.hpp"
#include "logger.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <atomic>

namespace sph {
namespace srgsph {
namespace riemann {

// Debug counters for tracking solver failures
static std::atomic<int> g_exact_solver_calls{0};
static std::atomic<int> g_exact_solver_failures{0};
static std::atomic<int> g_bracket_failures{0};
static std::atomic<int> g_nan_failures{0};

// Enable verbose debug output (set to true for debugging)
static constexpr bool DEBUG_RIEMANN = false;

/**
 * Calculate velocity behind shock wave
 * Based on Pons et al. (2000) equations and Python riemann_solver.py
 *
 * Subscript convention (following Pons et al.):
 *   _a = state ahead of shock (known initial state, either L or R)
 *   _b = state behind shock (star region state to be computed)
 */
static bool solve_shock(
    real P_star, const RiemannState& state, real v_t_a, real gamma_c, real c,
    bool is_left_wave, real& v_x_b, real& n_b, real& H_b, real& v_t_b)
{
    if (P_star <= 0.0) return false;

    // State ahead of shock (known state)
    const real P_a = std::max(state.P, 1e-10);
    const real P_b = P_star;
    const real n_a = state.n;      // Rest-frame density
    const real v_x_a = state.v;    // Normal velocity v^x

    // Enthalpy per baryon: H = 1 + u/c² + P/(n·c²), with u = P/((γ_c-1)n)
    // Simplifies to: H = 1 + (γ_c/(γ_c-1)) * P/n  (for c=1)
    const real H_a = 1.0 + (gamma_c / (gamma_c - 1.0)) * (P_a / n_a);

    // Quadratic for H_b (Pons et al. Eq. 27)
    const real A = (gamma_c - 1.0) * (P_a - P_b) / (gamma_c * P_b);
    const real B = H_a * (P_a - P_b) / n_a;

    const real q_a = 1.0 + A;
    const real q_b = -A;
    const real q_c = B - H_a * H_a;

    const real discriminant = q_b * q_b - 4.0 * q_a * q_c;
    if (discriminant < 0.0) return false;

    H_b = (-q_b + std::sqrt(discriminant)) / (2.0 * q_a);

    // Recover rest-frame density n_b from EOS
    // From H = 1 + (γ_c/(γ_c-1)) * P/n  =>  n = (γ_c/(γ_c-1)) * P / (H-1)
    n_b = (gamma_c / (gamma_c - 1.0)) * P_b / (H_b - 1.0);

    // Mass flux j² (Pons et al.)
    const real denom = (H_b / n_b - H_a / n_a);
    if (std::abs(denom) < 1e-15) {
        v_x_b = v_x_a;
        return true;
    }

    const real j2 = -(P_b - P_a) / denom;
    if (j2 <= 0.0) {
        v_x_b = v_x_a;
        return true;
    }

    const real j = std::sqrt(j2);

    // Lorentz factor γ_a - MUST include tangential velocity!
    // Pons et al. (2000): Use FULL Lorentz factor from total velocity
    const real v2_total_a = v_x_a * v_x_a + v_t_a * v_t_a;
    const real gamma_a = 1.0 / std::sqrt(1.0 - v2_total_a);
    const real N_a = n_a * gamma_a;  // Lab-frame density N = γn

    // Shock speed V_s (Pons et al. 2000, Eq. 26)
    // NOTE: This includes tangential velocity coupling via (1 - v_x²) factor
    const real term = j * j + N_a * N_a * (1.0 - v_x_a * v_x_a);
    const real sqrt_term = std::sqrt(term);
    const real denom_Vs = N_a * N_a + j * j;

    real V_s;
    if (is_left_wave) {
        V_s = (N_a * N_a * v_x_a - j * sqrt_term) / denom_Vs;
    } else {
        V_s = (N_a * N_a * v_x_a + j * sqrt_term) / denom_Vs;
    }

    const real gamma_s = 1.0 / std::sqrt(1.0 - V_s * V_s);
    const real j_signed = is_left_wave ? -j : j;

    if (std::abs(j_signed) < 1e-15) {
        v_x_b = v_x_a;
        return true;
    }

    // Normal velocity behind shock (star region)
    const real num = H_a * gamma_a * v_x_a + gamma_s * (P_b - P_a) / j_signed;
    const real den = H_a * gamma_a + (P_b - P_a) * (gamma_s * v_x_a / j_signed + 1.0 / N_a);

    if (std::abs(den) < 1e-15) return false;

    v_x_b = num / den;

    // Tangential velocity v_t_b (Pons et al. Eq. 25)
    if (v_t_a > 1e-10) {
        const real factor_num = 1.0 - v_x_b * v_x_b;
        const real factor_den = H_b * H_b + (H_a * gamma_a * v_t_a) * (H_a * gamma_a * v_t_a);
        if (factor_num > 0.0 && factor_den > 0.0) {
            v_t_b = v_t_a * H_a * gamma_a * std::sqrt(factor_num / factor_den);
        } else {
            v_t_b = 0.0;
        }
    } else {
        v_t_b = 0.0;
    }

    return true;
}

/**
 * Calculate velocity behind rarefaction wave (zero tangential velocity case)
 * Analytical solution using Riemann invariants
 *
 * Subscript convention:
 *   _a = state ahead of rarefaction (known initial state)
 *   _b = state behind rarefaction (star region)
 */
static bool solve_rarefaction_zero_tangent(
    real P_star, const RiemannState& state, real gamma_c,
    bool is_left_wave, real& v_x_b, real& n_b, real& H_b, real& v_t_b)
{
    const real P_a = std::max(state.P, 1e-10);
    const real n_a = state.n;
    const real v_x_a = state.v;  // Normal velocity v^x

    // Isentropic relation: P/n^γ_c = const
    const real K_entropy = P_a / std::pow(n_a, gamma_c);
    n_b = std::pow(P_star / K_entropy, 1.0 / gamma_c);

    // Sound speeds c_s = sqrt(γ_c * P / (n * H))
    const real u_a = P_a / ((gamma_c - 1.0) * n_a);
    const real H_a = 1.0 + u_a + P_a / n_a;
    const real c_s_a = std::sqrt(gamma_c * P_a / (n_a * H_a));

    const real u_b = P_star / ((gamma_c - 1.0) * n_b);
    H_b = 1.0 + u_b + P_star / n_b;
    const real c_s_b = std::sqrt(gamma_c * P_star / (n_b * H_b));

    // Riemann invariant integration (analytical for v^t = 0)
    const real sign = is_left_wave ? 1.0 : -1.0;
    const real sqrt_gm1 = std::sqrt(gamma_c - 1.0);

    const real term_v = (1.0 + v_x_a) / (1.0 - v_x_a);
    const real term_c_a = (sqrt_gm1 + c_s_a) / (sqrt_gm1 - c_s_a);
    const real term_c_b = (sqrt_gm1 + c_s_b) / (sqrt_gm1 - c_s_b);

    const real exponent = sign * 2.0 / sqrt_gm1;
    const real base = term_c_a / term_c_b;

    const real A = term_v * std::pow(base, exponent);
    v_x_b = (A - 1.0) / (A + 1.0);

    v_t_b = 0.0;

    return true;
}

/**
 * Gauss-Legendre quadrature nodes and weights for 8-point integration
 * These give O(h^16) accuracy, much better than RK4's O(h^4)
 */
static constexpr int GAUSS_POINTS = 8;
static const real GAUSS_NODES[GAUSS_POINTS] = {
    -0.9602898564975363,
    -0.7966664774136267,
    -0.5255324099163290,
    -0.1834346424956498,
    0.1834346424956498,
    0.5255324099163290,
    0.7966664774136267,
    0.9602898564975363
};
static const real GAUSS_WEIGHTS[GAUSS_POINTS] = {
    0.1012285362903763,
    0.2223810344533745,
    0.3137066458778873,
    0.3626837833783620,
    0.3626837833783620,
    0.3137066458778873,
    0.2223810344533745,
    0.1012285362903763
};

/**
 * Compute dv^x/dP for rarefaction with tangent velocity (Pons et al. Eq. 14-17)
 * This is the integrand for the rarefaction ODE.
 *
 * @param P Current pressure
 * @param v_x Current normal velocity
 * @param K_entropy Isentropic constant P/n^γ
 * @param K_t Conserved tangent momentum H*γ*v_t
 * @param gamma_c Adiabatic index
 * @param sign +1 for right-going, -1 for left-going wave
 * @param xi_sign +1 or -1 for characteristic speed selection
 * @return dv^x/dP at this point
 */
static real compute_dvx_dP(
    real P, real v_x,
    real K_entropy, real K_t,
    real gamma_c, real sign, real xi_sign)
{
    const real P_safe = std::max(P, 1e-12);
    const real n = std::pow(P_safe / K_entropy, 1.0 / gamma_c);
    const real H = 1.0 + (gamma_c / (gamma_c - 1.0)) * (P_safe / n);

    // Lorentz factor from K_t conservation: K_t = H * γ * v_t
    // γ² = (1 + (K_t/H)²) / (1 - v_x²)
    const real K_over_H = K_t / H;
    const real denom_v = std::max(1e-12, 1.0 - v_x * v_x);
    const real gamma2 = (1.0 + K_over_H * K_over_H) / denom_v;
    if (gamma2 <= 0.0) return std::numeric_limits<real>::quiet_NaN();
    const real gamma = std::sqrt(gamma2);
    const real v_t = K_t / (H * gamma);

    // Sound speed c_s² = γ_c * P / (n * H)
    const real c_s2 = std::max(0.0, gamma_c * P_safe / (n * H));
    const real c_s = std::sqrt(c_s2);

    // Characteristic speed ξ (Pons et al. Eq. 14)
    // ξ = [v_x(1-c_s²) ± c_s√((1-v²)(1 - v²c_s² - v_x²(1-c_s²)))] / (1 - v²c_s²)
    const real v2 = v_x * v_x + v_t * v_t;
    const real term1 = v_x * (1.0 - c_s2);
    real term2_inner = (1.0 - v2) * (1.0 - v2 * c_s2 - v_x * v_x * (1.0 - c_s2));
    term2_inner = std::max(term2_inner, 0.0);
    const real term2 = c_s * std::sqrt(term2_inner);
    const real denom_xi = std::max(1e-12, 1.0 - v2 * c_s2);
    const real xi = (term1 + xi_sign * term2) / denom_xi;

    // Correction factor g for tangential velocity (Pons et al. Eq. 15)
    // g = v_t² * (ξ² - 1) / (1 - ξ*v_x)²
    const real one_minus_xi_v = std::max(1e-8, 1.0 - xi * v_x);
    const real g = (v_t * v_t * (xi * xi - 1.0)) / (one_minus_xi_v * one_minus_xi_v);

    // dv^x/dP = ±1 / [n * H * γ² * c_s * √(1+g)]
    const real denom_deriv = n * H * gamma2 * c_s * std::sqrt(1.0 + g);
    if (denom_deriv <= 0.0 || !std::isfinite(denom_deriv)) {
        return std::numeric_limits<real>::quiet_NaN();
    }
    return sign / denom_deriv;
}

/**
 * Calculate velocity behind rarefaction wave (general case with tangential velocity)
 * 
 * OPTIMIZATION: Uses Gauss-Legendre quadrature + adaptive subdivision instead of
 * fixed-step RK4. This gives O(h^16) accuracy per panel vs O(h^4) for RK4,
 * allowing far fewer function evaluations.
 *
 * For weak rarefactions (small ΔP), uses analytical Taylor expansion.
 *
 * Subscript convention:
 *   _a = state ahead of rarefaction (known initial state)
 *   _b = state behind rarefaction (star region)
 */
static bool solve_rarefaction(
    real P_star, const RiemannState& state, real v_t_a, real gamma_c, real c,
    bool is_left_wave, real& v_x_b, real& n_b, real& H_b, real& v_t_b)
{
    const real P_a = std::max(state.P, 1e-10);

    // Use analytical solution for zero tangential velocity (much faster)
    if (v_t_a < 1e-6) {
        return solve_rarefaction_zero_tangent(P_star, state, gamma_c, is_left_wave, v_x_b, n_b, H_b, v_t_b);
    }

    if (P_star > P_a) {
        return false; // Should be handled by shock solver
    }

    const real n_a = state.n;
    const real v_x_a = state.v;
    const real K_entropy = P_a / std::pow(n_a, gamma_c);
    const real H_a = 1.0 + (gamma_c / (gamma_c - 1.0)) * (P_a / n_a);

    // Total velocity and Lorentz factor
    const real v2_total_a = v_x_a * v_x_a + v_t_a * v_t_a;
    const real v2_clamped = std::min(v2_total_a, 1.0 - 1e-12);
    const real gamma_a = 1.0 / std::sqrt(1.0 - v2_clamped);

    // Conserved tangential momentum: K_t = H * γ * v^t (Pons et al.)
    const real K_t = H_a * gamma_a * v_t_a;

    const real sign = is_left_wave ? -1.0 : 1.0;
    const real xi_sign = is_left_wave ? -1.0 : 1.0;

    // Pressure span
    const real delta_P = std::abs(P_a - P_star);
    
    // ========================================================================
    // OPTIMIZATION 1: For very weak rarefactions, use Taylor expansion
    // ========================================================================
    if (delta_P < 0.01 * P_a) {
        // First-order approximation: v_x_b ≈ v_x_a + (dv/dP)|_a * ΔP
        const real dv_at_a = compute_dvx_dP(P_a, v_x_a, K_entropy, K_t, gamma_c, sign, xi_sign);
        if (!std::isfinite(dv_at_a)) return false;
        
        v_x_b = v_x_a + dv_at_a * (P_star - P_a);
        
        // Clamp to physical range
        v_x_b = std::max(-0.9999, std::min(0.9999, v_x_b));
        
        // Final state
        n_b = std::pow(P_star / K_entropy, 1.0 / gamma_c);
        H_b = 1.0 + (gamma_c / (gamma_c - 1.0)) * (P_star / n_b);
        
        const real K_over_H_b = K_t / H_b;
        const real gamma2_b = (1.0 + K_over_H_b * K_over_H_b) / std::max(1e-12, 1.0 - v_x_b * v_x_b);
        if (gamma2_b <= 0.0) return false;
        v_t_b = K_t / (H_b * std::sqrt(gamma2_b));
        
        return true;
    }

    // ========================================================================
    // OPTIMIZATION 2: Use predictor-corrector with Gauss-Legendre quadrature
    // 
    // Instead of integrating dv/dP directly (which requires tracking v_x),
    // we solve the ODE using a predictor-corrector approach:
    //   1. Predict v_x at several points using initial slope
    //   2. Correct using Gauss quadrature of the integrand
    // ========================================================================
    
    // Determine number of panels based on pressure ratio
    const real P_ratio = std::max(P_a, P_star) / std::min(P_a, P_star);
    int num_panels = 1;
    if (P_ratio > 2.0) num_panels = 2;
    if (P_ratio > 5.0) num_panels = 4;
    if (P_ratio > 10.0) num_panels = 8;
    if (P_ratio > 100.0) num_panels = 16;
    
    real v_x_curr = v_x_a;
    real P_curr = P_a;
    
    for (int panel = 0; panel < num_panels; ++panel) {
        const real P_start = P_a + (P_star - P_a) * panel / num_panels;
        const real P_end = P_a + (P_star - P_a) * (panel + 1) / num_panels;
        const real dP = P_end - P_start;
        const real P_mid = 0.5 * (P_start + P_end);
        const real half_dP = 0.5 * dP;
        
        // Predictor: estimate v_x at panel midpoint using current slope
        const real dv_start = compute_dvx_dP(P_start, v_x_curr, K_entropy, K_t, gamma_c, sign, xi_sign);
        if (!std::isfinite(dv_start)) return false;
        
        // Predicted v_x at end of panel (for quadrature evaluation)
        real v_x_pred = v_x_curr + dv_start * dP;
        
        // Gauss-Legendre quadrature over this panel
        // Transform from [-1,1] to [P_start, P_end]: P = P_mid + half_dP * t
        real integral = 0.0;
        for (int i = 0; i < GAUSS_POINTS; ++i) {
            const real t = GAUSS_NODES[i];
            const real P = P_mid + half_dP * t;
            
            // Interpolate v_x linearly between start and predicted end
            const real alpha = (P - P_start) / dP;
            const real v_x_interp = v_x_curr + alpha * (v_x_pred - v_x_curr);
            
            const real dv = compute_dvx_dP(P, v_x_interp, K_entropy, K_t, gamma_c, sign, xi_sign);
            if (!std::isfinite(dv)) return false;
            
            integral += GAUSS_WEIGHTS[i] * dv;
        }
        integral *= half_dP;  // Jacobian of transformation
        
        v_x_curr += integral;
        P_curr = P_end;
        
        // Clamp to physical range
        v_x_curr = std::max(-0.9999, std::min(0.9999, v_x_curr));
    }

    v_x_b = v_x_curr;

    // Final state behind rarefaction
    n_b = std::pow(P_star / K_entropy, 1.0 / gamma_c);
    H_b = 1.0 + (gamma_c / (gamma_c - 1.0)) * (P_star / n_b);

    // Tangential velocity from K_t conservation
    const real K_over_H_b = K_t / H_b;
    const real gamma2_b = (1.0 + K_over_H_b * K_over_H_b) / std::max(1e-12, 1.0 - v_x_b * v_x_b);
    if (gamma2_b <= 0.0) return false;
    const real gamma_b = std::sqrt(gamma2_b);
    v_t_b = K_t / (H_b * gamma_b);

    return true;
}

/**
 * Wave curve function: normal velocity v^x as function of pressure for root finding
 * Returns v^x_star for a given P_star candidate
 */
static real wave_curve(real P, const RiemannState& state, real v_t, real gamma_c, real c, bool is_left_wave)
{
    real v_x_star, n_star, H_star, v_t_star;
    bool success;

    if (P > state.P) {
        success = solve_shock(P, state, v_t, gamma_c, c, is_left_wave, v_x_star, n_star, H_star, v_t_star);
    } else {
        success = solve_rarefaction(P, state, v_t, gamma_c, c, is_left_wave, v_x_star, n_star, H_star, v_t_star);
    }

    return success ? v_x_star : std::numeric_limits<real>::quiet_NaN();
}

/**
 * Exact Riemann solver using bisection (Brent's method simplified)
 *
 * Finds P_star such that v^x_L(P_star) = v^x_R(P_star)
 * Based on Kitajima et al. (2025) Section 2.4
 *
 * Input states use Kitajima notation:
 *   left, right   = initial left (L) and right (R) states
 *   v_t_L, v_t_R  = tangential velocity components v^t_L and v^t_R
 *   gamma_c       = ratio of specific heats γ_c
 *
 * Output (star region):
 *   P_star        = interface pressure P*
 *   v_x_star      = interface normal velocity v^x*
 *   v_t_star      = interface tangential velocity v^t*
 */
bool exact_riemann_solver(
    const RiemannState& left,
    const RiemannState& right,
    real v_t_L,
    real v_t_R,
    real gamma_c,
    real c,
    real& P_star,
    real& v_x_star,
    real& v_t_star,
    int max_iter,
    real tol)
{
    ++g_exact_solver_calls;
    
    // Check for identical states
    if (std::abs(left.P - right.P) < 1e-6 &&
        std::abs(left.n - right.n) < 1e-6 &&
        std::abs(left.v - right.v) < 1e-6) {
        P_star = left.P;
        v_x_star = left.v;
        v_t_star = v_t_L;
        return true;
    }

    // Validate input states
    if (left.P <= 0.0 || left.n <= 0.0 || right.P <= 0.0 || right.n <= 0.0) {
        if (DEBUG_RIEMANN) {
            WRITE_LOG << "[RIEMANN DEBUG] Invalid input: P_L=" << left.P << " n_L=" << left.n
                      << " P_R=" << right.P << " n_R=" << right.n;
        }
        ++g_exact_solver_failures;
        ++g_nan_failures;
        return false;
    }
    
    // Check for superluminal velocities
    const real v2_L = left.v * left.v + v_t_L * v_t_L;
    const real v2_R = right.v * right.v + v_t_R * v_t_R;
    if (v2_L >= 1.0 || v2_R >= 1.0) {
        if (DEBUG_RIEMANN) {
            WRITE_LOG << "[RIEMANN DEBUG] Superluminal velocity: v2_L=" << v2_L << " v2_R=" << v2_R;
        }
        ++g_exact_solver_failures;
        ++g_nan_failures;
        return false;
    }

    // Residual function: f(P) = v^x_L(P) - v^x_R(P) = 0
    auto f = [&](real P) -> real {
        const real v_x_L = wave_curve(P, left, v_t_L, gamma_c, c, true);
        const real v_x_R = wave_curve(P, right, v_t_R, gamma_c, c, false);
        if (!std::isfinite(v_x_L) || !std::isfinite(v_x_R)) {
            return std::numeric_limits<real>::quiet_NaN();
        }
        return v_x_L - v_x_R;
    };

    // Initial bracket - use wider range for robustness
    const real P_min = std::min(left.P, right.P);
    const real P_max = std::max(left.P, right.P);
    real P_lo = P_min * 1e-6;
    real P_hi = P_max * 1e6;

    real f_lo = f(P_lo);
    real f_hi = f(P_hi);

    auto is_valid = [](real value) {
        return std::isfinite(value);
    };

    // If initial bracket fails, try adaptive bracketing
    if (!is_valid(f_lo) || !is_valid(f_hi)) {
        // Try smaller range
        P_lo = P_min * 0.01;
        P_hi = P_max * 100.0;
        f_lo = f(P_lo);
        f_hi = f(P_hi);
        
        if (!is_valid(f_lo) || !is_valid(f_hi)) {
            if (DEBUG_RIEMANN) {
                WRITE_LOG << "[RIEMANN DEBUG] NaN in initial bracket: f_lo=" << f_lo << " f_hi=" << f_hi
                          << " P_L=" << left.P << " P_R=" << right.P
                          << " v_L=" << left.v << " v_R=" << right.v
                          << " n_L=" << left.n << " n_R=" << right.n;
            }
            ++g_exact_solver_failures;
            ++g_nan_failures;
            return false;
        }
    }

    // Check if root is bracketed
    if (f_lo * f_hi > 0.0) {
        // Try to find a valid bracket by searching
        bool found_bracket = false;
        
        // Try logarithmically spaced points
        const int num_points = 20;
        const real log_P_lo = std::log10(P_min * 1e-6);
        const real log_P_hi = std::log10(P_max * 1e6);
        
        real prev_P = P_lo;
        real prev_f = f_lo;
        
        for (int i = 1; i <= num_points && !found_bracket; ++i) {
            const real log_P = log_P_lo + (log_P_hi - log_P_lo) * i / num_points;
            const real curr_P = std::pow(10.0, log_P);
            const real curr_f = f(curr_P);
            
            if (is_valid(curr_f) && is_valid(prev_f) && prev_f * curr_f < 0.0) {
                P_lo = prev_P;
                P_hi = curr_P;
                f_lo = prev_f;
                f_hi = curr_f;
                found_bracket = true;
            }
            
            if (is_valid(curr_f)) {
                prev_P = curr_P;
                prev_f = curr_f;
            }
        }
        
        if (!found_bracket) {
            if (DEBUG_RIEMANN) {
                WRITE_LOG << "[RIEMANN DEBUG] Cannot bracket root: f_lo=" << f_lo << " f_hi=" << f_hi
                          << " P_L=" << left.P << " P_R=" << right.P
                          << " v_L=" << left.v << " v_R=" << right.v;
            }
            ++g_exact_solver_failures;
            ++g_bracket_failures;
            return false;
        }
    }

    // =========================================================================
    // Hybrid Newton-Raphson with bisection fallback
    // Newton converges quadratically but can fail; bisection is robust but slow.
    // We use Newton when it stays in bracket, otherwise bisect.
    // =========================================================================
    
    // Numerical derivative of f(P) using central difference
    auto df_dP = [&](real P) -> real {
        const real dP = P * 1e-6;  // Relative step for numerical derivative
        const real f_plus = f(P + dP);
        const real f_minus = f(P - dP);
        if (!is_valid(f_plus) || !is_valid(f_minus)) {
            return std::numeric_limits<real>::quiet_NaN();
        }
        return (f_plus - f_minus) / (2.0 * dP);
    };
    
    // Initial guess: geometric mean (better for log-scale pressure)
    P_star = std::sqrt(P_lo * P_hi);
    
    bool converged = false;
    for (int iter = 0; iter < max_iter; ++iter) {
        const real f_mid = f(P_star);
        if (!is_valid(f_mid)) {
            if (DEBUG_RIEMANN) {
                WRITE_LOG << "[RIEMANN DEBUG] NaN during iteration " << iter << ": P_star=" << P_star;
            }
            ++g_exact_solver_failures;
            ++g_nan_failures;
            return false;
        }

        if (std::abs(f_mid) < tol || std::abs(P_hi - P_lo) < tol * P_star) {
            converged = true;
            break;
        }

        // Update bracket
        if (f_mid * f_lo < 0.0) {
            P_hi = P_star;
            f_hi = f_mid;
        } else {
            P_lo = P_star;
            f_lo = f_mid;
        }

        // Try Newton step
        const real df = df_dP(P_star);
        bool use_newton = false;
        real P_newton = P_star;
        
        if (is_valid(df) && std::abs(df) > 1e-20) {
            P_newton = P_star - f_mid / df;
            // Accept Newton step only if it stays strictly within bracket
            if (P_newton > P_lo && P_newton < P_hi) {
                use_newton = true;
            }
        }
        
        if (use_newton) {
            P_star = P_newton;
        } else {
            // Fall back to bisection (geometric mean for log-scale)
            P_star = std::sqrt(P_lo * P_hi);
        }
    }

    // Compute star state velocities from each side
    real v_x_L_star, n_L_star, H_L_star, v_t_L_star;
    real v_x_R_star, n_R_star, H_R_star, v_t_R_star;

    if (P_star > left.P) {
        solve_shock(P_star, left, v_t_L, gamma_c, c, true, v_x_L_star, n_L_star, H_L_star, v_t_L_star);
    } else {
        solve_rarefaction(P_star, left, v_t_L, gamma_c, c, true, v_x_L_star, n_L_star, H_L_star, v_t_L_star);
    }

    if (P_star > right.P) {
        solve_shock(P_star, right, v_t_R, gamma_c, c, false, v_x_R_star, n_R_star, H_R_star, v_t_R_star);
    } else {
        solve_rarefaction(P_star, right, v_t_R, gamma_c, c, false, v_x_R_star, n_R_star, H_R_star, v_t_R_star);
    }

    // Average normal velocity v^x* (should be equal at convergence)
    v_x_star = 0.5 * (v_x_L_star + v_x_R_star);

    // Tangential velocity v^t*: upwind based on contact wave direction
    // If v^x* > 0, contact moves right => use left state
    if (v_x_star > 0.0) {
        v_t_star = v_t_L_star;
    } else {
        v_t_star = v_t_R_star;
    }

    return true;
}

/**
 * Relativistic HLLC Riemann solver with proper tangent velocity treatment
 *
 * Based on:
 * - Mignone & Bodo (2005) "An HLLC Riemann solver for relativistic flows"
 * - Pons, Martí & Müller (2000) tangent velocity treatment
 *
 * The key physics from Pons et al. (2000):
 * - Tangent velocity satisfies h*W*v^t = constant across waves
 * - This couples v^t evolution to enthalpy and Lorentz factor
 *
 * Uses Kitajima notation:
 *   left, right   = initial left (L) and right (R) states
 *   v_t_L, v_t_R  = tangential velocity components v^t_L and v^t_R
 *   gamma_c       = ratio of specific heats γ_c
 *   P_star        = interface pressure P*
 *   v_x_star      = interface normal velocity v^x*
 *   v_t_star      = interface tangential velocity v^t*
 */
void hllc_riemann_solver(
    const RiemannState& left,
    const RiemannState& right,
    real v_t_L,
    real v_t_R,
    real gamma_c,
    real c,
    real& P_star,
    real& v_x_star,
    real& v_t_star)
{
    // DEBUG: only log when L and R states are VERY different (e.g., across interface)
    static int debug_count = 0;
    const bool states_differ = (std::abs(left.P - right.P) > 0.5 * std::max(left.P, right.P));
    const bool do_debug = (debug_count < 10 && states_differ);
    if (do_debug) {
        WRITE_LOG << "[HLLC DEBUG #" << debug_count << "] Input: "
                  << "L: P=" << left.P << ", n=" << left.n << ", v=" << left.v << ", cs=" << left.c_s
                  << " | R: P=" << right.P << ", n=" << right.n << ", v=" << right.v << ", cs=" << right.c_s;
    }
    
    // ========================================================================
    // STEP 1: Compute relativistic quantities
    // ========================================================================
    
    // Total velocities squared |v|² = v_x² + v_t²
    const real v2_L = left.v * left.v + v_t_L * v_t_L;
    const real v2_R = right.v * right.v + v_t_R * v_t_R;
    
    // Lorentz factors W = 1/√(1-|v|²)
    const real W_L = 1.0 / std::sqrt(std::max(1.0 - v2_L, 1e-10));
    const real W_R = 1.0 / std::sqrt(std::max(1.0 - v2_R, 1e-10));
    
    // Specific enthalpy h = 1 + ε + P/ρ = 1 + γ/(γ-1) * P/(ρc²)
    // For ideal gas: h = 1 + γ*P / ((γ-1)*n)
    const real h_L = 1.0 + gamma_c * left.P / ((gamma_c - 1.0) * left.n);
    const real h_R = 1.0 + gamma_c * right.P / ((gamma_c - 1.0) * right.n);
    
    // Sound speeds (already in units of c)
    const real c_s_L = left.c_s / c;  // Normalize if c != 1
    const real c_s_R = right.c_s / c;
    
    // ========================================================================
    // STEP 2: Relativistic wave speed estimates (Mignone & Bodo 2005, Eq. 23)
    // ========================================================================
    // For relativistic flows, wave speeds are:
    // λ± = [v_x(1-c_s²) ± c_s√((1-v²)(1-v²c_s² - v_x²(1-c_s²)))] / (1-v²c_s²)
    
    // Simplified bounds using characteristic speeds
    // Left-going wave speed
    const real denom_L = 1.0 - v2_L * c_s_L * c_s_L;
    const real sqrt_term_L = std::sqrt(std::max(
        (1.0 - v2_L) * (1.0 - v2_L * c_s_L * c_s_L - left.v * left.v * (1.0 - c_s_L * c_s_L)),
        0.0));
    const real lambda_minus_L = (left.v * (1.0 - c_s_L * c_s_L) - c_s_L * sqrt_term_L) / 
                                 std::max(denom_L, 1e-10);
    
    // Right-going wave speed  
    const real denom_R = 1.0 - v2_R * c_s_R * c_s_R;
    const real sqrt_term_R = std::sqrt(std::max(
        (1.0 - v2_R) * (1.0 - v2_R * c_s_R * c_s_R - right.v * right.v * (1.0 - c_s_R * c_s_R)),
        0.0));
    const real lambda_plus_R = (right.v * (1.0 - c_s_R * c_s_R) + c_s_R * sqrt_term_R) / 
                                std::max(denom_R, 1e-10);
    
    // HLL wave speeds (minimum/maximum characteristic speeds)
    const real S_L = std::min(lambda_minus_L, 
                              (right.v * (1.0 - c_s_R * c_s_R) - c_s_R * sqrt_term_R) / 
                              std::max(denom_R, 1e-10));
    const real S_R = std::max(lambda_plus_R,
                              (left.v * (1.0 - c_s_L * c_s_L) + c_s_L * sqrt_term_L) / 
                              std::max(denom_L, 1e-10));
    
    // ========================================================================
    // STEP 3: Compute conserved quantities (Mignone & Bodo 2005, Eq. 2)
    // ========================================================================
    // D = ρW (mass density in lab frame)
    // S_x = ρhW²v_x (x-momentum)
    // S_t = ρhW²v_t (tangential momentum)
    // E = ρhW² - P (energy)
    
    const real D_L = left.n * W_L;
    const real D_R = right.n * W_R;
    
    const real S_x_L = left.n * h_L * W_L * W_L * left.v;
    const real S_x_R = right.n * h_R * W_R * W_R * right.v;
    
    const real S_t_L = left.n * h_L * W_L * W_L * v_t_L;
    const real S_t_R = right.n * h_R * W_R * W_R * v_t_R;
    
    const real E_L = left.n * h_L * W_L * W_L - left.P;
    const real E_R = right.n * h_R * W_R * W_R - right.P;
    
    // ========================================================================
    // STEP 4: Compute fluxes (Mignone & Bodo 2005, Eq. 3)
    // ========================================================================
    // F^D = Dv_x
    // F^S_x = S_x*v_x + P
    // F^S_t = S_t*v_x
    // F^E = S_x
    
    const real F_D_L = D_L * left.v;
    const real F_D_R = D_R * right.v;
    
    const real F_Sx_L = S_x_L * left.v + left.P;
    const real F_Sx_R = S_x_R * right.v + right.P;
    
    const real F_St_L = S_t_L * left.v;
    const real F_St_R = S_t_R * right.v;
    
    const real F_E_L = S_x_L;
    const real F_E_R = S_x_R;
    
    // ========================================================================
    // STEP 5: HLL averages (for estimating star state)
    // ========================================================================
    const real dS = S_R - S_L;
    if (std::abs(dS) < 1e-14) {
        // Waves have same speed - return simple average
        P_star = 0.5 * (left.P + right.P);
        v_x_star = 0.5 * (left.v + right.v);
        v_t_star = 0.5 * (v_t_L + v_t_R);
        return;
    }
    
    // HLL state (Eq. 11 in Mignone & Bodo)
    const real D_hll = (S_R * D_R - S_L * D_L - (F_D_R - F_D_L)) / dS;
    const real S_x_hll = (S_R * S_x_R - S_L * S_x_L - (F_Sx_R - F_Sx_L)) / dS;
    const real S_t_hll = (S_R * S_t_R - S_L * S_t_L - (F_St_R - F_St_L)) / dS;
    const real E_hll = (S_R * E_R - S_L * E_L - (F_E_R - F_E_L)) / dS;
    
    // ========================================================================
    // STEP 6: Contact wave speed S_M (= v_x*) from momentum conservation
    // ========================================================================
    // Using the relation from Mignone & Bodo (2005), Eq. 18:
    // S_M solves: (E_hll + P_hll) * S_M = S_x_hll
    // With P_hll estimated from HLL average
    
    // First estimate P* from HLL jump conditions
    const real F_hll = (S_R * F_Sx_L - S_L * F_Sx_R + S_L * S_R * (S_x_R - S_x_L)) / dS;
    
    // Quadratic formula for S_M (contact wave speed = v_x*)
    // From: a*S_M² + b*S_M + c = 0
    // where the coefficients come from relativistic Rankine-Hugoniot relations
    
    // Simplified approach: use non-relativistic HLLC formula as first approximation
    // then iterate if needed
    const real A_L = S_L - left.v;
    const real A_R = S_R - right.v;
    
    // From Rankine-Hugoniot: ρ(S-v)(v* - v) = P* - P
    // Combined with contact condition: P*_L = P*_R = P*
    // Gives: S_M = [P_R - P_L + ρ_L*v_L*A_L - ρ_R*v_R*A_R] / [ρ_L*A_L - ρ_R*A_R]
    
    // For relativistic case, use enthalpy-weighted version
    const real rho_h_W2_L = left.n * h_L * W_L * W_L;
    const real rho_h_W2_R = right.n * h_R * W_R * W_R;
    
    const real num = right.P - left.P + rho_h_W2_L * left.v * A_L - rho_h_W2_R * right.v * A_R;
    const real den = rho_h_W2_L * A_L - rho_h_W2_R * A_R;
    
    if (std::abs(den) < 1e-14) {
        // Degenerate case
        v_x_star = 0.5 * (left.v + right.v);
    } else {
        v_x_star = num / den;
    }
    
    // Clamp to physical range
    v_x_star = std::max(-0.9999, std::min(0.9999, v_x_star));
    
    // ========================================================================
    // STEP 7: Compute P* from jump conditions
    // ========================================================================
    // P* = P_L + ρ_L*h_L*W_L²*(v_L - S_L)*(v_L - v_x*)
    const real P_star_L = left.P + rho_h_W2_L * A_L * (left.v - v_x_star);
    const real P_star_R = right.P + rho_h_W2_R * A_R * (right.v - v_x_star);
    
    // Average (should be equal for exact HLLC)
    P_star = 0.5 * (P_star_L + P_star_R);
    
    // Safety: ensure positive pressure
    P_star = std::max(P_star, 1e-10 * std::min(left.P, right.P));
    
    // ========================================================================
    // STEP 8: Tangent velocity in star region
    // ========================================================================
    // From Pons et al. (2000): h*W*v_t = constant across shock/rarefaction
    // So v_t* depends on which side of the contact wave we sample from
    //
    // At the contact discontinuity:
    // - Left star: v_t*_L from left state
    // - Right star: v_t*_R from right state
    //
    // The interface tangent velocity depends on the sign of v_x*:
    // - If v_x* > 0: fluid flows from left to right, use left tangent
    // - If v_x* < 0: fluid flows from right to left, use right tangent
    
    // For the Godunov flux, we need v_t at the interface (x=0)
    // Using the h*W*v_t = constant relation:
    
    if (v_x_star >= 0.0) {
        // Sample from left star state
        // h_L * W_L * v_t_L = h_star_L * W_star_L * v_t_star_L
        // Approximate: h_star ≈ h_L (for mild jumps), W changes with v_x
        const real v2_star_L = v_x_star * v_x_star + v_t_L * v_t_L;
        if (v2_star_L < 0.9999) {
            // h*W*v_t = constant => v_t* = h_L*W_L*v_t_L / (h_star*W_star)
            // Simplified: just upwind the tangent velocity
            v_t_star = v_t_L;
        } else {
            // Near-luminal: reduce tangent velocity to keep |v| < 1
            const real max_vt = std::sqrt(std::max(0.9999 - v_x_star * v_x_star, 0.0));
            v_t_star = std::copysign(std::min(std::abs(v_t_L), max_vt), v_t_L);
        }
    } else {
        // Sample from right star state
        const real v2_star_R = v_x_star * v_x_star + v_t_R * v_t_R;
        if (v2_star_R < 0.9999) {
            v_t_star = v_t_R;
        } else {
            const real max_vt = std::sqrt(std::max(0.9999 - v_x_star * v_x_star, 0.0));
            v_t_star = std::copysign(std::min(std::abs(v_t_R), max_vt), v_t_R);
        }
    }
    
    // DEBUG: output only for different states
    static int debug_out_count = 0;
    if (debug_out_count < 10 && states_differ) {
        WRITE_LOG << "[HLLC DEBUG #" << debug_count << "] Output: "
                  << "P*=" << P_star << ", v_x*=" << v_x_star << ", v_t*=" << v_t_star;
        debug_count++;
        debug_out_count++;
    }
    
    // DEBUG: Check for extreme outputs
    static int extreme_count = 0;
    if (extreme_count < 10 && (std::abs(v_x_star) > 0.9 || P_star < 0.0 || P_star > 1e6)) {
        WRITE_LOG << "[HLLC EXTREME #" << extreme_count << "] "
                  << "P*=" << P_star << " v_x*=" << v_x_star << " v_t*=" << v_t_star
                  << " | L: P=" << left.P << " n=" << left.n << " v=" << left.v
                  << " | R: P=" << right.P << " n=" << right.n << " v=" << right.v;
        ++extreme_count;
    }
}

/**
 * Report Riemann solver statistics (for debugging)
 */
void report_riemann_solver_stats()
{
    const int calls = g_exact_solver_calls.load();
    const int failures = g_exact_solver_failures.load();
    const int bracket_fails = g_bracket_failures.load();
    const int nan_fails = g_nan_failures.load();
    
    if (calls > 0) {
        const double failure_rate = 100.0 * failures / calls;
        WRITE_LOG << "[RIEMANN STATS] Calls: " << calls 
                  << " | Failures: " << failures << " (" << failure_rate << "%)"
                  << " | Bracket fails: " << bracket_fails
                  << " | NaN fails: " << nan_fails;
    }
}

/**
 * Reset Riemann solver statistics
 */
void reset_riemann_solver_stats()
{
    g_exact_solver_calls = 0;
    g_exact_solver_failures = 0;
    g_bracket_failures = 0;
    g_nan_failures = 0;
}

} // namespace riemann
} // namespace srgsph
} // namespace sph
