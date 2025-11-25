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
 * Calculate velocity behind rarefaction wave (general case with tangential velocity)
 * Uses RK4 integration along characteristic (Pons et al. Eq. 14-17)
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

    // Use analytical solution for zero tangential velocity
    if (v_t_a < 1e-6) {
        return solve_rarefaction_zero_tangent(P_star, state, gamma_c, is_left_wave, v_x_b, n_b, H_b, v_t_b);
    }

    if (P_star > P_a) {
        return false; // Should be handled by shock solver
    }

    const real n_a = state.n;
    const real v_x_a = state.v;  // Normal velocity v^x
    const real K_entropy = P_a / std::pow(n_a, gamma_c);
    const real H_a = 1.0 + (gamma_c / (gamma_c - 1.0)) * (P_a / n_a);

    // Total velocity and Lorentz factor
    const real v2_total_a = v_x_a * v_x_a + v_t_a * v_t_a;
    const real v2_clamped = std::min(v2_total_a, 1.0 - 1e-12);
    const real gamma_a = 1.0 / std::sqrt(1.0 - v2_clamped);

    // Conserved tangential momentum: K_t = H * γ * v^t (Pons et al.)
    const real K_t = H_a * gamma_a * v_t_a;

    const real sign = is_left_wave ? -1.0 : 1.0;       // dv^x/dP sign
    const real xi_sign = is_left_wave ? -1.0 : 1.0;    // Characteristic speed sign

    // Check if pressure change is significant
    const real delta_P = std::abs(P_a - P_star);
    if (delta_P < 1e-12) {
        v_x_b = v_x_a;
        n_b = n_a;
        H_b = H_a;
        v_t_b = v_t_a;
        return true;
    }

    // Adaptive step count based on pressure span
    const int raw_steps = static_cast<int>(delta_P / std::max(P_a * 1e-3, 1e-6)) + 64;
    const int steps = std::max(64, std::min(2048, raw_steps));

    // RK4 derivative: dv^x/dP along rarefaction (Pons et al. Eq. 14-17)
    auto dv_x_dP = [&](real P, real v_x) -> real {
        const real P_safe = std::max(P, 1e-12);
        const real n = std::pow(P_safe / K_entropy, 1.0 / gamma_c);
        const real H = 1.0 + (gamma_c / (gamma_c - 1.0)) * (P_safe / n);

        // Lorentz factor including tangential velocity via K_t conservation
        const real K_over_H = K_t / H;
        const real denom_v = std::max(1e-12, 1.0 - v_x * v_x);
        const real gamma2 = (1.0 + K_over_H * K_over_H) / denom_v;
        if (gamma2 <= 0.0) return std::numeric_limits<real>::quiet_NaN();
        const real gamma = std::sqrt(gamma2);
        const real v_t = K_t / (H * gamma);

        // Sound speed c_s
        const real c_s2 = std::max(0.0, gamma_c * P_safe / (n * H));
        const real c_s = std::sqrt(c_s2);

        // Characteristic speed ξ (Pons et al. Eq. 14)
        const real v2 = v_x * v_x + v_t * v_t;
        const real term1 = v_x * (1.0 - c_s2);
        real term2_inner = (1.0 - v2) * (1.0 - v2 * c_s2 - v_x * v_x * (1.0 - c_s2));
        term2_inner = std::max(term2_inner, 0.0);
        const real term2 = c_s * std::sqrt(term2_inner);
        const real denom_xi = std::max(1e-12, 1.0 - v2 * c_s2);

        const real xi = (term1 + xi_sign * term2) / denom_xi;

        // Correction factor g for tangential velocity (Pons et al.)
        const real one_minus_xi_v = std::max(1e-8, 1.0 - xi * v_x);
        const real g = (v_t * v_t * (xi * xi - 1.0)) / (one_minus_xi_v * one_minus_xi_v);

        const real denom_deriv = n * H * gamma2 * c_s * std::sqrt(1.0 + g);
        if (denom_deriv <= 0.0 || !std::isfinite(denom_deriv)) {
            return std::numeric_limits<real>::quiet_NaN();
        }
        return sign / denom_deriv;
    };

    // RK4 integration from P_a to P_star
    real P_curr = P_a;
    real v_x_curr = v_x_a;

    for (int step = 0; step < steps; ++step) {
        const real P_next = P_a + (P_star - P_a) * static_cast<real>(step + 1) / steps;
        const real dP = P_next - P_curr;

        const real k1 = dv_x_dP(P_curr, v_x_curr);
        if (!std::isfinite(k1)) return false;
        const real k2 = dv_x_dP(P_curr + 0.5 * dP, v_x_curr + 0.5 * dP * k1);
        if (!std::isfinite(k2)) return false;
        const real k3 = dv_x_dP(P_curr + 0.5 * dP, v_x_curr + 0.5 * dP * k2);
        if (!std::isfinite(k3)) return false;
        const real k4 = dv_x_dP(P_curr + dP, v_x_curr + dP * k3);
        if (!std::isfinite(k4)) return false;

        v_x_curr += (dP / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        P_curr = P_next;
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

    // Bisection iteration
    P_star = 0.5 * (P_lo + P_hi);

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
            break;
        }

        if (f_mid * f_lo < 0.0) {
            P_hi = P_star;
            f_hi = f_mid;
        } else {
            P_lo = P_star;
            f_lo = f_mid;
        }

        P_star = 0.5 * (P_lo + P_hi);
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
 * HLLC Riemann solver (approximate but robust fallback)
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
    // Wave speed estimates
    const real c_s_L = left.c_s;
    const real c_s_R = right.c_s;

    // Left and right wave speeds S_L and S_R
    const real S_L = std::min(left.v - c_s_L, right.v - c_s_R);
    const real S_R = std::max(left.v + c_s_L, right.v + c_s_R);

    // Contact wave speed (star region normal velocity v^x*)
    const real num = right.P - left.P + left.n * left.v * (S_L - left.v) -
                     right.n * right.v * (S_R - right.v);
    const real den = left.n * (S_L - left.v) - right.n * (S_R - right.v);

    v_x_star = num / den;

    // Clamp velocity to physical range (|v| < c = 1)
    if (v_x_star >= 1.0) v_x_star = 0.99;
    if (v_x_star <= -1.0) v_x_star = -0.99;

    // Star region pressure P* (HLLC estimate)
    P_star = left.P + left.n * (left.v - S_L) * (left.v - v_x_star);

    // Safety bounds
    if (P_star < 0.0) P_star = 0.5 * (left.P + right.P);

    // Tangential velocity v^t* (upwinding based on contact wave direction)
    if (v_x_star > 0.0) {
        v_t_star = v_t_L;
    } else {
        v_t_star = v_t_R;
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
