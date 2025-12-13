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
 * Gauss-Legendre quadrature nodes and weights for 16-point integration
 * Used for Riemann invariant integration in rarefaction waves
 */
static constexpr int GAUSS_POINTS = 16;
static const real GAUSS_NODES[GAUSS_POINTS] = {
    -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
    -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.0950125098376374,
     0.0950125098376374,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
     0.7554044083550030,  0.8656312023878318,  0.9445750230732326,  0.9894009349916499
};
static const real GAUSS_WEIGHTS[GAUSS_POINTS] = {
    0.0271524594117541, 0.0622535239386479, 0.0951585116824928, 0.1246289712555339,
    0.1495959888165767, 0.1691565193950025, 0.1826034150449236, 0.1894506104550685,
    0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
    0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541
};

/**
 * Compute tangential velocity v_t from conserved invariant A and normal velocity v_x
 * 
 * Based on Python Util.py computeVt():
 *   v_t = A * sqrt((1 - v_x²) / (h² + A²))
 * 
 * where A = h * γ * v_t is the conserved tangential momentum invariant (Pons et al. 2000)
 * 
 * @param A Conserved invariant A = h * W * v_t
 * @param v_x Normal velocity
 * @param h Specific enthalpy
 * @return Tangential velocity v_t
 */
static real compute_vt_from_invariant(real A, real v_x, real h)
{
    if (std::abs(A) < 1e-15) return 0.0;
    
    const real one_minus_vx2 = std::max(0.0, 1.0 - v_x * v_x);
    const real denom = h * h + A * A;
    if (denom < 1e-30) return 0.0;
    
    return A * std::sqrt(one_minus_vx2 / denom);
}

/**
 * Isentropic state computation from EntropyConservation
 * Based on Python EquationOfState.py EntropyConservation class
 */
struct IsentropeState {
    real K;       // Isentropic constant K = P / n^γ
    real gamma_c; // Adiabatic index
    
    IsentropeState(real P, real n, real gamma) : gamma_c(gamma) {
        K = P / std::pow(std::max(n, 1e-15), gamma);
    }
    
    real computeRho(real P) const {
        return std::pow(std::max(P, 1e-15) / K, 1.0 / gamma_c);
    }
    
    real computeEnthalpy(real P) const {
        const real rho = computeRho(P);
        return 1.0 + (gamma_c / (gamma_c - 1.0)) * P / rho;
    }
    
    real computeEnthalpy(real P, real rho) const {
        return 1.0 + (gamma_c / (gamma_c - 1.0)) * P / rho;
    }
    
    real computeSpeedOfSound(real P) const {
        const real rho = computeRho(P);
        const real h = computeEnthalpy(P, rho);
        return std::sqrt(gamma_c * P / (rho * h));
    }
    
    real computeSpeedOfSound(real P, real rho, real h) const {
        return std::sqrt(gamma_c * P / (rho * h));
    }
};

/**
 * Compute the Riemann invariant integrand for rarefaction with tangential velocity
 * 
 * Based on Python Rarefaction.py computeVxb():
 *   integrand(P) = sqrt(h² + A²(1 - c_s²)) / ((h² + A²) * ρ * c_s)
 * 
 * This is the integrand for:
 *   B = 0.5 * ln((1+v_x)/(1-v_x)) + sign * ∫[P_a to P_b] integrand(P) dP
 *   v_x_b = tanh(B)
 * 
 * @param P Pressure
 * @param A_sqr Square of conserved invariant A² = (h*W*v_t)²
 * @param isentrope Isentropic state for computing thermodynamic quantities
 * @return Integrand value
 */
static real rarefaction_integrand(real P, real A_sqr, const IsentropeState& isentrope)
{
    const real P_safe = std::max(P, 1e-15);
    const real rho = isentrope.computeRho(P_safe);
    const real h = isentrope.computeEnthalpy(P_safe, rho);
    const real cs = isentrope.computeSpeedOfSound(P_safe, rho, h);
    
    if (cs < 1e-15 || rho < 1e-15) {
        return std::numeric_limits<real>::quiet_NaN();
    }
    
    const real h_sqr = h * h;
    const real cs_sqr = cs * cs;
    
    // Integrand = sqrt(h² + A²(1 - c_s²)) / ((h² + A²) * ρ * c_s)
    const real numerator_sqr = h_sqr + A_sqr * (1.0 - cs_sqr);
    if (numerator_sqr < 0.0) {
        return std::numeric_limits<real>::quiet_NaN();
    }
    
    const real numerator = std::sqrt(numerator_sqr);
    const real denominator = (h_sqr + A_sqr) * rho * cs;
    
    if (denominator < 1e-30) {
        return std::numeric_limits<real>::quiet_NaN();
    }
    
    return numerator / denominator;
}

/**
 * Calculate velocity behind rarefaction wave using arctanh transformation
 * 
 * EXACTLY follows Python Rarefaction.py computeVxb():
 *   B = 0.5 * ln((1 + v_x_a) / (1 - v_x_a)) + sign * ∫[P_a to P_b] integrand(P) dP
 *   v_x_b = tanh(B)
 * 
 * where:
 *   integrand(P) = sqrt(h² + A²(1 - c_s²)) / ((h² + A²) * ρ * c_s)
 *   A = h_a * W_a * v_t_a (conserved tangential momentum invariant)
 *   sign = -1 for left wave, +1 for right wave
 * 
 * The tangential velocity v_t_b is computed from the invariant:
 *   v_t_b = A * sqrt((1 - v_x_b²) / (h_b² + A²))
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
    if (std::abs(v_t_a) < 1e-10) {
        return solve_rarefaction_zero_tangent(P_star, state, gamma_c, is_left_wave, v_x_b, n_b, H_b, v_t_b);
    }

    if (P_star > P_a) {
        return false; // Should be handled by shock solver
    }

    const real n_a = state.n;
    const real v_x_a = state.v;
    
    // Create isentrope from initial state (Python: EntropyConservation.fromState)
    const IsentropeState isentrope(P_a, n_a, gamma_c);
    
    // Compute conserved invariant A = h * W * v_t (Python: computeA)
    const real H_a = isentrope.computeEnthalpy(P_a);
    const real v2_total_a = v_x_a * v_x_a + v_t_a * v_t_a;
    const real v2_clamped = std::min(v2_total_a, 1.0 - 1e-12);
    const real W_a = 1.0 / std::sqrt(1.0 - v2_clamped);  // Lorentz factor
    const real A = H_a * W_a * v_t_a;  // Conserved tangential momentum invariant
    const real A_sqr = A * A;
    
    // Sign convention: -1 for left wave, +1 for right wave
    const real sign = is_left_wave ? -1.0 : 1.0;
    
    // ========================================================================
    // ARCTANH TRANSFORMATION (exactly as Python Rarefaction.py)
    // ========================================================================
    // B = 0.5 * ln((1 + v_x_a) / (1 - v_x_a)) + sign * ∫[P_a to P_b] integrand(P) dP
    // v_x_b = tanh(B)
    
    // Initial term: arctanh(v_x_a) = 0.5 * ln((1+v_x)/(1-v_x))
    const real v_x_a_clamped = std::max(-0.9999, std::min(0.9999, v_x_a));
    const real arctanh_vxa = 0.5 * std::log((1.0 + v_x_a_clamped) / (1.0 - v_x_a_clamped));
    
    // Integrate from P_a to P_star using Gauss-Legendre quadrature
    // Transform integral from [P_a, P_star] to [-1, 1]
    const real P_mid = 0.5 * (P_a + P_star);
    const real half_dP = 0.5 * (P_star - P_a);
    
    real integral = 0.0;
    for (int i = 0; i < GAUSS_POINTS; ++i) {
        const real t = GAUSS_NODES[i];
        const real P = P_mid + half_dP * t;
        
        const real f = rarefaction_integrand(P, A_sqr, isentrope);
        if (!std::isfinite(f)) {
            return false;
        }
        
        integral += GAUSS_WEIGHTS[i] * f;
    }
    integral *= half_dP;  // Jacobian of transformation
    
    // B = arctanh(v_x_a) + sign * integral
    const real B = arctanh_vxa + sign * integral;
    
    // v_x_b = tanh(B)
    v_x_b = std::tanh(B);
    
    // Clamp to physical range
    v_x_b = std::max(-0.9999, std::min(0.9999, v_x_b));
    
    // ========================================================================
    // FINAL STATE COMPUTATION
    // ========================================================================
    
    // Density and enthalpy from isentropic relations
    n_b = isentrope.computeRho(P_star);
    H_b = isentrope.computeEnthalpy(P_star, n_b);
    
    // Tangential velocity from conserved invariant (Python: computeVt)
    // v_t_b = A * sqrt((1 - v_x_b²) / (h_b² + A²))
    v_t_b = compute_vt_from_invariant(A, v_x_b, H_b);
    
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
