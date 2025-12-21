/**
 * GR-GSPH Riemann Solver Implementation
 *
 * Wraps the SR exact Riemann solver for use in curved spacetime.
 * The key insight is that in the local Eulerian observer frame,
 * the metric is effectively Minkowski, so SR physics applies.
 */

#include "grgsph/gr_riemann.hpp"
#include "srgsph/sr_exact_riemann.hpp"
#include <cmath>
#include <algorithm>

namespace sph {
namespace grgsph {

/**
 * Estimate wave speeds for HLLC and timestep constraint
 *
 * Uses relativistic characteristic velocities:
 *   λ± = (v ± c_s) / (1 ± v·c_s)
 */
void GRRiemannSolver::estimate_wave_speeds(
    const GRRiemannState& left,
    const GRRiemannState& right,
    real& S_L, real& S_R)
{
    // Left-going characteristic from left state
    const real denom_Lm = 1.0 - left.V_x * left.c_s;
    const real lambda_L_minus = (std::abs(denom_Lm) > 1e-10) ?
        (left.V_x - left.c_s) / denom_Lm : -1.0;

    // Right-going characteristic from left state
    const real denom_Lp = 1.0 + left.V_x * left.c_s;
    const real lambda_L_plus = (std::abs(denom_Lp) > 1e-10) ?
        (left.V_x + left.c_s) / denom_Lp : 1.0;

    // Left-going characteristic from right state
    const real denom_Rm = 1.0 - right.V_x * right.c_s;
    const real lambda_R_minus = (std::abs(denom_Rm) > 1e-10) ?
        (right.V_x - right.c_s) / denom_Rm : -1.0;

    // Right-going characteristic from right state
    const real denom_Rp = 1.0 + right.V_x * right.c_s;
    const real lambda_R_plus = (std::abs(denom_Rp) > 1e-10) ?
        (right.V_x + right.c_s) / denom_Rp : 1.0;

    // Davis-type estimates
    const real v_avg = 0.5 * (left.V_x + right.V_x);
    const real c_max = std::max(left.c_s, right.c_s);

    S_L = std::min({lambda_L_minus, lambda_R_minus, v_avg - c_max});
    S_R = std::max({lambda_L_plus, lambda_R_plus, v_avg + c_max});

    // Clamp to subluminal
    S_L = std::max(S_L, -0.9999);
    S_R = std::min(S_R, 0.9999);
}

/**
 * HLLC approximate solver
 *
 * Based on Mignone & Bodo (2005) with K-invariant for tangent velocity.
 */
GRRiemannResult GRRiemannSolver::solve_hllc(
    const GRRiemannState& left,
    const GRRiemannState& right,
    real gamma_ad)
{
    GRRiemannResult result;
    result.converged = true;

    // Compute wave speeds
    estimate_wave_speeds(left, right, result.S_L, result.S_R);

    // Acoustic impedances
    const real Z_L = left.n * left.h * left.W * left.W * left.c_s;
    const real Z_R = right.n * right.h * right.W * right.W * right.c_s;
    const real Z_sum = Z_L + Z_R;

    if (Z_sum < 1e-20) {
        // Fallback to simple average
        result.P_star = 0.5 * (left.P + right.P);
        result.v_x_star = 0.5 * (left.V_x + right.V_x);
        result.v_t_star = 0.5 * (left.V_t + right.V_t);
        result.n_star = 0.5 * (left.n + right.n);
        return result;
    }

    // Acoustic approximation for contact speed and pressure
    result.v_x_star = (Z_L * left.V_x + Z_R * right.V_x + left.P - right.P) / Z_sum;
    result.P_star = (Z_L * right.P + Z_R * left.P + Z_L * Z_R * (left.V_x - right.V_x)) / Z_sum;

    // Clamp to physical range
    result.P_star = std::max(result.P_star, 1e-15);
    result.v_x_star = std::max(-0.9999, std::min(0.9999, result.v_x_star));

    // K-invariant for tangent velocity
    const real K_L = left.h * left.W * left.V_t;
    const real K_R = right.h * right.W * right.V_t;

    // Upwind K based on contact wave direction
    const real K_star = (result.v_x_star > 0) ? K_L : K_R;

    // Estimate star region enthalpy
    const real n_upwind = (result.v_x_star > 0) ? left.n : right.n;
    const real P_upwind = (result.v_x_star > 0) ? left.P : right.P;
    result.n_star = n_upwind * std::pow(result.P_star / std::max(P_upwind, 1e-15), 1.0 / gamma_ad);
    const real h_star = 1.0 + gamma_ad / (gamma_ad - 1.0) * result.P_star / std::max(result.n_star, 1e-15);

    // Solve for v_t*: K* = h* W* v_t*
    const real v_x_star2 = result.v_x_star * result.v_x_star;
    if (std::abs(K_star) > 1e-15 && v_x_star2 < 0.9999) {
        const real K2 = K_star * K_star;
        const real h2 = h_star * h_star;
        const real v_t_star2 = K2 * (1.0 - v_x_star2) / (K2 + h2);

        // Ensure subluminal
        const real v_t_max2 = std::max(0.0, 0.9999 - v_x_star2);
        result.v_t_star = std::copysign(std::sqrt(std::min(v_t_star2, v_t_max2)), K_star);
    } else {
        result.v_t_star = (result.v_x_star > 0) ? left.V_t : right.V_t;
    }

    return result;
}

/**
 * Main solve function
 *
 * Uses exact SR Riemann solver, falling back to HLLC if it fails.
 */
GRRiemannResult GRRiemannSolver::solve(
    const GRRiemannState& left,
    const GRRiemannState& right,
    real gamma_ad,
    bool use_exact)
{
    GRRiemannResult result;

    // Estimate wave speeds (needed for timestep regardless of solver)
    estimate_wave_speeds(left, right, result.S_L, result.S_R);

    if (!use_exact) {
        return solve_hllc(left, right, gamma_ad);
    }

    // Build RiemannState for SR solver
    srgsph::riemann::RiemannState sr_left{left.V_x, left.n, left.P, left.c_s};
    srgsph::riemann::RiemannState sr_right{right.V_x, right.n, right.P, right.c_s};

    // Call exact SR Riemann solver
    const real c = 1.0;  // Speed of light in code units
    result.converged = srgsph::riemann::exact_riemann_solver(
        sr_left, sr_right,
        left.V_t, right.V_t,
        gamma_ad, c,
        result.P_star, result.v_x_star, result.v_t_star,
        100, 1e-10
    );

    if (!result.converged) {
        // Fall back to HLLC
        return solve_hllc(left, right, gamma_ad);
    }

    // Estimate density in star region
    const real n_upwind = (result.v_x_star > 0) ? left.n : right.n;
    const real P_upwind = (result.v_x_star > 0) ? left.P : right.P;
    result.n_star = n_upwind * std::pow(result.P_star / std::max(P_upwind, 1e-15), 1.0 / gamma_ad);

    return result;
}

/**
 * Solve with metric transformation
 *
 * Transforms coordinate velocities to Eulerian frame, then solves.
 */
GRRiemannResult GRRiemannSolver::solve_with_metric(
    real n_left, real P_left, real c_s_left,
    real n_right, real P_right, real c_s_right,
    const real v_coord_left[3],
    const real v_coord_right[3],
    const real n_ij[3],
    const Metric31& metric_left,
    const Metric31& metric_right,
    real gamma_ad,
    bool use_exact)
{
    // Transform coordinate velocities to Eulerian frame
    real V_left[3], V_right[3];
    metric_left.coord_to_eulerian(v_coord_left, V_left);
    metric_right.coord_to_eulerian(v_coord_right, V_right);

    // Project onto line of sight (normal direction)
    real V_x_left = 0.0, V_x_right = 0.0;
    for (int i = 0; i < 3; ++i) {
        V_x_left += n_ij[i] * V_left[i];
        V_x_right += n_ij[i] * V_right[i];
    }

    // Compute tangent velocity magnitude
    // V_t² = |V|² - V_x²
    real V2_left = 0.0, V2_right = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            V2_left += metric_left.gamma_ij[i][j] * V_left[i] * V_left[j];
            V2_right += metric_right.gamma_ij[i][j] * V_right[i] * V_right[j];
        }
    }

    real V_t_left = std::sqrt(std::max(0.0, V2_left - V_x_left * V_x_left));
    real V_t_right = std::sqrt(std::max(0.0, V2_right - V_x_right * V_x_right));

    // Compute Lorentz factors
    const real W_left = 1.0 / std::sqrt(std::max(1e-10, 1.0 - V2_left));
    const real W_right = 1.0 / std::sqrt(std::max(1e-10, 1.0 - V2_right));

    // Compute specific enthalpies
    const real h_left = 1.0 + gamma_ad / (gamma_ad - 1.0) * P_left / std::max(n_left, 1e-15);
    const real h_right = 1.0 + gamma_ad / (gamma_ad - 1.0) * P_right / std::max(n_right, 1e-15);

    // Build GRRiemannState structures
    GRRiemannState left_state, right_state;

    left_state.n = n_left;
    left_state.P = P_left;
    left_state.c_s = c_s_left;
    left_state.V_x = V_x_left;
    left_state.V_t = V_t_left;
    left_state.h = h_left;
    left_state.W = W_left;
    left_state.metric = &metric_left;

    right_state.n = n_right;
    right_state.P = P_right;
    right_state.c_s = c_s_right;
    right_state.V_x = V_x_right;
    right_state.V_t = V_t_right;
    right_state.h = h_right;
    right_state.W = W_right;
    right_state.metric = &metric_right;

    return solve(left_state, right_state, gamma_ad, use_exact);
}

} // namespace grgsph
} // namespace sph
