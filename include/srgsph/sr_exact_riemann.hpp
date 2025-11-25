#pragma once

#include "defines.hpp"

namespace sph {
namespace srgsph {
namespace riemann {

/**
 * State for Riemann problem (rest-frame primitive variables)
 *
 * Variable naming follows Kitajima et al. (2025) notation:
 *   v    = normal velocity v^x (rest-frame)
 *   n    = baryon number density (rest-frame)
 *   P    = pressure
 *   c_s  = sound speed
 */
struct RiemannState {
    real v;      // Normal velocity v^x (rest-frame)
    real n;      // Baryon number density n (rest-frame)
    real P;      // Pressure P
    real c_s;    // Sound speed c_s
};

/**
 * Wave type in Riemann solution
 */
enum WaveType { RAREFACTION, SHOCK };

/**
 * Exact Riemann solver for special relativistic hydrodynamics
 * Based on Kitajima et al. (2025) and Pons et al. (2000)
 *
 * Solves the Riemann problem exactly using iterative bisection
 * to find the pressure P* in the star region, then computes velocity v^x*.
 *
 * Variable naming follows Kitajima notation (SRGSPH paper Section 2.1):
 *   N       = lab-frame baryon number density (N = γn)
 *   n       = rest-frame baryon number density
 *   γ       = Lorentz factor
 *   H       = enthalpy per baryon
 *   P       = pressure
 *   v^x     = normal velocity component
 *   v^t     = tangential velocity component
 *   c_s     = sound speed
 *   γ_c     = ratio of specific heats
 *
 * @param left Left state L (rest-frame primitives: v^x, n, P, c_s)
 * @param right Right state R (rest-frame primitives: v^x, n, P, c_s)
 * @param v_t_L Tangential velocity v^t_L (left state)
 * @param v_t_R Tangential velocity v^t_R (right state)
 * @param gamma_c Ratio of specific heats γ_c (e.g., 5/3 for ideal gas)
 * @param c Speed of light in code units (typically 1.0)
 * @param P_star [output] Pressure at interface P* (star region)
 * @param v_x_star [output] Normal velocity at interface v^x* (star region)
 * @param v_t_star [output] Tangential velocity at interface v^t* (star region)
 * @param max_iter Maximum iterations (default: 100)
 * @param tol Convergence tolerance (default: 1e-10)
 * @return true if converged, false if failed (caller should use HLLC fallback)
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
    int max_iter = 100,
    real tol = 1e-10
);

/**
 * HLLC Riemann solver (approximate but robust fallback)
 *
 * Uses same Kitajima notation as exact_riemann_solver.
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
    real& v_t_star
);

/**
 * Report Riemann solver statistics (for debugging)
 */
void report_riemann_solver_stats();

/**
 * Reset Riemann solver statistics
 */
void reset_riemann_solver_stats();

} // namespace riemann
} // namespace srgsph
} // namespace sph
