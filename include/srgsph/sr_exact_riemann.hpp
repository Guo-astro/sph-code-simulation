#pragma once

#include "defines.hpp"

namespace sph {
namespace srgsph {
namespace riemann {

/**
 * State for Riemann problem (rest-frame primitive variables)
 */
struct RiemannState {
    real v;      // velocity (rest-frame)
    real n;      // baryon number density (rest-frame)
    real P;      // pressure
    real c_s;    // sound speed
};

/**
 * Wave type in Riemann solution
 */
enum WaveType { RAREFACTION, SHOCK };

/**
 * Exact Riemann solver for special relativistic hydrodynamics
 * Based on Pons et al. (2000) J. Fluid Mech. 422, 125-139
 * 
 * Solves the Riemann problem exactly using iterative Newton-Raphson method
 * to find the pressure P* in the star region, then computes velocity v*.
 * 
 * @param left Left state (rest-frame primitives: v, n, P, c_s)
 * @param right Right state (rest-frame primitives: v, n, P, c_s)
 * @param gamma_c Ratio of specific heats (e.g., 5/3 for ideal gas)
 * @param c_speed Speed of light in code units (typically 1.0)
 * @param P_star [output] Pressure at interface (star region)
 * @param v_star [output] Velocity at interface (star region)
 * @param max_iter Maximum Newton-Raphson iterations (default: 100)
 * @param tol Convergence tolerance (default: 1e-10)
 * @return true if converged, false if failed (caller should use HLLC fallback)
 */
bool exact_riemann_solver(
    const RiemannState& left,
    const RiemannState& right,
    real gamma_c,
    real c_speed,
    real& P_star,
    real& v_star,
    int max_iter = 100,
    real tol = 1e-10
);

/**
 * HLLC Riemann solver (fallback for non-convergence or weak interactions)
 * Same as current implementation - approximate but fast
 */
void hllc_riemann_solver(
    const RiemannState& left,
    const RiemannState& right,
    real gamma_c,
    real c_speed,
    real& P_star,
    real& v_star
);

} // namespace riemann
} // namespace srgsph
} // namespace sph
