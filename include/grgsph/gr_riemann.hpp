#pragma once

#include "defines.hpp"
#include "grgsph/gr_metric.hpp"
#include "srgsph/sr_exact_riemann.hpp"

namespace sph {
namespace grgsph {

/**
 * Result of GR Riemann solver
 *
 * All velocities are in the Eulerian observer frame (physical velocities).
 */
struct GRRiemannResult {
    real P_star;      // Interface pressure
    real v_x_star;    // Interface normal velocity (Eulerian frame)
    real v_t_star;    // Interface tangent velocity (Eulerian frame)
    real n_star;      // Interface density (estimated)

    real S_L;         // Left wave speed (for timestep)
    real S_R;         // Right wave speed (for timestep)

    bool converged;   // Did the solver converge?
};

/**
 * State for GR Riemann problem
 *
 * Contains both primitive variables and metric information.
 */
struct GRRiemannState {
    // Primitive variables (rest frame)
    real n;           // Rest-frame baryon density
    real P;           // Pressure
    real c_s;         // Sound speed

    // Velocity components in Eulerian frame
    real V_x;         // Normal velocity (line of sight)
    real V_t;         // Tangent velocity magnitude

    // Derived quantities
    real h;           // Specific enthalpy h = 1 + ε + P/ρ
    real W;           // Lorentz factor

    // Local metric (optional, for full GR treatment)
    const Metric31* metric;

    GRRiemannState() : n(0), P(0), c_s(0), V_x(0), V_t(0), h(1), W(1), metric(nullptr) {}
};

/**
 * GR-aware Riemann Solver
 *
 * Solves the Riemann problem in curved spacetime by:
 * 1. Transforming to local Eulerian observer frame
 * 2. Solving SR Riemann problem (reusing existing sr_exact_riemann)
 * 3. Returning result in Eulerian frame
 *
 * Key insight: In the local Eulerian observer frame, the metric is
 * effectively Minkowski, so the SR Riemann solver is applicable.
 * The GR effects come from:
 * - Velocity transformation (coordinate → Eulerian)
 * - Gravitational source terms (handled separately in force calculation)
 */
class GRRiemannSolver {
public:
    /**
     * Solve Riemann problem in curved spacetime
     *
     * @param left Left state (primitives in Eulerian frame)
     * @param right Right state (primitives in Eulerian frame)
     * @param gamma_ad Adiabatic index
     * @param use_exact Use exact solver (true) or HLLC (false)
     * @return GRRiemannResult with interface state
     */
    static GRRiemannResult solve(
        const GRRiemannState& left,
        const GRRiemannState& right,
        real gamma_ad,
        bool use_exact = true
    );

    /**
     * Solve with coordinate velocities and metric
     *
     * Automatically transforms velocities to Eulerian frame before solving.
     *
     * @param left Left state (coordinate frame velocities)
     * @param right Right state (coordinate frame velocities)
     * @param v_coord_left Coordinate velocity vector (left)
     * @param v_coord_right Coordinate velocity vector (right)
     * @param n_ij Unit vector from left to right particle
     * @param metric_left Metric at left position
     * @param metric_right Metric at right position
     * @param gamma_ad Adiabatic index
     * @return GRRiemannResult with interface state
     */
    static GRRiemannResult solve_with_metric(
        real n_left, real P_left, real c_s_left,
        real n_right, real P_right, real c_s_right,
        const real v_coord_left[3],
        const real v_coord_right[3],
        const real n_ij[3],
        const Metric31& metric_left,
        const Metric31& metric_right,
        real gamma_ad,
        bool use_exact = true
    );

private:
    /**
     * HLLC approximate solver for GR
     *
     * Based on Mignone & Bodo (2005) with K-invariant for tangent velocity.
     */
    static GRRiemannResult solve_hllc(
        const GRRiemannState& left,
        const GRRiemannState& right,
        real gamma_ad
    );

    /**
     * Estimate wave speeds for HLLC and timestep
     */
    static void estimate_wave_speeds(
        const GRRiemannState& left,
        const GRRiemannState& right,
        real& S_L, real& S_R
    );
};

} // namespace grgsph
} // namespace sph
