#pragma once

#include "defines.hpp"
#include "grgsph/gr_metric.hpp"
#include <cmath>
#include <algorithm>

namespace sph {
namespace grgsph {

/**
 * GR Conserved Variables
 *
 * Following GRSPH (Liptai & Price 2019):
 *   ρ* = √γ · n · Γ          (conserved baryon density)
 *   p_i = Γ w γ_{ij} V^j     (covariant momentum per baryon)
 *   e = U^0 w - P√(-g)/ρ*    (total specific energy)
 *
 * Or using entropy formulation:
 *   K = P / n^γ              (entropy variable)
 */
struct GRConserved {
    real rho_star;     // Conserved density ρ* = √γ n Γ
    real p[3];         // Covariant momentum p_i
    real e;            // Total specific energy (or K for entropy)

    // Whether using entropy (K) or energy (e) formulation
    bool use_entropy;

    GRConserved() : rho_star(0), e(0), use_entropy(false) {
        p[0] = p[1] = p[2] = 0.0;
    }
};

/**
 * GR Primitive Variables
 *
 * Physical quantities in rest frame or Eulerian observer frame.
 */
struct GRPrimitive {
    real n;            // Rest-frame baryon density
    real P;            // Pressure
    real V[3];         // Eulerian velocity V^i
    real Gamma;        // Lorentz factor Γ
    real h;            // Specific enthalpy
    real c_s;          // Sound speed

    GRPrimitive() : n(0), P(0), Gamma(1), h(1), c_s(0) {
        V[0] = V[1] = V[2] = 0.0;
    }
};

/**
 * Conservative to Primitive Variable Recovery for GR
 *
 * The cons2prim problem is more involved in GR due to the implicit
 * relationship between conserved and primitive variables.
 *
 * Algorithm (from GRSPH Appendix B):
 * 1. Compute |p|² = γ^{ij} p_i p_j
 * 2. Initial guess: w₀ from previous timestep (or simple estimate)
 * 3. Newton iteration on f(w) = w(ρ,P) - w = 0
 *    where ρ and P are functions of w through the conserved variables
 * 4. Recover primitives from converged w
 */
class GRCons2Prim {
public:
    /**
     * Recover primitive variables from conserved
     *
     * @param U Conserved variables
     * @param metric Local metric
     * @param gamma_ad Adiabatic index
     * @param w_prev Previous enthalpy (initial guess)
     * @param prim [output] Recovered primitives
     * @return true if converged, false if failed
     */
    static bool solve(
        const GRConserved& U,
        const Metric31& metric,
        real gamma_ad,
        real w_prev,
        GRPrimitive& prim
    );

    /**
     * Compute conserved from primitive (forward direction)
     *
     * This is straightforward and used for initialization.
     */
    static void prim2cons(
        const GRPrimitive& prim,
        const Metric31& metric,
        real gamma_ad,
        GRConserved& U
    );

    /**
     * Compute enthalpy from EOS
     * w = 1 + ε + P/ρ = 1 + γ/(γ-1) P/ρ for ideal gas
     */
    static real compute_enthalpy(real n, real P, real gamma_ad) {
        return 1.0 + gamma_ad / (gamma_ad - 1.0) * P / std::max(n, 1e-15);
    }

    /**
     * Compute sound speed
     * c_s = √(γ P / (ρ h))
     */
    static real compute_sound_speed(real n, real P, real h, real gamma_ad) {
        return std::sqrt(gamma_ad * P / (std::max(n, 1e-15) * h));
    }

private:
    // Newton-Raphson tolerance and max iterations
    static constexpr real TOL = 1e-12;
    static constexpr int MAX_ITER = 50;
};

/**
 * GR Particle for GR-GSPH
 *
 * Extended particle structure with GR-specific fields.
 */
struct GRParticle {
    // Position and velocity
    vec_t pos;         // Position x^i
    vec_t vel;         // Coordinate velocity v^i = dx^i/dt
    real vel_t;        // Tangent velocity (for 1D tests)

    // Conserved variables
    GRConserved U;

    // Primitive variables
    GRPrimitive prim;

    // SPH quantities
    real mass;         // Particle mass m
    real sml;          // Smoothing length h
    real nu;           // Volume element ν = m/ρ*
    real N;            // Lab-frame density N = n Γ
    real gradh;        // Grad-h correction Ω

    // Derivatives (accelerations)
    vec_t dS;          // Momentum derivative dp_i/dt
    real de;           // Energy derivative de/dt

    // Previous timestep enthalpy (for cons2prim)
    real w_prev;

    // Metric at particle position
    Metric31 metric;

    // Flags
    bool is_ghost;

    GRParticle() : vel_t(0), mass(0), sml(0), nu(0), N(0), gradh(1),
                   de(0), w_prev(1), is_ghost(false) {}
};

} // namespace grgsph
} // namespace sph
