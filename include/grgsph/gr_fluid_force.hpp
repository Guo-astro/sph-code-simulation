#pragma once

#include "defines.hpp"
#include "fluid_force.hpp"
#include "parameters.hpp"
#include "particle.hpp"
#include "grgsph/gr_metric.hpp"
#include "grgsph/gr_riemann.hpp"
#include "grgsph/gr_cons2prim.hpp"
#include <memory>

namespace sph {
namespace grgsph {

/**
 * GR-GSPH Fluid Force
 *
 * Computes fluid forces using Godunov approach with exact Riemann solver.
 * Based on Liptai & Price (2019) and Kitajima et al. (2025) GRGSPH formulation.
 *
 * === GRGSPH Equations (Formulation §8) ===
 *
 * 1. Momentum equation (§8.1):
 *    ⟨ν_i Ṡ_k⟩ = -Σ_j P*_ij √(-g)_ij V²_ij [∇_k W_ij(√2 h_i) - ∇_k W_ij(√2 h_j)] + ν_i f_k
 *
 *    where:
 *    - P*_ij = interface pressure from Riemann solver
 *    - √(-g)_ij = (1/2)(√(-g_i) + √(-g_j))  (metric determinant average)
 *    - V²_ij = (1/2)(V²_ij(h_i) + V²_ij(h_j))  (volume interpolation)
 *    - f_k = gravitational source term
 *
 * 2. Energy equation (§8.2):
 *    ⟨ν_i ė⟩ = -Σ_j P*_ij v*^k_ij √(-g)_ij V²_ij [∇_k W_ij(√2 h_i) - ∇_k W_ij(√2 h_j)] + ν_i Λ
 *
 *    where:
 *    - v*_ij = interface velocity from Riemann solver
 *    - Λ = gravitational work (zero for stationary metrics)
 *
 * 3. Gravitational source terms (§5.3):
 *    f_i = (√-g / 2N*) T^μν ∂g_μν/∂x^i  (momentum source)
 *    Λ = -(√-g / 2N*) T^μν ∂g_μν/∂t     (energy source, zero for stationary metrics)
 *
 * Key differences from GRSPH:
 * - No artificial viscosity q
 * - Interface pressure P* from Riemann solver
 * - Proper entropy generation at shocks only
 * - √(-g) averaging in force calculation
 * - Kernel evaluated at √2 h (from Gaussian convolution)
 */
class GRFluidForce : public sph::FluidForce {
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;

    // Get/set metric
    void set_metric(std::unique_ptr<MetricBase> metric) { m_metric = std::move(metric); }
    const MetricBase* get_metric() const { return m_metric.get(); }

private:
    // Physics parameters
    real m_gamma;              // Adiabatic index
    real m_c_speed;            // Speed of light

    // Solver options
    RiemannSolverType m_riemann_solver;
    bool m_use_muscl;          // 2nd order reconstruction

    // Shock/CD detection thresholds
    real m_c_shock;
    real m_c_cd;

    // Geodesic test mode: skip pairwise SPH forces, only use metric source
    bool m_geodesic_mode;
    real m_bh_mass;  // Black hole mass for geodesic calculations

    // Metric for GR effects
    std::unique_ptr<MetricBase> m_metric;

    /**
     * Solve interface state using Riemann solver
     *
     * @param left_state [v^x, n, P, c_s, v^t] left state
     * @param right_state [v^x, n, P, c_s, v^t] right state
     * @param P_star [output] Interface pressure
     * @param v_x_star [output] Interface normal velocity
     * @param v_t_star [output] Interface tangent velocity
     */
    void solve_interface_state(
        const real left_state[5],
        const real right_state[5],
        real& P_star,
        real& v_x_star,
        real& v_t_star) const;

    /**
     * Compute gravitational momentum source term (Formulation §5.3)
     *
     * From the formulation, the gravitational force is:
     *   f_i = (√-g / 2N*) T^μν ∂g_μν/∂x^i
     *
     * Implemented using the equivalent Valencia formulation (L&P 2019 Eq. 15):
     *   S_i^{source} = α√γ Γ^α_{iβ} T_α^β
     *
     * where:
     *   Γ^α_{iβ} = Christoffel symbol with raised first index
     *   T_α^β = g_{αμ} T^{μβ} (mixed stress-energy tensor)
     *
     * For static metrics (zero shift):
     *   - Γ^0_{i0} = (1/α) ∂α/∂x^i (gravitational redshift)
     *   - Γ^j_{ik} = spatial Christoffel symbols (velocity-dependent terms)
     *
     * @param prim Primitive variables (n, P, V^i, Γ, h)
     * @param metric Local 3+1 metric (α, β^i, γ_ij)
     * @param derivs Metric derivatives (∂g_μν/∂x^i)
     * @param N_star Conserved baryon number density N* = √γ n Γ
     * @return Gravitational source force vector (per-baryon acceleration)
     */
    vec_t compute_metric_source(
        const GRPrimitive& prim,
        const Metric31& metric,
        const MetricDerivatives& derivs,
        real N_star) const;

    /**
     * Normalize SR derivatives (same as SR-GSPH)
     * Scale derivatives by baryon number so dS and de are per-baryon quantities.
     */
    static void normalize_derivatives(SPHParticle& p);

    /**
     * Compute gravitational work (energy source) for dynamical spacetimes
     *
     * From GRGSPH Formulation §5.3:
     *   Λ = -(√-g / 2N*) T^μν ∂g_μν/∂t
     *
     * For stationary metrics (Schwarzschild, Kerr, Minkowski), ∂g_μν/∂t = 0, so Λ = 0.
     * For dynamical spacetimes (cosmological, GW backgrounds, etc.), this term
     * represents the work done by the changing gravitational field on the fluid.
     *
     * The full expansion uses the stress-energy tensor:
     *   T^μν = ρh W² u^μ u^ν + P g^μν
     *
     * where u^μ = dx^μ/dτ is the 4-velocity.
     *
     * @param prim Primitive variables (n, P, V^i, Γ, h)
     * @param metric Local 3+1 metric (α, β^i, γ_ij)
     * @param time_derivs Time derivatives of metric (∂g_μν/∂t)
     * @param N_star Conserved baryon number density N* = √γ n Γ
     * @return Gravitational work Λ (energy source per baryon)
     */
    real compute_gravitational_work(
        const GRPrimitive& prim,
        const Metric31& metric,
        const MetricTimeDerivatives& time_derivs,
        real N_star) const;
};

} // namespace grgsph
} // namespace sph
