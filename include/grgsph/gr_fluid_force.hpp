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
 *
 * The force on particle a from particle b consists of:
 * 1. Pressure gradient (from Riemann solution):
 *    dp_i/dt = -Σ_b m_b [P*_{ab} V_a² Ω_a ∇W_a + P*_{ab} V_b² Ω_b ∇W_b]
 *
 * 2. Gravitational source term (from metric gradient):
 *    f_i = (√-g / 2ρ*) T^μν ∂g_μν/∂x^i
 *
 * Key differences from GRSPH:
 * - No artificial viscosity q
 * - Interface pressure P* from Riemann solver
 * - Proper entropy generation at shocks only
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
     * Compute gravitational source term f_i
     *
     * f_i = (√-g / 2ρ*) T^μν ∂g_μν/∂x^i
     *
     * This captures:
     * - Gravitational attraction from metric curvature
     * - Frame-dragging effects (for Kerr)
     *
     * @param prim Primitive variables
     * @param metric Local metric
     * @param derivs Metric derivatives
     * @param rho_star Conserved density
     * @return Gravitational source force vector
     */
    vec_t compute_metric_source(
        const GRPrimitive& prim,
        const Metric31& metric,
        const MetricDerivatives& derivs,
        real rho_star) const;

    /**
     * Normalize SR derivatives (same as SR-GSPH)
     * Scale derivatives by baryon number so dS and de are per-baryon quantities.
     */
    static void normalize_derivatives(SPHParticle& p);
};

} // namespace grgsph
} // namespace sph
