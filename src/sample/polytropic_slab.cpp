/**
 * @file polytropic_slab.cpp
 * @brief 1D Polytropic Self-Gravitating Slab using Planar Lane-Emden Solution
 *
 * This implements a self-gravitating slab in hydrostatic equilibrium for
 * studying diffusive instability in GSPH without grad-h correction.
 *
 * Physics:
 *   - Equation of state: P = K ρ^γ (polytropic)
 *   - Equilibrium: dP/dx = -ρ g(x) where g = 2πG Σ(x)
 *   - Solved by planar Lane-Emden equation: d²θ/dξ² = -θⁿ
 *   - Polytropic index: n = 1/(γ-1)
 *
 * Key difference from isothermal slab:
 *   - Isothermal (γ=1): ρ = ρ₀ sech²(x/H), infinite extent
 *   - Polytropic (γ>1): ρ = ρ₀ θⁿ, FINITE extent at surface ξ₁
 *
 * Reference: Chandrasekhar, "Stellar Structure" (1939), planar case
 */

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>

#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"

namespace sph
{

namespace {

/**
 * @brief Solve planar Lane-Emden equation: d²θ/dξ² = -θⁿ
 *
 * Initial conditions: θ(0) = 1, θ'(0) = 0
 * Integrates until θ → 0 (surface)
 *
 * @param n Polytropic index (n = 1/(γ-1))
 * @param xi_out Output: ξ values
 * @param theta_out Output: θ values
 * @param dtheta_out Output: dθ/dξ values
 * @return ξ₁ (surface location where θ = 0)
 */
double solvePlanarLaneEmden(double n,
                            std::vector<double>& xi_out,
                            std::vector<double>& theta_out,
                            std::vector<double>& dtheta_out)
{
    // RK4 integration
    const double dxi = 1e-4;
    const int max_steps = 1000000;

    xi_out.clear();
    theta_out.clear();
    dtheta_out.clear();

    double xi = 0.0;
    double theta = 1.0;
    double dtheta = 0.0;

    xi_out.push_back(xi);
    theta_out.push_back(theta);
    dtheta_out.push_back(dtheta);

    for (int step = 0; step < max_steps && theta > 0; ++step) {
        // RK4 for system: dθ/dξ = φ, dφ/dξ = -θⁿ (for θ ≥ 0)
        auto f1 = [](double /*xi*/, double /*theta*/, double phi) { return phi; };
        auto f2 = [n](double /*xi*/, double theta_val, double /*phi*/) {
            return (theta_val > 0) ? -std::pow(theta_val, n) : 0.0;
        };

        double k1_theta = dxi * f1(xi, theta, dtheta);
        double k1_phi = dxi * f2(xi, theta, dtheta);

        double k2_theta = dxi * f1(xi + 0.5*dxi, theta + 0.5*k1_theta, dtheta + 0.5*k1_phi);
        double k2_phi = dxi * f2(xi + 0.5*dxi, theta + 0.5*k1_theta, dtheta + 0.5*k1_phi);

        double k3_theta = dxi * f1(xi + 0.5*dxi, theta + 0.5*k2_theta, dtheta + 0.5*k2_phi);
        double k3_phi = dxi * f2(xi + 0.5*dxi, theta + 0.5*k2_theta, dtheta + 0.5*k2_phi);

        double k4_theta = dxi * f1(xi + dxi, theta + k3_theta, dtheta + k3_phi);
        double k4_phi = dxi * f2(xi + dxi, theta + k3_theta, dtheta + k3_phi);

        xi += dxi;
        theta += (k1_theta + 2*k2_theta + 2*k3_theta + k4_theta) / 6.0;
        dtheta += (k1_phi + 2*k2_phi + 2*k3_phi + k4_phi) / 6.0;

        if (theta < 0) {
            // Linear interpolation to find exact surface
            double xi_prev = xi_out.back();
            double theta_prev = theta_out.back();
            double xi_surface = xi_prev - theta_prev * dxi / (theta - theta_prev);
            xi = xi_surface;
            theta = 0.0;
        }

        xi_out.push_back(xi);
        theta_out.push_back(std::max(0.0, theta));
        dtheta_out.push_back(dtheta);
    }

    return xi_out.back();  // ξ₁
}

/**
 * @brief Compute cumulative mass function M(ξ) for equal-mass particle placement
 *
 * For planar geometry: M(ξ) = ∫₀^ξ ρ dξ' = ρ_c α ∫₀^ξ θⁿ dξ'
 * where α is the length scale.
 *
 * We return the normalized version: F(ξ) = M(ξ) / M_total
 */
std::vector<double> computeCumulativeMass(const std::vector<double>& xi,
                                          const std::vector<double>& theta,
                                          double n)
{
    std::vector<double> M(xi.size(), 0.0);

    // Trapezoidal integration
    for (size_t i = 1; i < xi.size(); ++i) {
        double dxi = xi[i] - xi[i-1];
        double rho_avg = 0.5 * (std::pow(theta[i-1], n) + std::pow(theta[i], n));
        M[i] = M[i-1] + rho_avg * dxi;
    }

    // Normalize to [0, 1]
    double M_total = M.back();
    for (auto& m : M) {
        m /= M_total;
    }

    return M;
}

/**
 * @brief Linearly interpolate to find ξ for given normalized mass F
 */
double interpolateXi(double F, const std::vector<double>& xi,
                     const std::vector<double>& M_norm)
{
    // Binary search for bracket
    auto it = std::lower_bound(M_norm.begin(), M_norm.end(), F);
    if (it == M_norm.begin()) return xi.front();
    if (it == M_norm.end()) return xi.back();

    size_t idx = std::distance(M_norm.begin(), it);
    double F0 = M_norm[idx-1];
    double F1 = M_norm[idx];
    double xi0 = xi[idx-1];
    double xi1 = xi[idx];

    // Linear interpolation
    double t = (F - F0) / (F1 - F0);
    return xi0 + t * (xi1 - xi0);
}

/**
 * @brief Linearly interpolate θ at given ξ
 */
double interpolateTheta(double xi_val, const std::vector<double>& xi,
                        const std::vector<double>& theta)
{
    auto it = std::lower_bound(xi.begin(), xi.end(), xi_val);
    if (it == xi.begin()) return theta.front();
    if (it == xi.end()) return theta.back();

    size_t idx = std::distance(xi.begin(), it);
    double xi0 = xi[idx-1];
    double xi1 = xi[idx];
    double theta0 = theta[idx-1];
    double theta1 = theta[idx];

    double t = (xi_val - xi0) / (xi1 - xi0);
    return theta0 + t * (theta1 - theta0);
}

} // anonymous namespace


void Solver::make_polytropic_slab()
{
#if DIM != 1
    THROW_ERROR("Polytropic slab requires DIM == 1");
#else

    // ========================================================================
    // Read parameters
    // ========================================================================
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real rho_center = boost::any_cast<double>(m_sample_parameters["rho_center"]);
    const real K = boost::any_cast<double>(m_sample_parameters["K"]);

    const real G = m_param->gravity.constant;
    const real gamma = m_param->physics.gamma;

    // Polytropic index
    const real n = 1.0 / (gamma - 1.0);

    std::cout << "\n=== Creating 1D Polytropic Self-Gravitating Slab ===" << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "  N = " << N << " particles" << std::endl;
    std::cout << "  ρ_center = " << rho_center << std::endl;
    std::cout << "  K = " << K << std::endl;
    std::cout << "  G = " << G << std::endl;
    std::cout << "  γ = " << gamma << std::endl;
    std::cout << "  n = 1/(γ-1) = " << n << std::endl;

    // ========================================================================
    // Solve planar Lane-Emden equation
    // ========================================================================
    std::vector<double> xi_arr, theta_arr, dtheta_arr;
    double xi_1 = solvePlanarLaneEmden(n, xi_arr, theta_arr, dtheta_arr);

    std::cout << "\nLane-Emden solution:" << std::endl;
    std::cout << "  Surface at ξ₁ = " << xi_1 << std::endl;
    std::cout << "  Data points: " << xi_arr.size() << std::endl;

    // ========================================================================
    // Physical scaling
    // ========================================================================
    // For planar polytrope:
    //   ξ = x / α where α² = (n+1) K ρ_c^(1/n - 1) / (4πG)
    //   But for planar geometry (not spherical), the relation is:
    //   α² = K (n+1) ρ_c^(1-n) / (2πG)
    //
    // Central pressure: P_c = K ρ_c^γ
    // Central sound speed: c_s² = γ P_c / ρ_c = γ K ρ_c^(γ-1)

    const real alpha_sq = K * (n + 1.0) * std::pow(rho_center, 1.0 - n) / (2.0 * M_PI * G);
    const real alpha = std::sqrt(alpha_sq);

    const real x_surface = alpha * xi_1;  // Physical half-width of slab
    const real P_center = K * std::pow(rho_center, gamma);
    const real c_s_center = std::sqrt(gamma * P_center / rho_center);

    std::cout << "\nPhysical parameters:" << std::endl;
    std::cout << "  Length scale α = " << alpha << std::endl;
    std::cout << "  Slab half-width x₁ = α·ξ₁ = " << x_surface << std::endl;
    std::cout << "  Central pressure P_c = " << P_center << std::endl;
    std::cout << "  Central sound speed c_s = " << c_s_center << std::endl;

    // ========================================================================
    // Total mass (per unit area, since 1D is really a slab)
    // M = 2 ∫₀^x₁ ρ dx = 2 ρ_c α ∫₀^ξ₁ θⁿ dξ
    // ========================================================================
    double mass_integral = 0.0;
    for (size_t i = 1; i < xi_arr.size(); ++i) {
        double dxi = xi_arr[i] - xi_arr[i-1];
        double rho_avg = 0.5 * (std::pow(theta_arr[i-1], n) + std::pow(theta_arr[i], n));
        mass_integral += rho_avg * dxi;
    }
    const real M_half = rho_center * alpha * mass_integral;  // Mass in x > 0
    const real M_total = 2.0 * M_half;
    const real mass_per_particle = M_total / N;

    std::cout << "  Total mass = " << M_total << std::endl;
    std::cout << "  Mass per particle = " << mass_per_particle << std::endl;

    // ========================================================================
    // Compute cumulative mass for equal-mass particle placement
    // ========================================================================
    std::vector<double> M_norm = computeCumulativeMass(xi_arr, theta_arr, n);

    // ========================================================================
    // Place particles with equal mass
    // ========================================================================
    std::vector<SPHParticle> particles(N);

    // Place N/2 particles on each side (symmetric about x=0)
    // For particle i, its cumulative mass fraction is F_i = (i + 0.5) / N
    // Map this to position via the cumulative mass function

    for (int i = 0; i < N; ++i) {
        auto& p_i = particles[i];

        // Normalized cumulative mass for this particle
        // Particles are numbered 0 to N-1, spanning the full slab
        // F = 0 at x = -x_surface, F = 1 at x = +x_surface
        const real F_i = (i + 0.5) / N;

        // Map F to ξ in the range [0, ξ₁] for each half
        // F ∈ [0, 0.5] → x ∈ [-x_surface, 0]
        // F ∈ [0.5, 1] → x ∈ [0, +x_surface]
        real xi_val, x;
        if (F_i < 0.5) {
            // Left half: F ∈ [0, 0.5] maps to ξ ∈ [ξ₁, 0] (reversed)
            real F_half = 0.5 - F_i;  // F_half ∈ [0, 0.5]
            real F_local = 2.0 * F_half;  // F_local ∈ [0, 1]
            xi_val = interpolateXi(F_local, xi_arr, M_norm);
            x = -alpha * xi_val;
        } else {
            // Right half: F ∈ [0.5, 1] maps to ξ ∈ [0, ξ₁]
            real F_half = F_i - 0.5;  // F_half ∈ [0, 0.5]
            real F_local = 2.0 * F_half;  // F_local ∈ [0, 1]
            xi_val = interpolateXi(F_local, xi_arr, M_norm);
            x = alpha * xi_val;
        }

        // Position
        p_i.pos[0] = x;

        // Velocity (start at rest)
        p_i.vel[0] = 0.0;

        // Density from Lane-Emden solution: ρ = ρ_c θⁿ
        real theta_val = interpolateTheta(std::abs(x) / alpha, xi_arr, theta_arr);
        theta_val = std::max(theta_val, 1e-10);  // Avoid exact zero
        p_i.dens = rho_center * std::pow(theta_val, n);

        // Pressure from polytropic EOS: P = K ρ^γ
        p_i.pres = K * std::pow(p_i.dens, gamma);

        // Internal energy for adiabatic EOS: P = (γ-1) ρ u
        p_i.ene = p_i.pres / ((gamma - 1.0) * p_i.dens);

        // Sound speed: c_s = √(γ P / ρ)
        p_i.sound = std::sqrt(gamma * p_i.pres / p_i.dens);

        // Mass
        p_i.mass = mass_per_particle;

        // ID
        p_i.id = i;

        // Initial smoothing length estimate
        // For 1D: h ≈ N_ngb * m / (2 * ρ)
        const int n_ngb = m_param->physics.neighbor_number;
        p_i.sml = n_ngb * p_i.mass / (2.0 * p_i.dens);
    }

    // ========================================================================
    // Compute initial gravitational acceleration
    // g(x) = -2πG ∫₀^|x| ρ dx' * sign(x) = -2πG ρ_c α ∫₀^ξ θⁿ dξ' * sign(x)
    // ========================================================================
    for (int i = 0; i < N; ++i) {
        auto& p_i = particles[i];
        const real x = p_i.pos[0];
        const real xi_pos = std::abs(x) / alpha;

        // Integrate θⁿ from 0 to ξ
        real Sigma = 0.0;
        for (size_t j = 1; j < xi_arr.size() && xi_arr[j] <= xi_pos; ++j) {
            double dxi = xi_arr[j] - xi_arr[j-1];
            double rho_avg = 0.5 * (std::pow(theta_arr[j-1], n) + std::pow(theta_arr[j], n));
            Sigma += rho_avg * dxi;
        }
        Sigma *= rho_center * alpha;  // Physical surface density

        // Gravitational acceleration (points toward center)
        real sign_x = (x > 0) ? 1.0 : ((x < 0) ? -1.0 : 0.0);
        p_i.grav_acc[0] = -2.0 * M_PI * G * Sigma * sign_x;

        // Gravitational potential (for diagnostics)
        // φ(x) = -2πG ∫ |x - x'| ρ(x') dx'
        // This is more complex; use a simpler approximation
        p_i.phi = -2.0 * M_PI * G * Sigma * std::abs(x);
    }

    // ========================================================================
    // Report setup and verify equilibrium
    // ========================================================================
    std::cout << "\nParticle distribution:" << std::endl;
    std::cout << "  x range: [" << particles[0].pos[0] << ", "
              << particles[N-1].pos[0] << "]" << std::endl;
    std::cout << "  ρ range: [" << particles[N/2].dens << " (center), "
              << particles[0].dens << " (edge)]" << std::endl;

    // Estimate dynamical time
    const real t_dyn = x_surface / c_s_center;
    std::cout << "\nTimescales:" << std::endl;
    std::cout << "  t_dyn = x₁/c_s = " << t_dyn << std::endl;

    // Estimate expected growth rate (theory)
    // Γ = ε * c_s * h / L² where L ~ x_surface
    const real h_avg = particles[N/2].sml;
    const real epsilon = 0.4;  // Expected force error without grad-h
    const real Gamma_theory = epsilon * c_s_center * h_avg / (x_surface * x_surface);
    std::cout << "  h_avg (center) = " << h_avg << std::endl;
    std::cout << "  Expected Γ (without grad-h) ≈ " << Gamma_theory << std::endl;
    std::cout << "  e-folding time ≈ " << 1.0 / Gamma_theory << " = "
              << (1.0 / Gamma_theory) / t_dyn << " t_dyn" << std::endl;

    // Verify hydrostatic equilibrium at a few points
    // Equilibrium: dP/dx = ρ * g (where g includes sign toward center)
    std::cout << "\nEquilibrium check (dP/dx vs ρg):" << std::endl;
    for (int idx : {N/4, N/2, 3*N/4}) {
        if (idx <= 0 || idx >= N-1) continue;
        const auto& p = particles[idx];
        const auto& p_prev = particles[idx-1];
        const auto& p_next = particles[idx+1];

        real dPdx = (p_next.pres - p_prev.pres) / (p_next.pos[0] - p_prev.pos[0]);
        real rho_g = p.dens * p.grav_acc[0];

        std::cout << "  x=" << p.pos[0] << ": dP/dx=" << dPdx
                  << ", ρg=" << rho_g
                  << ", error=" << std::abs(dPdx - rho_g) / std::max(std::abs(rho_g), 1e-10) * 100 << "%"
                  << std::endl;
    }

    std::cout << "\n=== Polytropic Slab Setup Complete ===" << std::endl;

    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());

#endif
}

} // namespace sph
