/**
 * @file polytropic_slab_2d.cpp
 * @brief 2D Polytropic Self-Gravitating Slab using Planar Lane-Emden Solution
 *
 * This implements a self-gravitating slab in hydrostatic equilibrium in 2D,
 * where the density varies in the y-direction according to the planar Lane-Emden
 * solution, and is uniform in the x-direction.
 *
 * Geometry:
 *   - Infinite slab in x-z plane (uniform in x)
 *   - Density profile ρ(y) from Lane-Emden planar equation
 *   - Gravity acts only in y-direction (1D gravity in 2D space)
 *
 * Physics:
 *   - Equation of state: P = K ρ^γ (polytropic)
 *   - Equilibrium: dP/dy = -ρ g(y) where g = 2πG Σ(y)
 *   - Solved by planar Lane-Emden equation: d²θ/dξ² = -θⁿ
 *   - Polytropic index: n = 1/(γ-1)
 *
 * For γ = 5/3, n = 1.5 (Lane-Emden n=3/2)
 *
 * Reference: Chandrasekhar, "Stellar Structure" (1939), planar case
 */

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>

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
double solvePlanarLaneEmden2D(double n,
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
 */
std::vector<double> computeCumulativeMass2D(const std::vector<double>& xi,
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
double interpolateXi2D(double F, const std::vector<double>& xi,
                       const std::vector<double>& M_norm)
{
    auto it = std::lower_bound(M_norm.begin(), M_norm.end(), F);
    if (it == M_norm.begin()) return xi.front();
    if (it == M_norm.end()) return xi.back();

    size_t idx = std::distance(M_norm.begin(), it);
    double F0 = M_norm[idx-1];
    double F1 = M_norm[idx];
    double xi0 = xi[idx-1];
    double xi1 = xi[idx];

    double t = (F - F0) / (F1 - F0);
    return xi0 + t * (xi1 - xi0);
}

/**
 * @brief Linearly interpolate θ at given ξ
 */
double interpolateTheta2D(double xi_val, const std::vector<double>& xi,
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


void Solver::make_polytropic_slab_2d()
{
#if DIM != 2
    THROW_ERROR("Polytropic slab 2D requires DIM == 2");
#else

    // ========================================================================
    // Read parameters
    // ========================================================================
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);  // N × N particles
    const real rho_center = boost::any_cast<double>(m_sample_parameters["rho_center"]);
    const real K = boost::any_cast<double>(m_sample_parameters["K"]);
    const real L_x = boost::any_cast<double>(m_sample_parameters["L_x"]);  // Width in x

    const real G = m_param->gravity.constant;
    const real gamma = m_param->physics.gamma;

    // Polytropic index
    const real n = 1.0 / (gamma - 1.0);

    std::cout << "\n=== Creating 2D Polytropic Self-Gravitating Slab ===" << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "  N = " << N << " × " << N << " = " << N*N << " particles" << std::endl;
    std::cout << "  ρ_center = " << rho_center << std::endl;
    std::cout << "  K = " << K << std::endl;
    std::cout << "  L_x = " << L_x << std::endl;
    std::cout << "  G = " << G << std::endl;
    std::cout << "  γ = " << gamma << std::endl;
    std::cout << "  n = 1/(γ-1) = " << n << std::endl;

    // ========================================================================
    // Solve planar Lane-Emden equation
    // ========================================================================
    std::vector<double> xi_arr, theta_arr, dtheta_arr;
    double xi_1 = solvePlanarLaneEmden2D(n, xi_arr, theta_arr, dtheta_arr);

    std::cout << "\nLane-Emden solution:" << std::endl;
    std::cout << "  Surface at ξ₁ = " << xi_1 << std::endl;
    std::cout << "  Data points: " << xi_arr.size() << std::endl;

    // ========================================================================
    // Physical scaling
    // α² = K (n+1) ρ_c^(1-n) / (2πG)
    // ========================================================================
    const real alpha_sq = K * (n + 1.0) * std::pow(rho_center, 1.0 - n) / (2.0 * M_PI * G);
    const real alpha = std::sqrt(alpha_sq);

    const real y_surface = alpha * xi_1;  // Physical half-width of slab in y
    const real P_center = K * std::pow(rho_center, gamma);
    const real c_s_center = std::sqrt(gamma * P_center / rho_center);

    std::cout << "\nPhysical parameters:" << std::endl;
    std::cout << "  Length scale α = " << alpha << std::endl;
    std::cout << "  Slab half-width y₁ = α·ξ₁ = " << y_surface << std::endl;
    std::cout << "  Central pressure P_c = " << P_center << std::endl;
    std::cout << "  Central sound speed c_s = " << c_s_center << std::endl;

    // ========================================================================
    // Total mass (per unit length in x)
    // M = 2 ∫₀^y₁ ρ dy = 2 ρ_c α ∫₀^ξ₁ θⁿ dξ
    // ========================================================================
    double mass_integral = 0.0;
    for (size_t i = 1; i < xi_arr.size(); ++i) {
        double dxi = xi_arr[i] - xi_arr[i-1];
        double rho_avg = 0.5 * (std::pow(theta_arr[i-1], n) + std::pow(theta_arr[i], n));
        mass_integral += rho_avg * dxi;
    }
    const real sigma_half = rho_center * alpha * mass_integral;  // Surface density in y > 0
    const real sigma_total = 2.0 * sigma_half;  // Total surface density
    const real M_total = sigma_total * L_x;  // Total mass
    
    const int N_total = N * N;
    const real mass_per_particle = M_total / N_total;

    std::cout << "  Surface density σ = " << sigma_total << std::endl;
    std::cout << "  Total mass = " << M_total << std::endl;
    std::cout << "  Mass per particle = " << mass_per_particle << std::endl;

    // ========================================================================
    // Compute cumulative mass for equal-mass particle placement in y
    // ========================================================================
    std::vector<double> M_norm = computeCumulativeMass2D(xi_arr, theta_arr, n);

    // ========================================================================
    // Place particles on a grid: uniform in x, Lane-Emden distributed in y
    // ========================================================================
    const int N_x = N;  // Particles in x
    const int N_y = N;  // Particles in y (spanning full slab from -y_surface to +y_surface)
    
    std::vector<SPHParticle> particles;
    particles.reserve(N_x * N_y);

    const real dx = L_x / N_x;  // Uniform spacing in x

    int id = 0;
    for (int iy = 0; iy < N_y; ++iy) {
        // Determine y-position using equal-mass distribution
        // F_i = (iy + 0.5) / N_y maps to position
        const real F_i = (iy + 0.5) / N_y;

        real y;
        if (F_i < 0.5) {
            // Left half: F ∈ [0, 0.5] maps to y ∈ [-y_surface, 0]
            real F_local = 2.0 * (0.5 - F_i);  // F_local ∈ [0, 1]
            real xi_val = interpolateXi2D(F_local, xi_arr, M_norm);
            y = -alpha * xi_val;
        } else {
            // Right half: F ∈ [0.5, 1] maps to y ∈ [0, +y_surface]
            real F_local = 2.0 * (F_i - 0.5);  // F_local ∈ [0, 1]
            real xi_val = interpolateXi2D(F_local, xi_arr, M_norm);
            y = alpha * xi_val;
        }

        // Density from Lane-Emden solution: ρ = ρ_c θⁿ
        real theta_val = interpolateTheta2D(std::abs(y) / alpha, xi_arr, theta_arr);
        theta_val = std::max(theta_val, 1e-10);  // Avoid exact zero
        real dens = rho_center * std::pow(theta_val, n);

        // Pressure from polytropic EOS: P = K ρ^γ
        real pres = K * std::pow(dens, gamma);

        // Internal energy for adiabatic EOS: P = (γ-1) ρ u
        real ene = pres / ((gamma - 1.0) * dens);

        // Sound speed: c_s = √(γ P / ρ)
        real sound = std::sqrt(gamma * pres / dens);

        // Initial smoothing length estimate
        // For 2D: h ≈ (N_ngb * m / (π ρ))^(1/2)
        const int n_ngb = m_param->physics.neighbor_number;
        real sml = std::sqrt(n_ngb * mass_per_particle / (M_PI * dens));

        // Place N_x particles along x at this y
        for (int ix = 0; ix < N_x; ++ix) {
            SPHParticle p;

            // Position: uniform in x, Lane-Emden in y
            p.pos[0] = -L_x / 2.0 + (ix + 0.5) * dx;
            p.pos[1] = y;

            // Velocity (start at rest)
            p.vel[0] = 0.0;
            p.vel[1] = 0.0;

            // Thermodynamic quantities
            p.dens = dens;
            p.pres = pres;
            p.ene = ene;
            p.sound = sound;

            // Mass
            p.mass = mass_per_particle;

            // ID
            p.id = id++;

            // Smoothing length
            p.sml = sml;

            particles.push_back(p);
        }
    }

    std::cout << "\nParticle distribution:" << std::endl;
    std::cout << "  Total particles: " << particles.size() << std::endl;
    std::cout << "  x range: [" << particles[0].pos[0] << ", "
              << particles[N_x-1].pos[0] << "]" << std::endl;
    std::cout << "  y range: [" << particles[0].pos[1] << ", "
              << particles.back().pos[1] << "]" << std::endl;

    // ========================================================================
    // Compute initial gravitational acceleration
    // g(y) = -2πG ∫₀^|y| ρ dy' * sign(y)  (only y-component, x=0)
    // ========================================================================
    for (auto& p_i : particles) {
        const real y = p_i.pos[1];
        const real xi_pos = std::abs(y) / alpha;

        // Integrate θⁿ from 0 to ξ
        real Sigma = 0.0;
        for (size_t j = 1; j < xi_arr.size() && xi_arr[j] <= xi_pos; ++j) {
            double dxi = xi_arr[j] - xi_arr[j-1];
            double rho_avg = 0.5 * (std::pow(theta_arr[j-1], n) + std::pow(theta_arr[j], n));
            Sigma += rho_avg * dxi;
        }
        Sigma *= rho_center * alpha;  // Physical surface density

        // Gravitational acceleration (points toward center, y-direction only)
        real sign_y = (y > 0) ? 1.0 : ((y < 0) ? -1.0 : 0.0);
        p_i.grav_acc[0] = 0.0;  // No x-gravity
        p_i.grav_acc[1] = -2.0 * M_PI * G * Sigma * sign_y;

        // Gravitational potential (for diagnostics)
        p_i.phi = -2.0 * M_PI * G * Sigma * std::abs(y);
    }

    // ========================================================================
    // Report setup
    // ========================================================================
    // Estimate dynamical time
    const real t_dyn = y_surface / c_s_center;
    std::cout << "\nTimescales:" << std::endl;
    std::cout << "  t_dyn = y₁/c_s = " << t_dyn << std::endl;

    std::cout << "\n=== 2D Polytropic Slab Setup Complete ===" << std::endl;
    std::cout << "NOTE: Use 'useKernelGravityPlanar2D: true' for correct gravity!" << std::endl;

    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());

#endif
}

} // namespace sph
