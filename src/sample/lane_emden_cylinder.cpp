/**
 * @file lane_emden_cylinder.cpp
 * @brief 3D Cylindrical Lane-Emden Solution (Infinite Cylinder)
 *
 * This implements a self-gravitating infinite cylinder in hydrostatic equilibrium,
 * where the density varies radially in the xy-plane according to the cylindrical
 * Lane-Emden solution, and is uniform in the z-direction.
 *
 * Geometry:
 *   - Infinite cylinder along z-axis
 *   - Density profile ρ(r) from cylindrical Lane-Emden equation
 *   - Gravity acts only radially in xy-plane (2D gravity in 3D space)
 *
 * Physics:
 *   - Equation of state: P = K ρ^γ (polytropic)
 *   - Equilibrium: dP/dr = -ρ g(r) where g = 2G M(r) / r
 *   - Cylindrical Lane-Emden: (1/ξ) d/dξ (ξ dθ/dξ) = -θⁿ
 *   - Polytropic index: n = 1/(γ-1)
 *
 * For γ = 5/3, n = 1.5 (Lane-Emden n=3/2)
 *
 * Key difference from spherical Lane-Emden:
 *   - Cylindrical: Laplacian is (1/r) d/dr (r d/dr) = d²/dr² + (1/r) d/dr
 *   - Spherical: Laplacian is (1/r²) d/dr (r² d/dr) = d²/dr² + (2/r) d/dr
 *
 * Reference: Chandrasekhar, "Stellar Structure" (1939), cylindrical case
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
 * @brief Solve cylindrical Lane-Emden equation: (1/ξ) d/dξ (ξ dθ/dξ) = -θⁿ
 *
 * This is equivalent to: d²θ/dξ² + (1/ξ) dθ/dξ = -θⁿ
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
double solveCylindricalLaneEmden(double n,
                                  std::vector<double>& xi_out,
                                  std::vector<double>& theta_out,
                                  std::vector<double>& dtheta_out)
{
    // RK4 integration
    const double dxi = 1e-4;
    const int max_steps = 1000000;
    const double xi_start = 1e-6;  // Small offset to avoid 1/ξ singularity

    xi_out.clear();
    theta_out.clear();
    dtheta_out.clear();

    // Initial conditions: θ(0) = 1, θ'(0) = 0
    // Near ξ = 0, series solution: θ ≈ 1 - ξ²/4 + O(ξ⁴) for n=1.5
    double xi = xi_start;
    double theta = 1.0 - xi_start * xi_start / 4.0;  // First correction term
    double dtheta = -xi_start / 2.0;  // Derivative at small ξ

    xi_out.push_back(0.0);      // Store exact ξ=0 point
    theta_out.push_back(1.0);
    dtheta_out.push_back(0.0);

    xi_out.push_back(xi);
    theta_out.push_back(theta);
    dtheta_out.push_back(dtheta);

    for (int step = 0; step < max_steps && theta > 0; ++step) {
        // RK4 for system:
        // dθ/dξ = φ
        // dφ/dξ = -θⁿ - φ/ξ  (from cylindrical Lane-Emden)
        auto f1 = [](double /*xi*/, double /*theta*/, double phi) { return phi; };
        auto f2 = [n](double xi_val, double theta_val, double phi_val) {
            if (theta_val <= 0) return 0.0;
            return -std::pow(theta_val, n) - phi_val / xi_val;
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
 * For cylindrical geometry: M(ξ) = 2π ∫₀^ξ ρ ξ' dξ' = 2π ρ_c α² ∫₀^ξ θⁿ ξ' dξ'
 */
std::vector<double> computeCumulativeMassCylinder(const std::vector<double>& xi,
                                                   const std::vector<double>& theta,
                                                   double n)
{
    std::vector<double> M(xi.size(), 0.0);

    // Trapezoidal integration of ξ θⁿ
    for (size_t i = 1; i < xi.size(); ++i) {
        double dxi = xi[i] - xi[i-1];
        double xi_mid = 0.5 * (xi[i] + xi[i-1]);
        double theta_i_n = (theta[i] > 0) ? std::pow(theta[i], n) : 0.0;
        double theta_im1_n = (theta[i-1] > 0) ? std::pow(theta[i-1], n) : 0.0;
        double integrand_avg = 0.5 * (xi[i-1] * theta_im1_n + xi[i] * theta_i_n);
        M[i] = M[i-1] + integrand_avg * dxi;
    }

    // Normalize to [0, 1]
    double M_total = M.back();
    if (M_total > 0) {
        for (auto& m : M) {
            m /= M_total;
        }
    }

    return M;
}

/**
 * @brief Linearly interpolate to find ξ for given normalized mass F
 */
double interpolateXiCylinder(double F, const std::vector<double>& xi,
                              const std::vector<double>& M_norm)
{
    if (F <= 0) return xi.front();
    if (F >= 1) return xi.back();

    auto it = std::lower_bound(M_norm.begin(), M_norm.end(), F);
    if (it == M_norm.begin()) return xi.front();
    if (it == M_norm.end()) return xi.back();

    size_t idx = std::distance(M_norm.begin(), it);
    if (idx == 0) return xi.front();
    
    double F0 = M_norm[idx-1];
    double F1 = M_norm[idx];
    double xi0 = xi[idx-1];
    double xi1 = xi[idx];

    double t = (F - F0) / (F1 - F0 + 1e-30);
    return xi0 + t * (xi1 - xi0);
}

/**
 * @brief Linearly interpolate θ at given ξ
 */
double interpolateThetaCylinder(double xi_val, const std::vector<double>& xi,
                                 const std::vector<double>& theta)
{
    if (xi_val <= xi.front()) return theta.front();
    if (xi_val >= xi.back()) return theta.back();

    auto it = std::lower_bound(xi.begin(), xi.end(), xi_val);
    if (it == xi.begin()) return theta.front();
    if (it == xi.end()) return theta.back();

    size_t idx = std::distance(xi.begin(), it);
    if (idx == 0) return theta.front();
    
    double xi0 = xi[idx-1];
    double xi1 = xi[idx];
    double theta0 = theta[idx-1];
    double theta1 = theta[idx];

    double t = (xi_val - xi0) / (xi1 - xi0 + 1e-30);
    return theta0 + t * (theta1 - theta0);
}

} // anonymous namespace


void Solver::make_lane_emden_cylinder()
{
#if DIM != 3
    THROW_ERROR("Lane-Emden cylinder requires DIM == 3");
#else

    // ========================================================================
    // Read parameters
    // ========================================================================
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);  // Nominal N³ particles
    const real R = boost::any_cast<double>(m_sample_parameters["R"]);  // Cylinder radius
    const real L_z = boost::any_cast<double>(m_sample_parameters["L_z"]);  // Length in z
    const real M_total = boost::any_cast<double>(m_sample_parameters["M_total"]);

    const real G = m_param->gravity.constant;
    const real gamma = m_param->physics.gamma;

    // Polytropic index
    const real n = 1.0 / (gamma - 1.0);

    std::cout << "\n=== Creating 3D Lane-Emden Cylinder ===" << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "  N = " << N << " (nominal N³ = " << N*N*N << " particles)" << std::endl;
    std::cout << "  R = " << R << std::endl;
    std::cout << "  L_z = " << L_z << std::endl;
    std::cout << "  M_total = " << M_total << std::endl;
    std::cout << "  G = " << G << std::endl;
    std::cout << "  γ = " << gamma << std::endl;
    std::cout << "  n = 1/(γ-1) = " << n << std::endl;

    // ========================================================================
    // Solve cylindrical Lane-Emden equation
    // ========================================================================
    std::vector<double> xi_arr, theta_arr, dtheta_arr;
    double xi_1 = solveCylindricalLaneEmden(n, xi_arr, theta_arr, dtheta_arr);

    std::cout << "\nCylindrical Lane-Emden solution:" << std::endl;
    std::cout << "  Surface at ξ₁ = " << xi_1 << std::endl;
    std::cout << "  |dθ/dξ|_{ξ₁} = " << std::abs(dtheta_arr.back()) << std::endl;
    std::cout << "  Data points: " << xi_arr.size() << std::endl;

    // ========================================================================
    // Physical scaling
    // For cylindrical polytrope:
    //   ξ = r / α where α² = K (n+1) ρ_c^(1-n) / (2πG)
    //   Mass per unit length: λ = 2π α² ρ_c ξ₁ |dθ/dξ|_{ξ₁}
    // ========================================================================
    const real alpha = R / xi_1;  // Choose α so that r(ξ₁) = R
    
    // Mass per unit length from Lane-Emden
    // λ = 2π α² ρ_c ∫₀^ξ₁ ξ θⁿ dξ = 2π α² ρ_c × ξ₁ |dθ/dξ|_{ξ₁}  (Lane-Emden identity)
    // We want total mass M_total = λ × L_z
    // So: ρ_c = M_total / (2π α² L_z × ξ₁ |dθ/dξ|_{ξ₁})
    const real lambda_factor = xi_1 * std::abs(dtheta_arr.back());
    const real rho_center = M_total / (2.0 * M_PI * alpha * alpha * L_z * lambda_factor);

    // Polytropic constant K from equilibrium condition
    // α² = K (n+1) ρ_c^(1-n) / (2πG)
    const real K = 2.0 * M_PI * G * alpha * alpha * std::pow(rho_center, n - 1.0) / (n + 1.0);

    const real P_center = K * std::pow(rho_center, gamma);
    const real c_s_center = std::sqrt(gamma * P_center / rho_center);

    std::cout << "\nPhysical parameters:" << std::endl;
    std::cout << "  Length scale α = " << alpha << std::endl;
    std::cout << "  ρ_center = " << rho_center << std::endl;
    std::cout << "  K = " << K << std::endl;
    std::cout << "  Central pressure P_c = " << P_center << std::endl;
    std::cout << "  Central sound speed c_s = " << c_s_center << std::endl;

    // ========================================================================
    // Compute cumulative mass for equal-mass particle placement
    // ========================================================================
    std::vector<double> M_norm = computeCumulativeMassCylinder(xi_arr, theta_arr, n);

    // ========================================================================
    // Place particles: random in azimuth and z, radii mapped by mass
    // ========================================================================
    const int N_particles = N * N * N;
    const real mass_per_particle = M_total / N_particles;

    std::cout << "\nParticle placement:" << std::endl;
    std::cout << "  Target particles: " << N_particles << std::endl;
    std::cout << "  Mass per particle: " << mass_per_particle << std::endl;

    std::vector<SPHParticle> particles;
    particles.reserve(N_particles);

    // Use fixed seed for reproducibility
    std::mt19937 gen(42);
    std::uniform_real_distribution<real> dis_phi(0.0, 2.0 * M_PI);  // Azimuthal angle
    std::uniform_real_distribution<real> dis_z(-L_z / 2.0, L_z / 2.0);  // z-position
    std::uniform_real_distribution<real> dis_mass(0.0, 1.0);  // Mass fraction for radius

    for (int i = 0; i < N_particles; ++i) {
        SPHParticle p;

        // Random azimuthal angle
        real phi = dis_phi(gen);

        // Random z-position (uniform along cylinder)
        real z = dis_z(gen);

        // Radius from inverse cumulative mass function
        // For equal-mass particles, we sample uniformly in mass
        real F_r = dis_mass(gen);
        real xi_val = interpolateXiCylinder(F_r, xi_arr, M_norm);
        real r = alpha * xi_val;

        // Cap at surface
        if (r > R) r = R * 0.99;

        // Convert to Cartesian
        p.pos[0] = r * std::cos(phi);
        p.pos[1] = r * std::sin(phi);
        p.pos[2] = z;

        // Velocity (start at rest)
        p.vel[0] = 0.0;
        p.vel[1] = 0.0;
        p.vel[2] = 0.0;

        // Density from Lane-Emden solution: ρ = ρ_c θⁿ
        real theta_val = interpolateThetaCylinder(xi_val, xi_arr, theta_arr);
        theta_val = std::max(theta_val, 1e-10);  // Avoid exact zero
        p.dens = rho_center * std::pow(theta_val, n);

        // Pressure from polytropic EOS: P = K ρ^γ
        p.pres = K * std::pow(p.dens, gamma);

        // Internal energy for adiabatic EOS: P = (γ-1) ρ u
        p.ene = p.pres / ((gamma - 1.0) * p.dens);

        // Sound speed: c_s = √(γ P / ρ)
        p.sound = std::sqrt(gamma * p.pres / p.dens);

        // Mass
        p.mass = mass_per_particle;

        // ID
        p.id = i;

        // Initial smoothing length estimate
        // For 3D: h ≈ (N_ngb * m / (4π/3 ρ))^(1/3)
        const int n_ngb = m_param->physics.neighbor_number;
        p.sml = std::pow(3.0 * n_ngb * p.mass / (4.0 * M_PI * p.dens), 1.0 / 3.0);

        particles.push_back(p);
    }

    // ========================================================================
    // Compute initial gravitational acceleration
    // For cylinder: g(r) = -2G λ(r) / r where λ(r) = enclosed mass per unit length
    // λ(r) = 2π α² ρ_c ∫₀^ξ ξ' θⁿ dξ'
    // ========================================================================
    for (auto& p_i : particles) {
        const real x = p_i.pos[0];
        const real y = p_i.pos[1];
        const real r = std::sqrt(x * x + y * y);

        if (r < 1e-30) {
            p_i.grav_acc[0] = 0.0;
            p_i.grav_acc[1] = 0.0;
            p_i.grav_acc[2] = 0.0;
            p_i.phi = 0.0;
            continue;
        }

        const real xi_pos = r / alpha;

        // Integrate ξ θⁿ from 0 to ξ
        real enclosed_integral = 0.0;
        for (size_t j = 1; j < xi_arr.size() && xi_arr[j] <= xi_pos; ++j) {
            double dxi = xi_arr[j] - xi_arr[j-1];
            double xi_mid = 0.5 * (xi_arr[j] + xi_arr[j-1]);
            double theta_j_n = (theta_arr[j] > 0) ? std::pow(theta_arr[j], n) : 0.0;
            double theta_jm1_n = (theta_arr[j-1] > 0) ? std::pow(theta_arr[j-1], n) : 0.0;
            double integrand_avg = 0.5 * (xi_arr[j-1] * theta_jm1_n + xi_arr[j] * theta_j_n);
            enclosed_integral += integrand_avg * dxi;
        }

        // Enclosed mass per unit length
        real lambda_enclosed = 2.0 * M_PI * alpha * alpha * rho_center * enclosed_integral;

        // Gravitational acceleration (radial, in xy-plane only)
        real g_mag = 2.0 * G * lambda_enclosed / r;
        p_i.grav_acc[0] = -g_mag * x / r;
        p_i.grav_acc[1] = -g_mag * y / r;
        p_i.grav_acc[2] = 0.0;  // No z-gravity for infinite cylinder

        // Gravitational potential (for diagnostics)
        // φ = -2G λ ln(r) for point mass at r
        p_i.phi = -2.0 * G * lambda_enclosed * std::log(r + 1e-30);
    }

    // ========================================================================
    // Report setup
    // ========================================================================
    std::cout << "\nParticle distribution:" << std::endl;
    std::cout << "  Total particles: " << particles.size() << std::endl;
    
    real x_min = particles[0].pos[0], x_max = particles[0].pos[0];
    real y_min = particles[0].pos[1], y_max = particles[0].pos[1];
    real z_min = particles[0].pos[2], z_max = particles[0].pos[2];
    real dens_min = particles[0].dens, dens_max = particles[0].dens;
    
    for (const auto& p : particles) {
        x_min = std::min(x_min, p.pos[0]);
        x_max = std::max(x_max, p.pos[0]);
        y_min = std::min(y_min, p.pos[1]);
        y_max = std::max(y_max, p.pos[1]);
        z_min = std::min(z_min, p.pos[2]);
        z_max = std::max(z_max, p.pos[2]);
        dens_min = std::min(dens_min, p.dens);
        dens_max = std::max(dens_max, p.dens);
    }
    
    std::cout << "  x range: [" << x_min << ", " << x_max << "]" << std::endl;
    std::cout << "  y range: [" << y_min << ", " << y_max << "]" << std::endl;
    std::cout << "  z range: [" << z_min << ", " << z_max << "]" << std::endl;
    std::cout << "  ρ range: [" << dens_min << " (edge), " << dens_max << " (center)]" << std::endl;

    // Estimate dynamical time
    const real t_dyn = R / c_s_center;
    std::cout << "\nTimescales:" << std::endl;
    std::cout << "  t_dyn = R/c_s = " << t_dyn << std::endl;

    std::cout << "\n=== 3D Lane-Emden Cylinder Setup Complete ===" << std::endl;
    std::cout << "NOTE: Use 'useKernelGravityCylinder3D: true' for correct gravity!" << std::endl;

    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());

#endif
}

} // namespace sph
