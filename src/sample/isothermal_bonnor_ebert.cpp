/**
 * @file isothermal_bonnor_ebert.cpp
 * @brief True isothermal Bonnor-Ebert sphere for self-gravitating hydrostatic equilibrium
 *
 * Solves the isothermal Lane-Emden equation:
 *   (1/xi^2) d/dxi(xi^2 dpsi/dxi) = exp(-psi)
 *
 * where:
 *   xi = r / r_0          (dimensionless radius)
 *   r_0 = c_s / sqrt(4*pi*G*rho_c)  (scale length)
 *   psi = -ln(rho/rho_c)  (dimensionless potential)
 *   rho = rho_c * exp(-psi)
 *
 * This gives TRUE hydrostatic equilibrium: dP/dr = -rho * g
 * where g = G*M(r)/r^2 is the gravitational acceleration.
 *
 * The critical Bonnor-Ebert sphere has xi_crit ~ 6.45
 * For xi < xi_crit: stable equilibrium
 * For xi > xi_crit: gravitationally unstable (collapse)
 */

#include "solver.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "sample/ghost_envelope.hpp"
#include "exception.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <random>

namespace sph
{

// Physical constants
static constexpr real k_B_cgs = 1.380649e-16;       // erg/K
static constexpr real m_proton_cgs = 1.6726219e-24; // g
static constexpr real mu_atomic = 1.27;             // Mean molecular weight for atomic HI + He
static constexpr real mu_molecular = 2.33;          // Mean molecular weight for H2 + He
static constexpr real Msun_cgs = 1.989e33;          // g
static constexpr real pc_cgs = 3.086e18;            // cm
static constexpr real kms_cgs = 1.0e5;              // cm/s

/**
 * @brief Solve isothermal Lane-Emden equation using RK4
 *
 * Equations:
 *   dpsi/dxi = z
 *   dz/dxi = exp(-psi) - 2*z/xi
 *
 * @param xi_max Maximum dimensionless radius
 * @param n_points Number of integration points
 * @param xi Output: dimensionless radius array
 * @param psi Output: dimensionless potential array
 * @param dpsi Output: dpsi/dxi array
 */
static void solve_isothermal_lane_emden(
    real xi_max, int n_points,
    std::vector<real>& xi,
    std::vector<real>& psi,
    std::vector<real>& dpsi)
{
    xi.resize(n_points);
    psi.resize(n_points);
    dpsi.resize(n_points);

    real dxi = xi_max / (n_points - 1);

    // Initial conditions at xi=0: psi=0, dpsi/dxi=0
    // Use series expansion near origin: psi ~ xi^2/6 for small xi
    xi[0] = 1e-6;  // Small but non-zero to avoid singularity
    psi[0] = xi[0] * xi[0] / 6.0;
    dpsi[0] = xi[0] / 3.0;

    // RK4 integration
    for (int i = 1; i < n_points; ++i) {
        real x = xi[i-1];
        real y = psi[i-1];
        real z = dpsi[i-1];

        // f1 = dy/dx = z
        // f2 = dz/dx = exp(-y) - 2*z/x
        auto f1 = [](real, real, real z_) { return z_; };
        auto f2 = [](real x_, real y_, real z_) {
            if (x_ < 1e-10) return 0.0;  // Avoid division by zero
            return std::exp(-y_) - 2.0 * z_ / x_;
        };

        // RK4 coefficients
        real k1_y = dxi * f1(x, y, z);
        real k1_z = dxi * f2(x, y, z);

        real k2_y = dxi * f1(x + 0.5*dxi, y + 0.5*k1_y, z + 0.5*k1_z);
        real k2_z = dxi * f2(x + 0.5*dxi, y + 0.5*k1_y, z + 0.5*k1_z);

        real k3_y = dxi * f1(x + 0.5*dxi, y + 0.5*k2_y, z + 0.5*k2_z);
        real k3_z = dxi * f2(x + 0.5*dxi, y + 0.5*k2_y, z + 0.5*k2_z);

        real k4_y = dxi * f1(x + dxi, y + k3_y, z + k3_z);
        real k4_z = dxi * f2(x + dxi, y + k3_y, z + k3_z);

        xi[i] = x + dxi;
        psi[i] = y + (k1_y + 2*k2_y + 2*k3_y + k4_y) / 6.0;
        dpsi[i] = z + (k1_z + 2*k2_z + 2*k3_z + k4_z) / 6.0;
    }
}

/**
 * @brief Create isothermal Bonnor-Ebert sphere IC
 *
 * Parameters from m_sample_parameters:
 *   - N: Grid resolution (N^3 target particles)
 *   - M_cloud: Cloud mass [M_sun]
 *   - T_cloud: Temperature [K]
 *   - xi_s: Dimensionless truncation radius (default 6.0, critical ~6.45)
 *   - P_ext: External pressure P/k_B [K cm^-3] (alternative to xi_s)
 */
void Solver::make_isothermal_bonnor_ebert()
{
#if DIM != 3
    THROW_ERROR("Isothermal Bonnor-Ebert requires DIM == 3");
#else
    std::cout << "==========================================================" << std::endl;
    std::cout << "  Isothermal Bonnor-Ebert Sphere (True Self-Gravitating)" << std::endl;
    std::cout << "==========================================================" << std::endl;

    // ========================================================================
    // READ PARAMETERS
    // ========================================================================

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real M_cloud = boost::any_cast<real>(m_sample_parameters["M_cloud"]);
    const real T_cloud = boost::any_cast<real>(m_sample_parameters["T_cloud"]);
    const real G_code = m_param->gravity.constant;
    const real gamma = m_param->physics.gamma;

    // Mean molecular weight (default: atomic HI + He = 1.27 for Option A)
    real mu = mu_atomic;  // Default to atomic gas
    if (m_sample_parameters.count("mu")) {
        mu = boost::any_cast<real>(m_sample_parameters["mu"]);
    }

    // Dimensionless truncation radius (critical BE sphere: xi_crit ~ 6.45)
    // For stability, use xi_s < 6.45
    real xi_s = 6.0;  // Default: slightly sub-critical
    if (m_sample_parameters.count("xi_s")) {
        xi_s = boost::any_cast<real>(m_sample_parameters["xi_s"]);
    }

    // Ghost envelope parameters
    const bool use_envelope = m_sample_parameters.count("useEnvelope") ?
                              boost::any_cast<bool>(m_sample_parameters["useEnvelope"]) : true;
    const int envelope_layers = m_sample_parameters.count("envelopeLayers") ?
                                boost::any_cast<int>(m_sample_parameters["envelopeLayers"]) : 4;

    std::cout << "Input parameters:" << std::endl;
    std::cout << "  M_cloud = " << M_cloud << " M_sun" << std::endl;
    std::cout << "  T_cloud = " << T_cloud << " K" << std::endl;
    std::cout << "  mu = " << mu << " (atomic=1.27, molecular=2.33)" << std::endl;
    std::cout << "  xi_s = " << xi_s << " (critical ~6.45)" << std::endl;
    std::cout << "  G = " << G_code << " (code units)" << std::endl;
    std::cout << "  gamma = " << gamma << std::endl;

    // ========================================================================
    // COMPUTE PHYSICAL SCALES
    // ========================================================================

    // Sound speed: c_s^2 = k_B T / (mu m_H)
    // In code units (km/s): c_s^2 [km/s]^2 = k_B T / (mu m_H) / (km/s)^2
    const real c_s_sq_cgs = k_B_cgs * T_cloud / (mu * m_proton_cgs);  // cm^2/s^2
    const real c_s_sq = c_s_sq_cgs / (kms_cgs * kms_cgs);  // (km/s)^2
    const real c_s = std::sqrt(c_s_sq);  // km/s

    std::cout << "\nSound speed:" << std::endl;
    std::cout << "  c_s = " << c_s << " km/s" << std::endl;
    std::cout << "  c_s^2 = " << c_s_sq << " (km/s)^2" << std::endl;

    // ========================================================================
    // SOLVE ISOTHERMAL LANE-EMDEN EQUATION
    // ========================================================================

    std::cout << "\nSolving isothermal Lane-Emden equation..." << std::endl;

    const int n_profile = 1000;
    std::vector<real> xi_arr, psi_arr, dpsi_arr;
    solve_isothermal_lane_emden(xi_s * 1.1, n_profile, xi_arr, psi_arr, dpsi_arr);

    // Find index closest to xi_s
    int i_s = 0;
    for (int i = 0; i < n_profile; ++i) {
        if (xi_arr[i] >= xi_s) {
            i_s = i;
            break;
        }
    }

    real psi_s = psi_arr[i_s];
    real dpsi_s = dpsi_arr[i_s];

    std::cout << "  At xi_s = " << xi_s << ":" << std::endl;
    std::cout << "    psi = " << psi_s << std::endl;
    std::cout << "    dpsi/dxi = " << dpsi_s << std::endl;
    std::cout << "    rho_edge/rho_center = " << std::exp(-psi_s) << std::endl;

    // ========================================================================
    // COMPUTE PHYSICAL PARAMETERS
    // ========================================================================
    // Two modes:
    // 1. Specify n_center (central number density) -> compute M_cloud
    // 2. Specify M_cloud -> compute n_center (results in low density for M=40, T=10K)
    //
    // For mode 1: Given rho_c, compute r_0 = c_s / sqrt(4*pi*G*rho_c)
    //             Then M = 4*pi*rho_c*r_0^3 * m_s
    // For mode 2: Given M, solve for rho_c (results in very diffuse cloud)

    real m_s = xi_s * xi_s * dpsi_s;  // Dimensionless enclosed mass at xi_s
    std::cout << "  Dimensionless mass m(xi_s) = " << m_s << std::endl;

    // Unit conversion for density
    const real density_code_to_cgs = Msun_cgs / std::pow(pc_cgs, 3);
    const real density_to_n = density_code_to_cgs / (mu * m_proton_cgs);

    real rho_c;
    real M_actual;

    // Check if n_center is specified (preferred mode for realistic clouds)
    if (m_sample_parameters.count("n_center")) {
        real n_center_input = boost::any_cast<real>(m_sample_parameters["n_center"]);
        rho_c = n_center_input / density_to_n;
        std::cout << "\nMode: Central density specified" << std::endl;
        std::cout << "  n_center (input) = " << n_center_input << " cm^-3" << std::endl;

        // Compute r_0 and M from rho_c
        real r_0_temp = c_s / std::sqrt(4.0 * M_PI * G_code * rho_c);
        M_actual = 4.0 * M_PI * rho_c * std::pow(r_0_temp, 3) * m_s;
        std::cout << "  Computed M_cloud = " << M_actual << " M_sun" << std::endl;
    } else {
        // Original mode: specify M_cloud, compute rho_c
        // From M = 4*pi*rho_c*r_0^3 * m_s where r_0 = c_s/sqrt(4*pi*G*rho_c)
        // => M = c_s^3 * m_s / (sqrt(4*pi) * G^(3/2) * sqrt(rho_c))
        // => sqrt(rho_c) = c_s^3 * m_s / (M * sqrt(4*pi) * G^(3/2))
        std::cout << "\nMode: Mass specified (WARNING: may result in very diffuse cloud)" << std::endl;
        real factor = std::sqrt(4.0 * M_PI) * std::pow(G_code, 1.5);
        real sqrt_rho_c = c_s * c_s * c_s * m_s / (M_cloud * factor);
        rho_c = sqrt_rho_c * sqrt_rho_c;
        M_actual = M_cloud;
    }

    // Scale length: r_0 = c_s / sqrt(4*pi*G*rho_c)
    real r_0 = c_s / std::sqrt(4.0 * M_PI * G_code * rho_c);

    // Physical cloud radius
    real R_cloud = xi_s * r_0;

    // Edge density
    real rho_edge = rho_c * std::exp(-psi_s);

    // Convert to number density (using already-defined density_to_n)
    real n_center = rho_c * density_to_n;
    real n_edge = rho_edge * density_to_n;

    // External pressure at truncation
    real P_ext = n_edge * T_cloud;  // K cm^-3

    std::cout << "\nDerived parameters:" << std::endl;
    std::cout << "  r_0 (scale length) = " << r_0 << " pc" << std::endl;
    std::cout << "  R_cloud = " << R_cloud << " pc" << std::endl;
    std::cout << "  rho_center = " << rho_c << " M_sun/pc^3" << std::endl;
    std::cout << "  n_center = " << n_center << " cm^-3" << std::endl;
    std::cout << "  rho_edge = " << rho_edge << " M_sun/pc^3" << std::endl;
    std::cout << "  n_edge = " << n_edge << " cm^-3" << std::endl;
    std::cout << "  P_ext/k_B = " << P_ext << " K cm^-3" << std::endl;

    // ========================================================================
    // BUILD CUMULATIVE MASS PROFILE FOR PARTICLE PLACEMENT
    // ========================================================================

    // Physical radius and density arrays
    std::vector<real> r_arr(n_profile);
    std::vector<real> rho_arr(n_profile);
    std::vector<real> M_enc_arr(n_profile);

    for (int i = 0; i < n_profile; ++i) {
        r_arr[i] = xi_arr[i] * r_0;
        rho_arr[i] = rho_c * std::exp(-psi_arr[i]);
        // Enclosed mass: M(r) = 4*pi*rho_c*r_0^3 * xi^2 * dpsi/dxi
        M_enc_arr[i] = 4.0 * M_PI * rho_c * r_0 * r_0 * r_0 * xi_arr[i] * xi_arr[i] * dpsi_arr[i];
    }

    // Normalize cumulative mass
    real M_total = M_enc_arr[i_s];
    std::vector<real> M_cumulative(n_profile);
    for (int i = 0; i < n_profile; ++i) {
        M_cumulative[i] = M_enc_arr[i] / M_total;
    }

    std::cout << "  Total enclosed mass = " << M_total << " M_sun" << std::endl;

    // ========================================================================
    // CREATE PARTICLES USING GLASS-LIKE RANDOM DISTRIBUTION
    // ========================================================================
    //
    // Random distribution eliminates lattice artifacts (shell banding).
    // Cumulative mass mapping ensures correct radial mass distribution.
    // Brief GLASS relaxation smooths local neighbor distribution.
    // ========================================================================

    const int N_target = N * N * N;
    std::cout << "\nCreating particles (Glass-like random + mass mapping):" << std::endl;
    std::cout << "  N_target = " << N_target << std::endl;

    // Generate random points uniformly in a unit sphere
    // Using rejection sampling for uniform distribution
    std::vector<vec_t> random_points;
    random_points.reserve(N_target);

    // Seed with fixed value for reproducibility
    std::mt19937 rng(42);
    std::uniform_real_distribution<real> uniform(-1.0, 1.0);

    while (static_cast<int>(random_points.size()) < N_target) {
        real x = uniform(rng);
        real y = uniform(rng);
        real z = uniform(rng);
        real r_sq = x*x + y*y + z*z;

        // Reject points outside unit sphere or at exact center
        if (r_sq < 1.0 && r_sq > 1e-10) {
            random_points.push_back(vec_t(x, y, z));
        }
    }

    std::cout << "  Random points: " << random_points.size() << std::endl;

    // Sort by radius for cumulative mass mapping
    std::vector<std::pair<real, int>> sorted_points;
    for (size_t i = 0; i < random_points.size(); ++i) {
        real r = std::abs(random_points[i]);
        sorted_points.push_back({r, static_cast<int>(i)});
    }
    std::sort(sorted_points.begin(), sorted_points.end());

    // Interpolation functions
    auto interpolate_radius = [&](real mass_frac) -> real {
        for (int i = 1; i <= i_s; ++i) {
            if (mass_frac <= M_cumulative[i]) {
                real f = (mass_frac - M_cumulative[i-1]) / (M_cumulative[i] - M_cumulative[i-1]);
                return r_arr[i-1] + f * (r_arr[i] - r_arr[i-1]);
            }
        }
        return R_cloud;
    };

    auto interpolate_density = [&](real r) -> real {
        for (int i = 1; i <= i_s; ++i) {
            if (r <= r_arr[i]) {
                real f = (r - r_arr[i-1]) / (r_arr[i] - r_arr[i-1]);
                return rho_arr[i-1] + f * (rho_arr[i] - rho_arr[i-1]);
            }
        }
        return rho_edge;
    };

    // Map random points to BE density profile using cumulative mass
    std::vector<vec_t> mapped_positions;
    std::vector<real> local_densities;
    int N_points = static_cast<int>(sorted_points.size());

    for (int i = 0; i < N_points; ++i) {
        real mass_frac = (i + 0.5) / N_points;
        real r_be = interpolate_radius(mass_frac);

        int idx = sorted_points[i].second;
        real r_orig = sorted_points[i].first;

        if (r_orig > 1e-10) {
            real scale = r_be / r_orig;
            vec_t new_pos = random_points[idx] * scale;
            mapped_positions.push_back(new_pos);
            local_densities.push_back(interpolate_density(r_be));
        }
    }

    std::cout << "  Mapped " << mapped_positions.size() << " particles to BE profile" << std::endl;

    // ========================================================================
    // CREATE SPH PARTICLES
    // ========================================================================

    const real particle_mass = M_actual / mapped_positions.size();
    std::cout << "  Particle mass = " << particle_mass << " M_sun" << std::endl;

    // Internal energy: u = c_s^2 / (gamma - 1) for isothermal
    const real u_cloud = c_s_sq / (gamma - 1.0);
    std::cout << "  u_cloud = " << u_cloud << " (code units)" << std::endl;

    auto& particles = m_sim->get_particles();
    particles.clear();
    particles.reserve(mapped_positions.size() + 5000);  // Extra for envelope

    int particle_id = 0;
    const int N_neighbor = m_param->physics.neighbor_number;
    constexpr real A = 4.0 * M_PI / 3.0;

    for (size_t idx = 0; idx < mapped_positions.size(); ++idx) {
        const auto& pos = mapped_positions[idx];
        real rho_local = local_densities[idx];

        // Smoothing length from local density
        real h_est = std::pow(N_neighbor * particle_mass / (rho_local * A), 1.0 / 3.0);

        SPHParticle p;
        p.pos = pos;
        p.vel = 0.0;
        p.mass = particle_mass;
        p.dens = rho_local;
        p.ene = u_cloud;
        p.pres = (gamma - 1.0) * rho_local * u_cloud;
        p.sml = h_est;
        p.id = particle_id++;
        p.is_ghost = false;

        particles.push_back(p);
    }

    std::cout << "  Created " << particles.size() << " cloud particles" << std::endl;

    // ========================================================================
    // CREATE GHOST ENVELOPE
    // ========================================================================

    if (use_envelope) {
        GhostEnvelopeConfig env_config;
        env_config.R_cloud = R_cloud;
        env_config.rho_edge = rho_edge;
        env_config.u_envelope = u_cloud;
        env_config.particle_mass = particle_mass;
        env_config.N_neighbor = m_param->physics.neighbor_number;
        env_config.num_layers = envelope_layers;

        auto envelope_particles = GhostEnvelopeGenerator::generate(env_config);

        for (auto& p : envelope_particles) {
            p.id = particle_id++;
        }

        particles.insert(particles.end(), envelope_particles.begin(), envelope_particles.end());
        GhostEnvelopeGenerator::print_summary(env_config, envelope_particles.size());
    }

    // ========================================================================
    // STORE PARAMETERS FOR OUTPUT/ANALYSIS
    // ========================================================================

    m_sample_parameters["R_cloud"] = R_cloud;
    m_sample_parameters["R_cloud_code"] = R_cloud;  // For relaxation compatibility
    m_sample_parameters["r_0"] = r_0;
    m_sample_parameters["r_c_code"] = r_0;  // For relaxation: use r_0 as core radius
    m_sample_parameters["rho_center"] = rho_c;
    m_sample_parameters["rho_center_code"] = rho_c;  // For relaxation compatibility
    m_sample_parameters["rho_edge"] = rho_edge;
    m_sample_parameters["c_s"] = c_s;
    m_sample_parameters["P_ext"] = P_ext;
    m_sample_parameters["T_cloud"] = T_cloud;  // For relaxation
    m_sample_parameters["density_to_n"] = density_to_n;
    m_sample_parameters["mu"] = mu;
    m_sample_parameters["xi_s"] = xi_s;

    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());

    std::cout << "\n=== Isothermal Bonnor-Ebert IC Complete ===" << std::endl;
    std::cout << "Total particles: " << particles.size() << std::endl;
    std::cout << "This configuration is in TRUE hydrostatic equilibrium with self-gravity." << std::endl;
#endif  // DIM == 3
}

} // namespace sph
