/**
 * @file uniform_cloud.cpp
 * @brief Super-critical Bonnor-Ebert sphere for IMBH-cloud tidal interaction
 *
 * Creates a super-critical isothermal Bonnor-Ebert sphere for studying
 * IMBH-cloud interactions where tidal forces create density perturbations
 * and drive gravitational collapse.
 *
 * Key features:
 *   - Bonnor-Ebert density profile (centrally concentrated)
 *   - M = 40 M_sun at T = 15 K gives R ~ 3 pc
 *   - Super-critical: M/M_crit ~ 20-40 (gravitationally unstable)
 *   - ZERO initial velocity (no turbulence injection)
 *   - IMBH tidal forces create density gradients naturally
 *   - Compression ratio: R: 3pc -> 0.3pc gives density x1000
 *
 * Physics:
 *   - Solves isothermal Lane-Emden equation for BE profile
 *   - T = 15 K (KI2000 equilibrium at high density)
 *   - Thermal support << Gravitational energy (super-Jeans)
 *   - IMBH present from start at distance ~10-15 pc
 *
 * Reference:
 *   - Oka et al. (2017) - HVCC CO-0.40-0.22 candidate IMBH
 *   - Koyama & Inutsuka (2000) - ISM thermal equilibrium
 *   - Bonnor (1956), Ebert (1955) - Isothermal sphere equilibrium
 */

#include "solver.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "exception.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace sph {

// ============================================================================
// Physical constants (CGS)
// ============================================================================
namespace {
    constexpr real k_B_cgs = 1.380649e-16;       // erg/K
    constexpr real m_proton_cgs = 1.6726219e-24; // g
    constexpr real mu_mol = 2.33;                // Mean molecular weight (molecular H2)
    constexpr real m_n_cgs = mu_mol * m_proton_cgs;  // Mean particle mass

    // Unit conversions (code units: M[M_sun], L[pc], V[km/s])
    constexpr real Msun_cgs = 1.989e33;          // g
    constexpr real pc_cgs = 3.086e18;            // cm
    constexpr real kms_cgs = 1.0e5;              // cm/s

    // Code density [M_sun/pc^3] to CGS [g/cm^3]
    constexpr real density_code_to_cgs = Msun_cgs / (pc_cgs * pc_cgs * pc_cgs);

    // Code density to number density [cm^-3]
    constexpr real density_to_n_factor = density_code_to_cgs / m_n_cgs;

    /**
     * @brief Isothermal sound speed
     * @param T Temperature [K]
     * @param mu Mean molecular weight
     * @return Sound speed [km/s]
     */
    real isothermal_sound_speed(real T, real mu)
    {
        real c_s_cgs = std::sqrt(k_B_cgs * T / (mu * m_proton_cgs));
        return c_s_cgs / kms_cgs;  // km/s
    }

    /**
     * @brief Convert code density to number density
     * @param rho Density [M_sun/pc^3]
     * @return Number density [cm^-3]
     */
    real code_density_to_n(real rho)
    {
        return rho * density_to_n_factor;  // cm^-3
    }

    /**
     * @brief Solve isothermal Lane-Emden equation using RK4
     *
     * Equations:
     *   dpsi/dxi = z
     *   dz/dxi = exp(-psi) - 2*z/xi
     *
     * where:
     *   xi = r / r_0          (dimensionless radius)
     *   r_0 = c_s / sqrt(4*pi*G*rho_c)  (scale length)
     *   psi = -ln(rho/rho_c)  (dimensionless potential)
     */
    void solve_isothermal_lane_emden(
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
                if (x_ < 1e-10) return real(0.0);
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

#if DIM == 3
    /**
     * @brief Generate HCP (hexagonal close-packed) lattice points in a sphere
     *
     * Creates a structured lattice that fills a unit sphere with approximately
     * N_target particles. Returns positions normalized to unit sphere.
     */
    std::vector<vec_t> generate_hcp_sphere(int N_target)
    {
        std::vector<vec_t> points;
        points.reserve(N_target * 2);

        // HCP lattice parameters
        real V_sphere = 4.0 * M_PI / 3.0;
        real a = std::pow(V_sphere / (0.74 * N_target), 1.0/3.0);

        // HCP layer spacing
        real dz = a * std::sqrt(2.0/3.0);
        real dy = a * std::sqrt(3.0) / 2.0;

        // Generate lattice layers
        int nz = static_cast<int>(2.0 / dz) + 2;
        int ny = static_cast<int>(2.0 / dy) + 2;
        int nx = static_cast<int>(2.0 / a) + 2;

        for (int iz = -nz; iz <= nz; ++iz) {
            real z = iz * dz;
            real y_offset = (iz % 2 == 0) ? 0.0 : dy / 3.0;
            real x_offset = (iz % 2 == 0) ? 0.0 : a / 2.0;

            for (int iy = -ny; iy <= ny; ++iy) {
                real y = iy * dy + y_offset;
                real row_x_offset = (iy % 2 == 0) ? 0.0 : a / 2.0;

                for (int ix = -nx; ix <= nx; ++ix) {
                    real x = ix * a + x_offset + row_x_offset;

                    vec_t pos(x, y, z);
                    real r = std::abs(pos);

                    // Keep only points inside unit sphere
                    if (r <= 1.0 && r > 1.0e-6) {
                        points.push_back(pos);
                    }
                }
            }
        }

        std::cout << "  HCP lattice: generated " << points.size()
                  << " points (target " << N_target << ")" << std::endl;
        return points;
    }
#endif // DIM == 3

} // anonymous namespace

// ============================================================================
// Configuration structure
// ============================================================================
struct BECloudConfig {
    // Cloud parameters
    real M_cloud = 40.0;         // Cloud mass [M_sun]
    real T_cloud = 15.0;         // Temperature [K] (KI2000 equilibrium)
    real xi_s = 10.0;            // Dimensionless truncation radius (super-critical)
    int N_particles = 10000;     // Number of particles

    // IMBH parameters
    real M_BH = 1.0e5;           // Black hole mass [M_sun]
    real d_BH = 10.0;            // Initial BH distance [pc] (along +x axis)
    real epsilon_BH = 0.01;      // BH softening length [pc]

    // Envelope parameters
    bool use_envelope = true;
    int envelope_layers = 4;
};

// ============================================================================
// Main IC generator: make_uniform_cloud()
// Now creates super-critical Bonnor-Ebert sphere
// ============================================================================

void Solver::make_uniform_cloud()
{
#if DIM != 3
    THROW_ERROR("BE sphere cloud requires DIM == 3");
#else
    // ========================================================================
    // READ PARAMETERS
    // ========================================================================

    BECloudConfig config;

    // Override defaults from m_sample_parameters
    if (m_sample_parameters.count("N")) {
        int N = boost::any_cast<int>(m_sample_parameters["N"]);
        config.N_particles = N * N * N;  // N is grid resolution
    }
    if (m_sample_parameters.count("M_cloud"))
        config.M_cloud = boost::any_cast<real>(m_sample_parameters["M_cloud"]);
    if (m_sample_parameters.count("T_cloud"))
        config.T_cloud = boost::any_cast<real>(m_sample_parameters["T_cloud"]);
    if (m_sample_parameters.count("xi_s"))
        config.xi_s = boost::any_cast<real>(m_sample_parameters["xi_s"]);
    if (m_sample_parameters.count("M_BH"))
        config.M_BH = boost::any_cast<real>(m_sample_parameters["M_BH"]);
    if (m_sample_parameters.count("d_BH"))
        config.d_BH = boost::any_cast<real>(m_sample_parameters["d_BH"]);
    if (m_sample_parameters.count("epsilon_BH"))
        config.epsilon_BH = boost::any_cast<real>(m_sample_parameters["epsilon_BH"]);
    if (m_sample_parameters.count("useEnvelope"))
        config.use_envelope = boost::any_cast<bool>(m_sample_parameters["useEnvelope"]);
    if (m_sample_parameters.count("envelopeLayers"))
        config.envelope_layers = boost::any_cast<int>(m_sample_parameters["envelopeLayers"]);

    const real gamma = m_param->physics.gamma;
    const real G = m_param->gravity.constant;

    std::cout << "==========================================================" << std::endl;
    std::cout << "  Super-Critical Bonnor-Ebert Sphere for IMBH Interaction" << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << std::endl;
    std::cout << "Physical model:" << std::endl;
    std::cout << "  - BONNOR-EBERT density profile (centrally concentrated)" << std::endl;
    std::cout << "  - T = " << config.T_cloud << " K (KI2000 equilibrium)" << std::endl;
    std::cout << "  - ZERO initial velocity (no turbulence injection)" << std::endl;
    std::cout << "  - IMBH tidal forces create density perturbations" << std::endl;
    std::cout << "  - Super-critical: M/M_crit ~ 20-40" << std::endl;
    std::cout << std::endl;

    // ========================================================================
    // SOLVE ISOTHERMAL LANE-EMDEN EQUATION
    // ========================================================================

    // Isothermal sound speed
    const real c_s = isothermal_sound_speed(config.T_cloud, mu_mol);
    const real c_s_sq = c_s * c_s;

    std::cout << "Sound speed:" << std::endl;
    std::cout << "  c_s = " << c_s << " km/s (at T = " << config.T_cloud << " K)" << std::endl;
    std::cout << std::endl;

    std::cout << "Solving isothermal Lane-Emden equation..." << std::endl;

    const int n_profile = 1000;
    std::vector<real> xi_arr, psi_arr, dpsi_arr;
    solve_isothermal_lane_emden(config.xi_s * 1.1, n_profile, xi_arr, psi_arr, dpsi_arr);

    // Find index closest to xi_s
    int i_s = 0;
    for (int i = 0; i < n_profile; ++i) {
        if (xi_arr[i] >= config.xi_s) {
            i_s = i;
            break;
        }
    }

    real psi_s = psi_arr[i_s];
    real dpsi_s = dpsi_arr[i_s];
    real m_s = config.xi_s * config.xi_s * dpsi_s;  // Dimensionless enclosed mass

    std::cout << "  At xi_s = " << config.xi_s << " (critical ~6.45):" << std::endl;
    std::cout << "    psi = " << psi_s << std::endl;
    std::cout << "    rho_edge/rho_center = " << std::exp(-psi_s) << std::endl;
    std::cout << "    Dimensionless mass m(xi_s) = " << m_s << std::endl;
    std::cout << std::endl;

    // ========================================================================
    // COMPUTE PHYSICAL PARAMETERS FROM MASS
    // ========================================================================

    // From M = 4*pi*rho_c*r_0^3 * m_s where r_0 = c_s/sqrt(4*pi*G*rho_c)
    // => M = c_s^3 * m_s / (sqrt(4*pi) * G^(3/2) * sqrt(rho_c))
    // => sqrt(rho_c) = c_s^3 * m_s / (M * sqrt(4*pi) * G^(3/2))

    real factor = std::sqrt(4.0 * M_PI) * std::pow(G, 1.5);
    real sqrt_rho_c = c_s * c_s * c_s * m_s / (config.M_cloud * factor);
    real rho_c = sqrt_rho_c * sqrt_rho_c;

    // Scale length: r_0 = c_s / sqrt(4*pi*G*rho_c)
    real r_0 = c_s / std::sqrt(4.0 * M_PI * G * rho_c);

    // Physical cloud radius
    real R_cloud = config.xi_s * r_0;

    // Edge density
    real rho_edge = rho_c * std::exp(-psi_s);

    // Convert to number density
    real n_center = code_density_to_n(rho_c);
    real n_edge = code_density_to_n(rho_edge);

    // Average density for Jeans mass calculation
    real V_cloud = (4.0 / 3.0) * M_PI * R_cloud * R_cloud * R_cloud;
    real rho_avg = config.M_cloud / V_cloud;

    // Jeans mass at central density and temperature
    real M_J_center = std::pow(M_PI, 2.5) * std::pow(c_s, 3) /
                      (6.0 * std::pow(G, 1.5) * std::sqrt(rho_c));

    // Critical BE mass at this temperature
    // M_crit = 1.18 * c_s^4 / (G^(3/2) * sqrt(P_ext))
    // where P_ext = rho_edge * c_s^2
    real P_ext = rho_edge * c_s_sq;
    real M_crit = 1.18 * std::pow(c_s, 4) / (std::pow(G, 1.5) * std::sqrt(P_ext));

    // Free-fall time at central density
    real t_ff = std::sqrt(3.0 * M_PI / (32.0 * G * rho_c));

    std::cout << "Derived parameters:" << std::endl;
    std::cout << "  r_0 (scale length) = " << r_0 << " pc" << std::endl;
    std::cout << "  R_cloud = " << R_cloud << " pc" << std::endl;
    std::cout << "  rho_center = " << rho_c << " M_sun/pc^3" << std::endl;
    std::cout << "  n_center = " << n_center << " cm^-3" << std::endl;
    std::cout << "  rho_edge = " << rho_edge << " M_sun/pc^3" << std::endl;
    std::cout << "  n_edge = " << n_edge << " cm^-3" << std::endl;
    std::cout << std::endl;
    std::cout << "Stability analysis:" << std::endl;
    std::cout << "  M_cloud = " << config.M_cloud << " M_sun" << std::endl;
    std::cout << "  M_crit (BE) ~ " << M_crit << " M_sun" << std::endl;
    std::cout << "  M/M_crit = " << config.M_cloud / M_crit << " (SUPER-CRITICAL!)" << std::endl;
    std::cout << "  M_J (center) = " << M_J_center << " M_sun" << std::endl;
    std::cout << "  t_ff (center) = " << t_ff << " Myr" << std::endl;
    std::cout << std::endl;

    // ========================================================================
    // BUILD CUMULATIVE MASS PROFILE FOR PARTICLE PLACEMENT
    // ========================================================================

    std::vector<real> r_arr(n_profile);
    std::vector<real> rho_arr(n_profile);
    std::vector<real> M_enc_arr(n_profile);

    for (int i = 0; i < n_profile; ++i) {
        r_arr[i] = xi_arr[i] * r_0;
        rho_arr[i] = rho_c * std::exp(-psi_arr[i]);
        // Enclosed mass: M(r) = 4*pi*rho_c*r_0^3 * xi^2 * dpsi/dxi
        M_enc_arr[i] = 4.0 * M_PI * rho_c * r_0 * r_0 * r_0 * xi_arr[i] * xi_arr[i] * dpsi_arr[i];
    }

    real M_total = M_enc_arr[i_s];
    std::vector<real> M_cumulative(n_profile);
    for (int i = 0; i < n_profile; ++i) {
        M_cumulative[i] = M_enc_arr[i] / M_total;
    }

    // ========================================================================
    // CREATE PARTICLES USING HCP LATTICE MAPPED TO BE PROFILE
    // ========================================================================

    std::cout << "Creating particles:" << std::endl;
    std::cout << "  N_target = " << config.N_particles << std::endl;

    // Generate HCP lattice in unit sphere
    std::vector<vec_t> lattice_points = generate_hcp_sphere(config.N_particles);

    // Sort by radius for BE mapping
    std::vector<std::pair<real, int>> sorted_points;
    for (size_t i = 0; i < lattice_points.size(); ++i) {
        real r = std::abs(lattice_points[i]);
        sorted_points.push_back({r, static_cast<int>(i)});
    }
    std::sort(sorted_points.begin(), sorted_points.end());

    // Interpolation function for radius from cumulative mass
    auto interpolate_radius = [&](real mass_frac) -> real {
        for (int i = 1; i < i_s; ++i) {
            if (mass_frac <= M_cumulative[i]) {
                real f = (mass_frac - M_cumulative[i-1]) / (M_cumulative[i] - M_cumulative[i-1]);
                return r_arr[i-1] + f * (r_arr[i] - r_arr[i-1]);
            }
        }
        return R_cloud;
    };

    // Interpolation function for density
    auto interpolate_density = [&](real r) -> real {
        for (int i = 1; i < i_s; ++i) {
            if (r <= r_arr[i]) {
                real f = (r - r_arr[i-1]) / (r_arr[i] - r_arr[i-1]);
                return rho_arr[i-1] + f * (rho_arr[i] - rho_arr[i-1]);
            }
        }
        return rho_edge;
    };

    // Map uniform lattice to BE profile
    std::vector<vec_t> mapped_positions;
    int N_lattice = static_cast<int>(sorted_points.size());

    for (int i = 0; i < N_lattice; ++i) {
        real mass_frac = (i + 0.5) / N_lattice;
        real r_be = interpolate_radius(mass_frac);

        int idx = sorted_points[i].second;
        real r_orig = sorted_points[i].first;

        if (r_orig > 1e-10) {
            real scale = r_be / r_orig;
            mapped_positions.push_back(lattice_points[idx] * scale);
        }
    }

    std::cout << "  Mapped " << mapped_positions.size() << " particles to BE profile" << std::endl;

    // ========================================================================
    // CREATE SPH PARTICLES
    // ========================================================================

    const int N_actual = static_cast<int>(mapped_positions.size());
    const real particle_mass = config.M_cloud / N_actual;

    std::cout << "  Particle mass = " << particle_mass << " M_sun" << std::endl;

    // Internal energy: u = c_s^2 / (gamma - 1)
    const real gamma_eff = std::max(gamma - 1.0, 1.0e-4);
    const real u_cloud = c_s_sq / gamma_eff;

    const int N_neighbor = m_param->physics.neighbor_number;
    constexpr real A_vol = 4.0 * M_PI / 3.0;

    std::vector<SPHParticle> particles;
    particles.reserve(N_actual + 5000);

    int particle_id = 0;
    for (const auto& pos : mapped_positions) {
        real r = std::abs(pos);
        real rho_local = interpolate_density(r);

        // Smoothing length from local density
        real h_local = std::pow(N_neighbor * particle_mass / (rho_local * A_vol), 1.0 / 3.0);

        SPHParticle p;
        p.pos = pos;
        p.vel = vec_t{0.0, 0.0, 0.0};  // ZERO velocity!
        p.mass = particle_mass;
        p.dens = rho_local;
        p.ene = u_cloud;
        p.pres = rho_local * c_s_sq;  // Isothermal: P = rho * c_s^2
        p.sml = h_local;
        p.is_ghost = false;
        p.id = particle_id++;

        particles.push_back(p);
    }

    std::cout << "  Created " << particles.size() << " cloud particles" << std::endl;
    std::cout << "  u_cloud = " << u_cloud << " (code units)" << std::endl;
    std::cout << std::endl;

    // ========================================================================
    // CREATE GHOST ENVELOPE
    // ========================================================================

    if (config.use_envelope) {
        std::cout << "Creating ghost envelope:" << std::endl;
        std::cout << "  Envelope layers: " << config.envelope_layers << std::endl;
        std::cout << "  Envelope density: " << rho_edge << " M_sun/pc^3" << std::endl;

        real h_envelope = std::pow(N_neighbor * particle_mass / (rho_edge * A_vol), 1.0 / 3.0);
        real layer_spacing = 2.0 * h_envelope;

        int envelope_count = 0;
        for (int layer = 1; layer <= config.envelope_layers; ++layer) {
            real R_layer = R_cloud + layer * layer_spacing;
            real shell_area = 4.0 * M_PI * R_layer * R_layer;
            int N_shell = static_cast<int>(shell_area / (layer_spacing * layer_spacing));

            for (int i = 0; i < N_shell; ++i) {
                // Fibonacci sphere distribution
                real golden_ratio = (1.0 + std::sqrt(5.0)) / 2.0;
                real theta_angle = 2.0 * M_PI * i / golden_ratio;
                real z = 1.0 - 2.0 * (static_cast<real>(i) + 0.5) / N_shell;
                real radius_xy = std::sqrt(1.0 - z * z);

                SPHParticle p;
                p.pos = vec_t{
                    R_layer * radius_xy * std::cos(theta_angle),
                    R_layer * radius_xy * std::sin(theta_angle),
                    R_layer * z
                };
                p.vel = vec_t{0.0, 0.0, 0.0};
                p.mass = particle_mass;
                p.dens = rho_edge;
                p.ene = u_cloud;
                p.pres = rho_edge * c_s_sq;
                p.sml = h_envelope;
                p.is_ghost = true;
                p.id = particle_id++;

                particles.push_back(p);
                envelope_count++;
            }
        }

        std::cout << "  Created " << envelope_count << " envelope particles" << std::endl;
        std::cout << std::endl;
    }

    // ========================================================================
    // STORE PARAMETERS FOR EXTERNAL BH AND OUTPUT
    // ========================================================================

    m_sample_parameters["BH_position_x"] = config.d_BH;
    m_sample_parameters["BH_position_y"] = real(0.0);
    m_sample_parameters["BH_position_z"] = real(0.0);
    m_sample_parameters["BH_mass"] = config.M_BH;
    m_sample_parameters["BH_softening"] = config.epsilon_BH;
    m_sample_parameters["T_cloud_K"] = config.T_cloud;
    m_sample_parameters["c_s_kms"] = c_s;
    m_sample_parameters["R_cloud"] = R_cloud;
    m_sample_parameters["M_cloud"] = config.M_cloud;
    m_sample_parameters["rho_center"] = rho_c;
    m_sample_parameters["rho_edge"] = rho_edge;
    m_sample_parameters["n_center"] = n_center;
    m_sample_parameters["n_edge"] = n_edge;
    m_sample_parameters["M_crit"] = M_crit;
    m_sample_parameters["t_ff_Myr"] = t_ff;

    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());

    // ========================================================================
    // PRINT SUMMARY
    // ========================================================================

    real r_tidal = R_cloud * std::pow(config.M_cloud / (3.0 * config.M_BH), 1.0 / 3.0);
    real a_tidal = G * config.M_BH * R_cloud / (config.d_BH * config.d_BH * config.d_BH);

    std::cout << "==========================================================" << std::endl;
    std::cout << "  Super-Critical BE Sphere IC Complete" << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << "Cloud (Bonnor-Ebert profile):" << std::endl;
    std::cout << "  R = " << R_cloud << " pc" << std::endl;
    std::cout << "  M = " << config.M_cloud << " M_sun" << std::endl;
    std::cout << "  T = " << config.T_cloud << " K" << std::endl;
    std::cout << "  n_center = " << n_center << " cm^-3" << std::endl;
    std::cout << "  n_edge = " << n_edge << " cm^-3" << std::endl;
    std::cout << "  M/M_crit = " << config.M_cloud / M_crit << " (SUPER-CRITICAL)" << std::endl;
    std::cout << std::endl;
    std::cout << "Particles: " << particles.size() << std::endl;
    std::cout << std::endl;
    std::cout << "IMBH:" << std::endl;
    std::cout << "  M_BH = " << config.M_BH << " M_sun" << std::endl;
    std::cout << "  Position: (" << config.d_BH << ", 0, 0) pc" << std::endl;
    std::cout << "  r_tidal ~ " << r_tidal << " pc" << std::endl;
    std::cout << "  a_tidal ~ " << a_tidal << " km^2/s^2/pc" << std::endl;
    std::cout << std::endl;
    std::cout << "Timescales:" << std::endl;
    std::cout << "  t_ff (center) = " << t_ff << " Myr" << std::endl;
    std::cout << std::endl;
    std::cout << "Expected behavior:" << std::endl;
    std::cout << "  - Cloud collapses (super-critical BE sphere)" << std::endl;
    std::cout << "  - IMBH tidal forces create density gradients" << std::endl;
    std::cout << "  - Compression: R " << R_cloud << " pc -> ~0.3 pc (x" << (int)(R_cloud/0.3) << ")" << std::endl;
    std::cout << "  - Density increase: ~x" << (int)(std::pow(R_cloud/0.3, 3)) << " (to ~10^6 cm^-3)" << std::endl;
    std::cout << "==========================================================" << std::endl;
#endif // DIM == 3
}

} // namespace sph
