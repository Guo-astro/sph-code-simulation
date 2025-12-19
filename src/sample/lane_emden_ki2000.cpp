/**
 * @file lane_emden_ki2000.cpp
 * @brief Lane-Emden density structure with K&I 2000 thermal equilibrium temperatures
 *
 * This creates a "hybrid" IC that is:
 *   - Lane-Emden n=3/2 density profile (self-gravitating sphere)
 *   - K&I 2000 equilibrium temperature at each density: T = T_eq(n)
 *   - NOT in hydrostatic equilibrium (pressure gradient ≠ gravity)
 *
 * Physics rationale:
 *   - Lane-Emden gives proper mass distribution for self-gravitating cloud
 *   - K&I temperatures ensure thermal equilibrium with cooling/heating
 *   - The cloud will start to adjust when gravity is enabled
 *   - For tidal disruption: t_tidal << t_collapse, so dynamics will be captured
 *
 * Use case:
 *   - IMBH-cloud tidal disruption simulations
 *   - Large (~1 pc) molecular cloud ICs
 *   - Simulations with K&I cooling enabled from the start
 */

#include "solver.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "thermal/koyama_inutsuka_cooling.hpp"
#include "exception.hpp"
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <sstream>
#include <iostream>

namespace sph
{

// Physical constants
static constexpr real k_B_cgs = 1.380649e-16;       // erg/K
static constexpr real m_proton_cgs = 1.6726219e-24; // g
static constexpr real mu_mol = 2.33;                // Mean molecular weight (molecular gas)
static constexpr real m_n_cgs = mu_mol * m_proton_cgs;

// Unit conversions
static constexpr real Msun_cgs = 1.989e33;          // g
static constexpr real pc_cgs = 3.086e18;            // cm
static constexpr real kms_cgs = 1.0e5;              // cm/s

/**
 * @brief Create Lane-Emden structure with K&I 2000 temperatures
 *
 * Parameters:
 *   - N: Number of particles (N³ in glass, or total for shell distribution)
 *   - R_cloud_pc: Cloud radius [pc]
 *   - M_cloud_Msun: Cloud mass [M_☉]
 *   - N_H_cm2: Column density for K&I tables [cm⁻²]
 *   - use_glass: If true, use glass-making method; if false, use shell distribution
 */
void Solver::make_lane_emden_ki2000()
{
#if DIM != 3
    THROW_ERROR("Lane-Emden K&I 2000 test requires DIM == 3");
#endif

    // ========================================================================
    // READ PARAMETERS
    // ========================================================================

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real R_cloud = boost::any_cast<real>(m_sample_parameters.count("R_cloud_pc") ?
                            m_sample_parameters["R_cloud_pc"] : boost::any(1.0));
    const real M_cloud = boost::any_cast<real>(m_sample_parameters.count("M_cloud_Msun") ?
                            m_sample_parameters["M_cloud_Msun"] : boost::any(1000.0));
    const real N_H_cm2 = boost::any_cast<real>(m_sample_parameters.count("N_H_cm2") ?
                            m_sample_parameters["N_H_cm2"] : boost::any(1.0e19));
    const bool use_glass = boost::any_cast<bool>(m_sample_parameters.count("use_glass") ?
                            m_sample_parameters["use_glass"] : boost::any(true));

    const real gamma = m_param->physics.gamma;
    const real G = m_param->gravity.constant;

    std::cout << "==========================================================" << std::endl;
    std::cout << "  Lane-Emden Structure + K&I 2000 Temperatures" << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << "Cloud parameters:" << std::endl;
    std::cout << "  R_cloud = " << R_cloud << " pc" << std::endl;
    std::cout << "  M_cloud = " << M_cloud << " M_☉" << std::endl;
    std::cout << "  N_H = " << N_H_cm2 << " cm⁻²" << std::endl;
    std::cout << "  γ = " << gamma << std::endl;
    std::cout << "  G = " << G << std::endl;
    std::cout << std::endl;

    // ========================================================================
    // UNIT CONVERSIONS
    // ========================================================================

    // Code units: L[pc], M[M_☉], V[km/s]
    const real density_code_to_cgs = Msun_cgs / (pc_cgs * pc_cgs * pc_cgs);
    const real density_to_n = density_code_to_cgs / m_n_cgs;  // code → n_H [cm⁻³]

    std::cout << "Unit conversions:" << std::endl;
    std::cout << "  density_code → n_H: " << density_to_n << " cm⁻³" << std::endl;

    // ========================================================================
    // LOAD LANE-EMDEN SOLUTION
    // ========================================================================

    std::string filename = "data/lane_emden/n1.5_3d.dat";
    std::ifstream infile(filename);

    if (!infile.is_open()) {
        // Try alternative path
        filename = "lane_emden/data/numerical_solutions/3d/n1.5.dat";
        infile.open(filename);
    }

    if (!infile.is_open()) {
        THROW_ERROR("Cannot open Lane-Emden solution file. Tried:\n"
                    "  - data/lane_emden/n1.5_3d.dat\n"
                    "  - lane_emden/data/numerical_solutions/3d/n1.5.dat");
    }

    real xi_1 = 3.65375;  // Default for n=1.5
    real dtheta_1 = -0.20330;
    std::vector<real> xi_array, theta_array;

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') {
            // Parse header for xi_1 and dtheta_1
            if (line.find("xi_1") != std::string::npos) {
                sscanf(line.c_str(), "# xi_1 = %lf", &xi_1);
            } else if (line.find("dtheta_1") != std::string::npos) {
                sscanf(line.c_str(), "# dtheta_1 = %lf", &dtheta_1);
            }
            continue;
        }

        real xi, theta, dtheta;
        std::istringstream iss(line);
        if (iss >> xi >> theta >> dtheta) {
            xi_array.push_back(xi);
            theta_array.push_back(theta);
        }
    }
    infile.close();

    if (xi_array.empty()) {
        THROW_ERROR("Failed to read Lane-Emden data");
    }

    std::cout << "Lane-Emden solution loaded:" << std::endl;
    std::cout << "  ξ₁ = " << xi_1 << std::endl;
    std::cout << "  |dθ/dξ|_{ξ₁} = " << std::abs(dtheta_1) << std::endl;
    std::cout << "  Data points: " << xi_array.size() << std::endl;

    // ========================================================================
    // COMPUTE PHYSICAL SCALING
    // ========================================================================

    // Lane-Emden: r = α × ξ, where R_cloud = α × ξ₁
    const real alpha = R_cloud / xi_1;

    // Central density from mass normalization:
    // M = 4π α³ ρ_c × ξ₁² |dθ/dξ|_{ξ₁}
    const real rho_center = M_cloud / (4.0 * M_PI * alpha * alpha * alpha *
                                        xi_1 * xi_1 * std::abs(dtheta_1));

    // Central number density
    const real n_center = rho_center * density_to_n;

    std::cout << "\nPhysical scaling:" << std::endl;
    std::cout << "  α = " << alpha << " pc" << std::endl;
    std::cout << "  ρ_center = " << rho_center << " M_☉/pc³" << std::endl;
    std::cout << "  n_center = " << n_center << " cm⁻³" << std::endl;

    // ========================================================================
    // INITIALIZE K&I 2000 COOLING FOR T_eq LOOKUP
    // ========================================================================

    thermal::KoyamaInutsukaCooling ki_cooling;

    real T_center = ki_cooling.equilibrium_temperature(n_center, N_H_cm2);
    real T_edge = ki_cooling.equilibrium_temperature(0.1, N_H_cm2);  // Low density limit

    std::cout << "\nK&I 2000 thermal equilibrium:" << std::endl;
    std::cout << "  T_eq(center, n=" << n_center << ") = " << T_center << " K" << std::endl;
    std::cout << "  T_eq(edge, n~0.1) = " << T_edge << " K" << std::endl;

    // ========================================================================
    // WARNING: NOT IN HYDROSTATIC EQUILIBRIUM
    // ========================================================================

    std::cout << "\n*** WARNING ***" << std::endl;
    std::cout << "This IC is NOT in hydrostatic equilibrium!" << std::endl;
    std::cout << "Lane-Emden assumes T ∝ ρ^(2/3), but K&I gives T_eq(n)." << std::endl;
    std::cout << "The cloud will start to adjust when simulation runs." << std::endl;
    std::cout << "For tidal disruption: t_tidal << t_collapse, so this is OK." << std::endl;
    std::cout << "***************" << std::endl;

    // ========================================================================
    // INTERPOLATION HELPERS
    // ========================================================================

    auto get_theta = [&](real xi) -> real {
        if (xi >= xi_1) return 0.0;
        if (xi <= 0.0) return 1.0;

        for (size_t i = 1; i < xi_array.size(); ++i) {
            if (xi_array[i] >= xi) {
                real w = (xi - xi_array[i-1]) / (xi_array[i] - xi_array[i-1]);
                return theta_array[i-1] * (1.0 - w) + theta_array[i] * w;
            }
        }
        return 0.0;
    };

    // ========================================================================
    // BUILD CUMULATIVE MASS PROFILE
    // ========================================================================

    std::vector<real> M_cumulative(xi_array.size());
    M_cumulative[0] = 0.0;

    for (size_t i = 1; i < xi_array.size(); ++i) {
        real dxi = xi_array[i] - xi_array[i-1];
        real xi_mid = 0.5 * (xi_array[i] + xi_array[i-1]);
        real theta_mid = 0.5 * (theta_array[i] + theta_array[i-1]);
        real rho_mid = std::pow(theta_mid, 1.5);  // ρ/ρ_c = θ^1.5

        // dM = 4π α³ ρ_c × ξ² θ^1.5 dξ
        real dM = 4.0 * M_PI * alpha * alpha * alpha * rho_center *
                  xi_mid * xi_mid * rho_mid * dxi;
        M_cumulative[i] = M_cumulative[i-1] + dM;
    }

    // Normalize
    real M_total = M_cumulative.back();
    for (auto& m : M_cumulative) {
        m /= M_total;
    }

    // ========================================================================
    // CREATE PARTICLES
    // ========================================================================

    const int N_particles = use_glass ? (N * N * N) : N;
    const real particle_mass = M_cloud / N_particles;

    std::cout << "\nCreating particles:" << std::endl;
    std::cout << "  N_particles = " << N_particles << std::endl;
    std::cout << "  particle_mass = " << particle_mass << " M_☉" << std::endl;

    // Energy conversion: T[K] → u[code]
    // u = k_B T / ((γ-1) m_n) [CGS], then convert to code units
    const real energy_factor = k_B_cgs / ((gamma - 1.0) * m_n_cgs * kms_cgs * kms_cgs);

    auto& particles = m_sim->get_particles();
    particles.clear();
    particles.reserve(N_particles);

    if (use_glass) {
        // Glass-making method: random positions → radius mapping
        std::mt19937 gen(42);
        std::uniform_real_distribution<real> dis(-1.0, 1.0);

        std::vector<vec_t> positions;
        positions.reserve(N_particles);

        int count = 0;
        while (count < N_particles) {
#if DIM == 1
            vec_t pos(dis(gen));
#elif DIM == 2
            vec_t pos(dis(gen), dis(gen));
#else
            vec_t pos(dis(gen), dis(gen), dis(gen));
#endif
            real r_norm = std::abs(pos);
            if (r_norm <= 1.0 && r_norm > 0.0) {
                positions.push_back(pos);
                count++;
            }
        }

        // Map radii using cumulative mass
        auto interpolate_radius = [&](real mass_frac) -> real {
            for (size_t i = 1; i < M_cumulative.size(); ++i) {
                if (mass_frac >= M_cumulative[i-1] && mass_frac <= M_cumulative[i]) {
                    real frac = (mass_frac - M_cumulative[i-1]) /
                                (M_cumulative[i] - M_cumulative[i-1]);
                    real xi_interp = xi_array[i-1] + frac * (xi_array[i] - xi_array[i-1]);
                    return xi_interp * alpha;
                }
            }
            return R_cloud;
        };

        int particle_id = 0;
        for (const auto& pos : positions) {
            real r_uniform = std::abs(pos);
            real mass_frac = r_uniform * r_uniform * r_uniform;
            real r_le = interpolate_radius(mass_frac);

            vec_t mapped_pos = pos * (r_le / r_uniform);
            real r = std::abs(mapped_pos);

            if (r > R_cloud) continue;

            // Lane-Emden density
            real xi = r / alpha;
            real theta = get_theta(xi);
            if (theta <= 0.0) continue;

            real rho_local = rho_center * std::pow(theta, 1.5);
            if (rho_local <= 1e-20) continue;

            // K&I 2000 equilibrium temperature
            real n_local = rho_local * density_to_n;
            n_local = std::max(n_local, 0.1);  // Clamp to table range
            n_local = std::min(n_local, 1.0e4);

            real T_local = ki_cooling.equilibrium_temperature(n_local, N_H_cm2);

            // Internal energy from K&I temperature (NOT polytropic!)
            real u_local = T_local * energy_factor;

            // Pressure from ideal gas: P = (γ-1) ρ u
            real pres_local = (gamma - 1.0) * rho_local * u_local;

            // Smoothing length estimate
            constexpr real A = 4.0 * M_PI / 3.0;
            const int N_neighbor = m_param->physics.neighbor_number;
            real h_est = std::pow(N_neighbor * particle_mass / (rho_local * A), 1.0 / 3.0);

            SPHParticle p;
            p.pos = mapped_pos;
            p.vel = 0.0;
            p.mass = particle_mass;
            p.dens = rho_local;
            p.ene = u_local;
            p.pres = pres_local;
            p.sml = h_est;
            p.id = particle_id++;

            particles.push_back(p);
        }
    } else {
        // Shell distribution method
        const int N_shells = static_cast<int>(std::cbrt(N_particles));
        int particle_id = 0;

        for (int shell = 0; shell < N_shells; ++shell) {
            real r_shell = R_cloud * (shell + 0.5) / N_shells;
            int particles_this_shell = N_particles / N_shells;

            for (int i = 0; i < particles_this_shell; ++i) {
                // Fibonacci sphere
                const real golden = (1.0 + std::sqrt(5.0)) / 2.0;
                real phi = 2.0 * M_PI * i / golden;
                real z = 1.0 - 2.0 * (i + 0.5) / particles_this_shell;
                real r_xy = std::sqrt(1.0 - z * z);

#if DIM == 1
                vec_t pos(r_shell * (i % 2 == 0 ? 1.0 : -1.0));
#elif DIM == 2
                vec_t pos(r_shell * r_xy * std::cos(phi),
                          r_shell * r_xy * std::sin(phi));
#else
                vec_t pos(r_shell * r_xy * std::cos(phi),
                          r_shell * r_xy * std::sin(phi),
                          r_shell * z);
#endif

                real r = std::abs(pos);
                real xi = r / alpha;
                real theta = get_theta(xi);

                if (theta <= 0.0) continue;

                real rho_local = rho_center * std::pow(theta, 1.5);
                if (rho_local <= 1e-20) continue;

                real n_local = std::max(0.1, std::min(rho_local * density_to_n, 1.0e4));
                real T_local = ki_cooling.equilibrium_temperature(n_local, N_H_cm2);
                real u_local = T_local * energy_factor;
                real pres_local = (gamma - 1.0) * rho_local * u_local;

                constexpr real A = 4.0 * M_PI / 3.0;
                real h_est = std::pow(m_param->physics.neighbor_number * particle_mass /
                                      (rho_local * A), 1.0 / 3.0);

                SPHParticle p;
                p.pos = pos;
                p.vel = 0.0;
                p.mass = particle_mass;
                p.dens = rho_local;
                p.ene = u_local;
                p.pres = pres_local;
                p.sml = h_est;
                p.id = particle_id++;

                particles.push_back(p);
            }
        }
    }

    std::cout << "  Created " << particles.size() << " particles" << std::endl;

    // ========================================================================
    // STORE PARAMETERS FOR LATER USE
    // ========================================================================

    m_sample_parameters["R_cloud"] = R_cloud;
    m_sample_parameters["M_cloud"] = M_cloud;
    m_sample_parameters["rho_center"] = rho_center;
    m_sample_parameters["n_center"] = n_center;
    m_sample_parameters["alpha"] = alpha;
    m_sample_parameters["density_to_n"] = density_to_n;
    m_sample_parameters["N_H_cm2"] = N_H_cm2;

    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());

    std::cout << "\n=== Lane-Emden + K&I 2000 IC Complete ===" << std::endl;
    std::cout << "Ready for simulation with gravity + K&I cooling." << std::endl;
}

} // namespace sph
