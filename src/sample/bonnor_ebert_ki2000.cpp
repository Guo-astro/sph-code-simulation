/**
 * @file bonnor_ebert_ki2000.cpp
 * @brief Bonnor-Ebert sphere IC using Koyama & Inutsuka (2000) barotropic EOS
 *
 * Creates a pressure-truncated isothermal sphere with T_eq(n) from K&I 2000 tables.
 * Designed for analytical relaxation using KoyamaInutsukaRelaxation module.
 *
 * Physical model:
 *   - P(ρ) = n k_B T_eq(n)  [barotropic EOS from K&I 2000]
 *   - c_eff² = (k_B T / μ m_H) × (1 + d ln T / d ln n)
 *   - Hydrostatic: dP/dr = -ρ G M(r) / r²
 *   - Truncation: P(R_cloud) = P_ext
 *
 * Unit system:
 *   - Mass: M_☉ (solar mass)
 *   - Length: pc (parsec)
 *   - Velocity: km/s
 *   - G = 1.0 in these units (approximately)
 */

#include "solver.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "thermal/koyama_inutsuka_cooling.hpp"
#include "relaxation/koyama_inutsuka_relaxation.hpp"
#include "sample/ghost_envelope.hpp"
#include "exception.hpp"
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <algorithm>
#include <fstream>

namespace sph
{

// Physical constants
static constexpr real k_B_cgs = 1.380649e-16;      // erg/K
static constexpr real m_proton_cgs = 1.6726219e-24; // g
static constexpr real mu = 1.27;                    // Mean molecular weight for ISM
static constexpr real m_n_cgs = mu * m_proton_cgs;  // Mean particle mass [g]

// Unit conversions (code → CGS)
// Code units: M[M_☉], L[pc], V[km/s], t[pc/(km/s)]
static constexpr real Msun_cgs = 1.989e33;          // g
static constexpr real pc_cgs = 3.086e18;            // cm
static constexpr real kms_cgs = 1.0e5;              // cm/s

/**
 * @brief Create Bonnor-Ebert sphere IC with K&I 2000 thermal equilibrium
 *
 * Parameters from m_sample_parameters:
 *   - N: Grid resolution (N³ particles in initial cube)
 *   - R_cloud_pc: Cloud radius [pc]
 *   - M_cloud_Msun: Cloud mass [M_☉]
 *   - rho_center_nH: Central density [cm^-3] (number density)
 *   - N_H_cm2: Column density [cm^-2] for K&I tables (1e19 or 1e20)
 *   - P_ext_K_cm3: External pressure P/k_B [K cm^-3]
 */
void Solver::make_bonnor_ebert_ki2000()
{
#if DIM != 3
    THROW_ERROR("Bonnor-Ebert K&I 2000 test requires DIM == 3");
#endif

    // ========================================================================
    // READ PARAMETERS
    // ========================================================================

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);

    // Cloud physical parameters (with defaults for typical ISM cloud)
    const real R_cloud_pc = boost::any_cast<real>(m_sample_parameters.count("R_cloud_pc") ?
                                m_sample_parameters["R_cloud_pc"] : boost::any(1.0));  // 1 pc
    const real M_cloud_Msun = boost::any_cast<real>(m_sample_parameters.count("M_cloud_Msun") ?
                                m_sample_parameters["M_cloud_Msun"] : boost::any(1000.0));  // 1000 M_☉
    const real rho_center_nH = boost::any_cast<real>(m_sample_parameters.count("rho_center_nH") ?
                                m_sample_parameters["rho_center_nH"] : boost::any(1000.0));  // 1000 cm^-3
    const real N_H_cm2 = boost::any_cast<real>(m_sample_parameters.count("N_H_cm2") ?
                                m_sample_parameters["N_H_cm2"] : boost::any(1.0e19));  // 10^19 cm^-2
    const real P_ext_K_cm3 = boost::any_cast<real>(m_sample_parameters.count("P_ext_K_cm3") ?
                                m_sample_parameters["P_ext_K_cm3"] : boost::any(3000.0));  // 3000 K cm^-3

    // Ghost envelope parameters (for pressure confinement during relaxation)
    const bool use_envelope = m_sample_parameters.count("useEnvelope") ?
                              boost::any_cast<bool>(m_sample_parameters["useEnvelope"]) : false;
    const int envelope_layers = m_sample_parameters.count("envelopeLayers") ?
                                boost::any_cast<int>(m_sample_parameters["envelopeLayers"]) : 4;

    const real gamma = m_param->physics.gamma;
    const real G_code = m_param->gravity.constant;

    std::cout << "==========================================================" << std::endl;
    std::cout << "  Bonnor-Ebert K&I 2000 Initial Conditions" << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << "Cloud parameters:" << std::endl;
    std::cout << "  R_cloud = " << R_cloud_pc << " pc" << std::endl;
    std::cout << "  M_cloud = " << M_cloud_Msun << " M_☉" << std::endl;
    std::cout << "  ρ_center = " << rho_center_nH << " cm^-3 (n_H)" << std::endl;
    std::cout << "  N_H = " << N_H_cm2 << " cm^-2" << std::endl;
    std::cout << "  P_ext/k_B = " << P_ext_K_cm3 << " K cm^-3" << std::endl;
    std::cout << "  γ = " << gamma << std::endl;
    std::cout << "  G (code) = " << G_code << std::endl;
    std::cout << std::endl;

    // ========================================================================
    // UNIT CONVERSIONS
    // ========================================================================

    // Code units: L[pc], M[M_☉], V[km/s]
    // Density: code → CGS
    // ρ[g/cm³] = M[M_☉] / L[pc]³ × (Msun/pc³)
    const real density_code_to_cgs = Msun_cgs / pow3(pc_cgs);  // g/cm³ per code

    // n_H[cm^-3] = ρ[g/cm³] / m_n[g]
    // density_to_n = (density_code_to_cgs / m_n_cgs)
    const real density_to_n = density_code_to_cgs / m_n_cgs;  // code density → n_H [cm^-3]

    std::cout << "Unit conversions:" << std::endl;
    std::cout << "  density_code → g/cm³: " << density_code_to_cgs << std::endl;
    std::cout << "  density_code → n_H[cm^-3]: " << density_to_n << std::endl;

    // Central density in code units
    const real rho_center_code = rho_center_nH / density_to_n;
    std::cout << "  ρ_center (code) = " << rho_center_code << std::endl;

    // ========================================================================
    // INITIALIZE K&I 2000 COOLING FOR T_eq(n) LOOKUP
    // ========================================================================

    thermal::KoyamaInutsukaCooling ki_cooling;

    // Get equilibrium temperature at center
    real T_center = ki_cooling.equilibrium_temperature(rho_center_nH, N_H_cm2);
    std::cout << "  T_eq(center) = " << T_center << " K" << std::endl;

    // Central pressure
    real P_center_cgs = rho_center_nH * k_B_cgs * T_center;  // dyn/cm²
    real P_center_K_cm3 = rho_center_nH * T_center;  // K cm^-3
    std::cout << "  P(center)/k_B = " << P_center_K_cm3 << " K cm^-3" << std::endl;

    // ========================================================================
    // BUILD EQUILIBRIUM PROFILE USING RELAXATION MODULE (SSOT)
    // ========================================================================
    // Use KoyamaInutsukaRelaxation to compute the profile - this ensures
    // IC and relaxation use EXACTLY the same profile (Single Source of Truth)

    std::cout << "\nComputing Bonnor-Ebert profile using SSOT relaxation module..." << std::endl;

    // Initialize relaxation module to compute profile
    KIRelaxationParams ki_params;
    ki_params.R_cloud = R_cloud_pc;
    ki_params.M_cloud = M_cloud_Msun;
    ki_params.P_ext = P_ext_K_cm3;  // K cm^-3
    ki_params.rho_center = rho_center_code;
    ki_params.N_H = N_H_cm2;
    ki_params.G = G_code;
    ki_params.gamma = gamma;
    ki_params.density_to_n = density_to_n;

    KoyamaInutsukaRelaxation ki_profile;
    ki_profile.initialize(ki_params);

    // Get actual truncation radius and mass from profile
    const real R_actual = ki_profile.get_R_cloud();
    const real M_actual = ki_profile.get_M_cloud();

    // Get radius table for building cumulative mass profile
    const auto& r_table = ki_profile.get_r_table();
    const int n_profile = r_table.size();

    std::cout << "Profile from SSOT relaxation module:" << std::endl;
    std::cout << "  R_actual = " << R_actual << " [code]" << std::endl;
    std::cout << "  M_actual = " << M_actual << " [code]" << std::endl;
    std::cout << "  Profile points: " << n_profile << std::endl;

    // ========================================================================
    // BUILD CUMULATIVE MASS PROFILE FOR PARTICLE PLACEMENT
    // ========================================================================

    std::vector<real> M_cumulative(n_profile);
    M_cumulative[0] = 0.0;
    for (int i = 1; i < n_profile; ++i) {
        M_cumulative[i] = ki_profile.get_M_enclosed(r_table[i]);
    }

    // Normalize to actual mass
    real M_total = M_cumulative.back();
    for (auto& m : M_cumulative) {
        m /= M_total;
    }

    // ========================================================================
    // CREATE PARTICLES USING GLASS-MAKING METHOD
    // ========================================================================

    const int N_particles = N * N * N;
    const real particle_mass = M_actual / N_particles;

    std::cout << "\nCreating particles:" << std::endl;
    std::cout << "  Target N = " << N_particles << std::endl;
    std::cout << "  Particle mass = " << particle_mass << " M_☉" << std::endl;

    // Step 1: Generate random points in unit sphere
    std::vector<vec_t> positions;
    positions.reserve(N_particles);

    std::mt19937 gen(42);  // Fixed seed for reproducibility
    std::uniform_real_distribution<real> dis(-1.0, 1.0);

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

    // Step 2: Map radii to Bonnor-Ebert profile using cumulative mass
    std::vector<vec_t> mapped_positions;
    mapped_positions.reserve(N_particles);

    auto interpolate_radius = [&](real mass_frac) -> real {
        // Find r where M(<r)/M_total = mass_frac
        for (int i = 1; i < n_profile; ++i) {
            if (mass_frac >= M_cumulative[i-1] && mass_frac <= M_cumulative[i]) {
                real frac = (mass_frac - M_cumulative[i-1]) / (M_cumulative[i] - M_cumulative[i-1]);
                return r_table[i-1] + frac * (r_table[i] - r_table[i-1]);
            }
        }
        return r_table.back();
    };

    // Use ki_profile interpolation methods for density and temperature
    auto interpolate_density = [&](real r) -> real {
        return ki_profile.get_rho_eq(r);
    };

    auto interpolate_temperature = [&](real r) -> real {
        return ki_profile.get_T_eq(r);
    };

    for (const auto& pos : positions) {
        real r_uniform = std::abs(pos);
        real mass_frac = r_uniform * r_uniform * r_uniform;  // M ∝ r³ for uniform sphere

        real r_be = interpolate_radius(mass_frac);
        real scale = r_be / r_uniform;

        vec_t mapped_pos = pos * scale;
        mapped_positions.push_back(mapped_pos);
    }

    // ========================================================================
    // CREATE SPH PARTICLES
    // ========================================================================

    auto& particles = m_sim->get_particles();
    particles.clear();
    particles.reserve(N_particles);

    // Energy conversion: T[K] → u[code]
    // u = (k_B T) / ((γ-1) m_n) [CGS: erg/g]
    // u[code] = u[CGS] / (km/s)² = u[CGS] / 1e10
    const real energy_factor = k_B_cgs / ((gamma - 1.0) * m_n_cgs * kms_cgs * kms_cgs);

    std::cout << "  Energy factor (T → u): " << energy_factor << std::endl;

    int particle_id = 0;
    for (const auto& pos : mapped_positions) {
        real r = std::abs(pos);
        if (r > R_actual) continue;

        // Get equilibrium properties
        real rho_local = interpolate_density(r);
        real T_local = interpolate_temperature(r);

        // Internal energy from equilibrium temperature
        real u_local = T_local * energy_factor;

        // Smoothing length estimate
        constexpr real A = 4.0 * M_PI / 3.0;  // 3D sphere
        const int N_neighbor = m_param->physics.neighbor_number;
        real h_est = std::pow(N_neighbor * particle_mass / (rho_local * A), 1.0 / 3.0);

        SPHParticle p;
        p.pos = pos;
        p.vel = 0.0;  // Initially at rest
        p.mass = particle_mass;
        p.dens = rho_local;
        p.ene = u_local;
        p.pres = (gamma - 1.0) * rho_local * u_local;
        p.sml = h_est;
        p.id = particle_id++;

        particles.push_back(p);
    }

    std::cout << "  Created " << particles.size() << " cloud (GAS) particles" << std::endl;

    // ========================================================================
    // CREATE GHOST ENVELOPE (optional, for pressure confinement)
    // ========================================================================

    if (use_envelope) {
        // Get edge density from profile at R_actual
        real rho_edge_code = ki_profile.get_rho_eq(R_actual);

        // Get edge temperature and compute internal energy
        real T_edge = ki_profile.get_T_eq(R_actual);
        real u_edge = T_edge * energy_factor;

        // Configure envelope using SSOT module
        GhostEnvelopeConfig env_config;
        env_config.R_cloud = R_actual;
        env_config.rho_edge = rho_edge_code;
        env_config.u_envelope = u_edge;
        env_config.particle_mass = particle_mass;
        env_config.N_neighbor = m_param->physics.neighbor_number;
        env_config.num_layers = envelope_layers;

        // Generate ghost envelope particles
        auto envelope_particles = GhostEnvelopeGenerator::generate(env_config);

        // Assign IDs continuing from cloud particles
        for (auto& p : envelope_particles) {
            p.id = particle_id++;
        }

        // Append envelope to particles
        particles.insert(particles.end(), envelope_particles.begin(), envelope_particles.end());

        // Print summary
        GhostEnvelopeGenerator::print_summary(env_config, envelope_particles.size());
    }

    // ========================================================================
    // STORE PARAMETERS FOR RELAXATION MODULE
    // ========================================================================

    m_sample_parameters["R_cloud"] = R_actual;
    m_sample_parameters["M_cloud"] = M_actual;
    m_sample_parameters["P_ext_K_cm3"] = P_ext_K_cm3;
    m_sample_parameters["rho_center_code"] = rho_center_code;
    m_sample_parameters["N_H_cm2"] = N_H_cm2;
    m_sample_parameters["density_to_n"] = density_to_n;

    // ========================================================================
    // SAVE PROFILE FOR RELAXATION MODULE (SSOT - Single Source Of Truth)
    // ========================================================================
    // Save the SAME profile that was used for IC particle placement
    // This ensures IC and relaxation use EXACTLY the same profile
    {
        // Save profile to output directory for relaxation to load
        std::string profile_file = m_output_dir + "/ki2000_profile.dat";
        ki_profile.save_profile_to_file(profile_file);

        // Store profile path for relaxation module
        m_sample_parameters["ki2000_profile_file"] = profile_file;

        std::cout << "\n=== SSOT Profile Saved ===" << std::endl;
        std::cout << "Profile file: " << profile_file << std::endl;
        std::cout << "IC and Relaxation use same profile with R_cloud = " << R_actual << " [code]" << std::endl;
    }

    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());

    std::cout << "\n=== Bonnor-Ebert K&I 2000 IC Complete ===" << std::endl;
    std::cout << "Total particles: " << particles.size() << std::endl;
    std::cout << "Ready for analytical relaxation phase." << std::endl;
}

} // namespace sph
