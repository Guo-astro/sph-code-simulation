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
#include "exception.hpp"
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <algorithm>

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
    // BUILD EQUILIBRIUM PROFILE BY INTEGRATING HYDROSTATIC ODE
    // ========================================================================
    //
    // Hydrostatic equilibrium: dP/dr = -ρ G M(r) / r²
    // With barotropic EOS: P(ρ) = n k_B T_eq(n)
    //
    // Rewrite as: dρ/dr = -(dP/dρ)⁻¹ × ρ G M(r) / r²
    //           = -ρ G M(r) / (r² c_eff²)
    // where c_eff² = dP/dρ = (k_B T / m_n) × (1 + d ln T / d ln n)

    std::cout << "\nIntegrating Bonnor-Ebert profile..." << std::endl;

    // Profile arrays
    std::vector<real> r_profile, rho_profile, M_profile, T_profile;
    r_profile.reserve(2000);
    rho_profile.reserve(2000);
    M_profile.reserve(2000);
    T_profile.reserve(2000);

    // Helper function: get effective sound speed squared
    auto get_c_eff_squared = [&](real rho_code) -> real {
        real n_H = rho_code * density_to_n;
        n_H = std::max(n_H, 0.1);
        n_H = std::min(n_H, 1.0e4);

        real T_eq = ki_cooling.equilibrium_temperature(n_H, N_H_cm2);

        // Numerical derivative d ln T / d ln n
        const real eps = 0.02;
        real n_lo = n_H * (1.0 - eps);
        real n_hi = n_H * (1.0 + eps);
        n_lo = std::max(n_lo, 0.1);
        n_hi = std::min(n_hi, 1.0e4);

        real T_lo = ki_cooling.equilibrium_temperature(n_lo, N_H_cm2);
        real T_hi = ki_cooling.equilibrium_temperature(n_hi, N_H_cm2);

        real dlnT_dlnn = 0.0;
        if (T_lo > 0.0 && T_hi > 0.0) {
            dlnT_dlnn = (std::log(T_hi) - std::log(T_lo)) / (std::log(n_hi) - std::log(n_lo));
        }

        // c_eff² in CGS [cm²/s²]
        real c_eff_sq_cgs = (k_B_cgs * T_eq / m_n_cgs) * (1.0 + dlnT_dlnn);

        // Convert to code units: [cm/s]² → [km/s]²
        real c_eff_sq_code = c_eff_sq_cgs / (kms_cgs * kms_cgs);

        return std::max(c_eff_sq_code, 1.0e-6);
    };

    // Initial conditions at center
    const real dr = 0.002 * R_cloud_pc;
    real r = 1.0e-6 * R_cloud_pc;  // Start at small r
    real rho = rho_center_code;
    real M_enc = (4.0 / 3.0) * M_PI * rho * r * r * r;

    // External pressure for truncation
    real P_ext_cgs = P_ext_K_cm3 * k_B_cgs;  // dyn/cm²

    real R_actual = R_cloud_pc;
    real M_actual = M_cloud_Msun;
    bool truncated = false;

    // RK4 integration
    const int max_steps = 5000;
    const real r_max = 5.0 * R_cloud_pc;

    for (int step = 0; step < max_steps && r < r_max; ++step) {
        // Store current state
        real n_H = rho * density_to_n;
        real T_eq = ki_cooling.equilibrium_temperature(std::max(n_H, 0.1), N_H_cm2);
        real P_cgs = n_H * k_B_cgs * T_eq;

        r_profile.push_back(r);
        rho_profile.push_back(rho);
        M_profile.push_back(M_enc);
        T_profile.push_back(T_eq);

        // Check truncation
        if (P_cgs <= P_ext_cgs && step > 10) {
            truncated = true;
            R_actual = r;
            M_actual = M_enc;
            break;
        }

        // RK4 step for (rho, M)
        real c_eff_sq = get_c_eff_squared(rho);

        // k1
        real drho_dr = -rho * G_code * M_enc / (r * r * c_eff_sq);
        real dM_dr = 4.0 * M_PI * r * r * rho;
        real k1_rho = dr * drho_dr;
        real k1_M = dr * dM_dr;

        // k2
        real r2 = r + 0.5 * dr;
        real rho2 = std::max(rho + 0.5 * k1_rho, 1.0e-30);
        real M2 = M_enc + 0.5 * k1_M;
        real c_eff_sq_2 = get_c_eff_squared(rho2);
        real k2_rho = dr * (-rho2 * G_code * M2 / (r2 * r2 * c_eff_sq_2));
        real k2_M = dr * (4.0 * M_PI * r2 * r2 * rho2);

        // k3
        real rho3 = std::max(rho + 0.5 * k2_rho, 1.0e-30);
        real M3 = M_enc + 0.5 * k2_M;
        real c_eff_sq_3 = get_c_eff_squared(rho3);
        real k3_rho = dr * (-rho3 * G_code * M3 / (r2 * r2 * c_eff_sq_3));
        real k3_M = dr * (4.0 * M_PI * r2 * r2 * rho3);

        // k4
        real r4 = r + dr;
        real rho4 = std::max(rho + k3_rho, 1.0e-30);
        real M4 = M_enc + k3_M;
        real c_eff_sq_4 = get_c_eff_squared(rho4);
        real k4_rho = dr * (-rho4 * G_code * M4 / (r4 * r4 * c_eff_sq_4));
        real k4_M = dr * (4.0 * M_PI * r4 * r4 * rho4);

        // Update
        rho += (k1_rho + 2.0 * k2_rho + 2.0 * k3_rho + k4_rho) / 6.0;
        M_enc += (k1_M + 2.0 * k2_M + 2.0 * k3_M + k4_M) / 6.0;
        r += dr;

        rho = std::max(rho, 1.0e-30);

        // Stop if density becomes too low
        if (rho * density_to_n < 0.01) {
            R_actual = r;
            M_actual = M_enc;
            break;
        }
    }

    if (!truncated) {
        std::cout << "Warning: Profile not truncated at P_ext, using final values" << std::endl;
    }

    std::cout << "Profile integration complete:" << std::endl;
    std::cout << "  R_actual = " << R_actual << " pc" << std::endl;
    std::cout << "  M_actual = " << M_actual << " M_☉" << std::endl;
    std::cout << "  Profile points: " << r_profile.size() << std::endl;

    // ========================================================================
    // BUILD CUMULATIVE MASS PROFILE FOR PARTICLE PLACEMENT
    // ========================================================================

    std::vector<real> M_cumulative(r_profile.size());
    M_cumulative[0] = 0.0;
    for (size_t i = 1; i < r_profile.size(); ++i) {
        M_cumulative[i] = M_profile[i];
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
        for (size_t i = 1; i < M_cumulative.size(); ++i) {
            if (mass_frac >= M_cumulative[i-1] && mass_frac <= M_cumulative[i]) {
                real frac = (mass_frac - M_cumulative[i-1]) / (M_cumulative[i] - M_cumulative[i-1]);
                return r_profile[i-1] + frac * (r_profile[i] - r_profile[i-1]);
            }
        }
        return r_profile.back();
    };

    auto interpolate_density = [&](real r) -> real {
        for (size_t i = 1; i < r_profile.size(); ++i) {
            if (r >= r_profile[i-1] && r <= r_profile[i]) {
                real frac = (r - r_profile[i-1]) / (r_profile[i] - r_profile[i-1]);
                return rho_profile[i-1] + frac * (rho_profile[i] - rho_profile[i-1]);
            }
        }
        return rho_profile.back();
    };

    auto interpolate_temperature = [&](real r) -> real {
        for (size_t i = 1; i < r_profile.size(); ++i) {
            if (r >= r_profile[i-1] && r <= r_profile[i]) {
                real frac = (r - r_profile[i-1]) / (r_profile[i] - r_profile[i-1]);
                return T_profile[i-1] + frac * (T_profile[i] - T_profile[i-1]);
            }
        }
        return T_profile.back();
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

    std::cout << "  Created " << particles.size() << " particles" << std::endl;

    // ========================================================================
    // STORE PARAMETERS FOR RELAXATION MODULE
    // ========================================================================

    m_sample_parameters["R_cloud"] = R_actual;
    m_sample_parameters["M_cloud"] = M_actual;
    m_sample_parameters["P_ext_K_cm3"] = P_ext_K_cm3;
    m_sample_parameters["rho_center_code"] = rho_center_code;
    m_sample_parameters["N_H_cm2"] = N_H_cm2;
    m_sample_parameters["density_to_n"] = density_to_n;

    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());

    std::cout << "\n=== Bonnor-Ebert K&I 2000 IC Complete ===" << std::endl;
    std::cout << "Ready for analytical relaxation phase." << std::endl;
}

} // namespace sph
