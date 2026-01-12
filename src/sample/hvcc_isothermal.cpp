/**
 * @file hvcc_isothermal.cpp
 * @brief HVCC initial conditions with pure isothermal Bonnor-Ebert sphere
 *
 * Implements isothermal cloud IC at constant T = 10K for IMBH-cloud interaction.
 * Uses proper Bonnor-Ebert density profile from hydrostatic equilibrium:
 *   dP/dr = -ρ G M(r)/r²  with  P = ρ c_s²  (constant c_s)
 *
 * This differs from K&I barotropic EOS which has T = T_eq(n).
 * Here T is CONSTANT, giving the classical isothermal BE sphere.
 *
 * Reference:
 *   - Bonnor (1956), MNRAS 116, 351
 *   - Ebert (1955), ZA 37, 217
 *   - Oka et al. (2017) - HVCC CO-0.40-0.22 candidate IMBH
 */

#include "solver.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "thermal/koyama_inutsuka_data.hpp"
#include "sample/ghost_envelope.hpp"
#include "sample/hvcc_isothermal.hpp"
#include "exception.hpp"
#include <vector>
#include <cmath>
#include <random>
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
}

// ============================================================================
// Helper function implementations
// ============================================================================
namespace sample {

real find_equilibrium_density(real T_target, real N_H)
{
    // For pure isothermal, this is not used - included for interface compatibility
    return 100.0;  // Default value
}

real isothermal_sound_speed(real T, real mu)
{
    // c_s = sqrt(k_B T / (mu m_H)) in CGS, then convert to km/s
    real c_s_cgs = std::sqrt(k_B_cgs * T / (mu * m_proton_cgs));
    return c_s_cgs / kms_cgs;  // km/s
}

real n_to_code_density(real n_H)
{
    return n_H / density_to_n_factor;  // M_sun/pc^3
}

real code_density_to_n(real rho)
{
    return rho * density_to_n_factor;  // cm^-3
}

} // namespace sample

// ============================================================================
// Isothermal Bonnor-Ebert sphere profile integration
// ============================================================================
namespace {

struct BEProfile {
    std::vector<real> r;      // Radius [code]
    std::vector<real> rho;    // Density [code]
    std::vector<real> M_enc;  // Enclosed mass [code]
    real R_cloud;             // Cloud radius where P = P_ext
    real M_cloud;             // Total cloud mass
};

/**
 * @brief Integrate isothermal Bonnor-Ebert hydrostatic equilibrium
 *
 * Equations:
 *   dρ/dr = -ρ G M(r) / (r² c_s²)
 *   dM/dr = 4π r² ρ
 *
 * Boundary conditions:
 *   ρ(0) = ρ_center
 *   M(0) = 0
 *
 * Truncation:
 *   P(R) = ρ(R) c_s² = P_ext
 *   → ρ(R) = P_ext / c_s²
 *
 * @param rho_center Central density [code units]
 * @param c_s Isothermal sound speed [code units = km/s]
 * @param G Gravitational constant [code units]
 * @param P_ext_code External pressure [code units]
 * @param dr Integration step [code units]
 * @return BEProfile structure with density profile
 */
BEProfile integrate_isothermal_BE(real rho_center, real c_s, real G,
                                   real P_ext_code, real dr)
{
    BEProfile profile;
    profile.r.reserve(5000);
    profile.rho.reserve(5000);
    profile.M_enc.reserve(5000);

    const real c_s_sq = c_s * c_s;
    const real rho_edge = P_ext_code / c_s_sq;  // Density at truncation

    // Start from small radius (avoid r=0 singularity)
    real r = 1.0e-6;
    real rho = rho_center;
    real M = (4.0 / 3.0) * M_PI * rho * r * r * r;

    profile.r.push_back(r);
    profile.rho.push_back(rho);
    profile.M_enc.push_back(M);

    const int max_steps = 100000;
    const real r_max = 10.0;  // Maximum radius to integrate [pc]

    for (int step = 0; step < max_steps && r < r_max; ++step) {
        // Check truncation condition: ρ <= ρ_edge
        if (rho <= rho_edge && step > 10) {
            profile.R_cloud = r;
            profile.M_cloud = M;
            break;
        }

        // RK4 integration for (ρ, M)
        // dρ/dr = -ρ G M / (r² c_s²)
        // dM/dr = 4π r² ρ

        auto drho_dr = [&](real r_val, real rho_val, real M_val) -> real {
            if (r_val < 1.0e-10) return 0.0;
            return -rho_val * G * M_val / (r_val * r_val * c_s_sq);
        };

        auto dM_dr = [](real r_val, real rho_val) -> real {
            return 4.0 * M_PI * r_val * r_val * rho_val;
        };

        // k1
        real k1_rho = dr * drho_dr(r, rho, M);
        real k1_M = dr * dM_dr(r, rho);

        // k2
        real r2 = r + 0.5 * dr;
        real rho2 = std::max(rho + 0.5 * k1_rho, 1.0e-30);
        real M2 = M + 0.5 * k1_M;
        real k2_rho = dr * drho_dr(r2, rho2, M2);
        real k2_M = dr * dM_dr(r2, rho2);

        // k3
        real rho3 = std::max(rho + 0.5 * k2_rho, 1.0e-30);
        real M3 = M + 0.5 * k2_M;
        real k3_rho = dr * drho_dr(r2, rho3, M3);
        real k3_M = dr * dM_dr(r2, rho3);

        // k4
        real r4 = r + dr;
        real rho4 = std::max(rho + k3_rho, 1.0e-30);
        real M4 = M + k3_M;
        real k4_rho = dr * drho_dr(r4, rho4, M4);
        real k4_M = dr * dM_dr(r4, rho4);

        // Update
        rho += (k1_rho + 2.0 * k2_rho + 2.0 * k3_rho + k4_rho) / 6.0;
        M += (k1_M + 2.0 * k2_M + 2.0 * k3_M + k4_M) / 6.0;
        r += dr;

        rho = std::max(rho, 1.0e-30);

        profile.r.push_back(r);
        profile.rho.push_back(rho);
        profile.M_enc.push_back(M);

        // Stop if density becomes very low
        if (rho < 1.0e-10 * rho_center) {
            profile.R_cloud = r;
            profile.M_cloud = M;
            break;
        }
    }

    // If not truncated, use final values
    if (profile.R_cloud == 0.0) {
        profile.R_cloud = profile.r.back();
        profile.M_cloud = profile.M_enc.back();
    }

    return profile;
}

/**
 * @brief Interpolate density at given radius from BE profile
 */
real interpolate_density(const BEProfile& profile, real r)
{
    if (r <= profile.r.front()) return profile.rho.front();
    if (r >= profile.r.back()) return profile.rho.back();

    for (size_t i = 1; i < profile.r.size(); ++i) {
        if (r <= profile.r[i]) {
            real frac = (r - profile.r[i-1]) / (profile.r[i] - profile.r[i-1]);
            return profile.rho[i-1] + frac * (profile.rho[i] - profile.rho[i-1]);
        }
    }
    return profile.rho.back();
}

/**
 * @brief Interpolate radius for given mass fraction from BE profile
 */
real interpolate_radius_from_mass(const BEProfile& profile, real mass_frac)
{
    real M_total = profile.M_cloud;
    real M_target = mass_frac * M_total;

    for (size_t i = 1; i < profile.M_enc.size(); ++i) {
        if (M_target <= profile.M_enc[i]) {
            real frac = (M_target - profile.M_enc[i-1]) /
                        (profile.M_enc[i] - profile.M_enc[i-1]);
            return profile.r[i-1] + frac * (profile.r[i] - profile.r[i-1]);
        }
    }
    return profile.R_cloud;
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
    points.reserve(N_target * 2);  // Over-allocate

    // HCP lattice parameters
    // In HCP: a = spacing in x-y plane, c = 2*sqrt(2/3)*a (z spacing)
    // Estimate spacing from target particle count
    // Volume of unit sphere = 4π/3, each particle "owns" ~V/N volume
    // For HCP packing efficiency ~0.74, effective spacing: a³ ≈ 4π/(3*N*0.74)
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
        // Offset every other z layer
        real y_offset = (iz % 2 == 0) ? 0.0 : dy / 3.0;
        real x_offset = (iz % 2 == 0) ? 0.0 : a / 2.0;

        for (int iy = -ny; iy <= ny; ++iy) {
            real y = iy * dy + y_offset;
            // Offset every other y row
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

    // If we have too many points, randomly subsample
    if (points.size() > static_cast<size_t>(N_target * 1.1)) {
        std::mt19937 gen(42);
        std::shuffle(points.begin(), points.end(), gen);
        points.resize(N_target);
    }

    std::cout << "  HCP lattice: generated " << points.size() << " points (target " << N_target << ")" << std::endl;
    return points;
}

/**
 * @brief Map uniform lattice to density profile using cumulative particle count
 *
 * For each lattice point ranked by radius:
 *   mass_frac = particle_rank / N_total
 *   r_profile = radius in BE profile with same enclosed mass fraction
 *
 * This is more accurate than using r³ because it uses the actual
 * cumulative distribution of lattice particles, not the theoretical
 * uniform sphere assumption.
 */
void stretch_to_density_profile(std::vector<vec_t>& points, const BEProfile& profile)
{
    // Sort points by radius to get cumulative rank
    std::vector<std::pair<real, size_t>> radii_indices;
    radii_indices.reserve(points.size());

    for (size_t i = 0; i < points.size(); ++i) {
        real r = std::abs(points[i]);
        radii_indices.push_back({r, i});
    }

    std::sort(radii_indices.begin(), radii_indices.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    // Map each particle using its cumulative rank
    const real N_total = static_cast<real>(points.size());

    for (size_t rank = 0; rank < radii_indices.size(); ++rank) {
        size_t idx = radii_indices[rank].second;
        real r_unit = radii_indices[rank].first;

        if (r_unit < 1.0e-10) continue;

        // Mass fraction = cumulative particle count / total
        // This ensures each particle represents equal mass
        real mass_frac = (rank + 0.5) / N_total;  // +0.5 for bin center

        // Find radius in BE profile with same enclosed mass fraction
        real r_profile = interpolate_radius_from_mass(profile, mass_frac);

        // Scale position
        real scale = r_profile / r_unit;
        points[idx] *= scale;
    }
}
#endif // DIM == 3

} // anonymous namespace

// ============================================================================
// Main IC generator: make_hvcc_isothermal_10k()
// ============================================================================

void Solver::make_hvcc_isothermal_10k()
{
#if DIM != 3
    THROW_ERROR("HVCC isothermal 10K requires DIM == 3");
#else
    // ========================================================================
    // READ PARAMETERS
    // ========================================================================

    sample::HVCCIsothermalConfig config;

    // Override defaults from m_sample_parameters
    // N is grid resolution (same convention as bonnor_ebert_ki2000): N³ particles
    if (m_sample_parameters.count("N")) {
        int N = boost::any_cast<int>(m_sample_parameters["N"]);
        config.N_particles = N * N * N;
    }
    if (m_sample_parameters.count("M_cloud"))
        config.M_cloud = boost::any_cast<real>(m_sample_parameters["M_cloud"]);
    if (m_sample_parameters.count("R_cloud"))
        config.R_cloud = boost::any_cast<real>(m_sample_parameters["R_cloud"]);
    if (m_sample_parameters.count("T_cloud"))
        config.T_cloud = boost::any_cast<real>(m_sample_parameters["T_cloud"]);
    if (m_sample_parameters.count("n_center"))
        config.n_center = boost::any_cast<real>(m_sample_parameters["n_center"]);
    if (m_sample_parameters.count("N_H_cm2"))
        config.N_H = boost::any_cast<real>(m_sample_parameters["N_H_cm2"]);
    if (m_sample_parameters.count("M_BH"))
        config.M_BH = boost::any_cast<real>(m_sample_parameters["M_BH"]);
    if (m_sample_parameters.count("b_impact"))
        config.b_impact = boost::any_cast<real>(m_sample_parameters["b_impact"]);
    if (m_sample_parameters.count("v_rel"))
        config.v_rel = boost::any_cast<real>(m_sample_parameters["v_rel"]);
    if (m_sample_parameters.count("epsilon_BH"))
        config.epsilon_BH = boost::any_cast<real>(m_sample_parameters["epsilon_BH"]);
    if (m_sample_parameters.count("useEnvelope"))
        config.use_envelope = boost::any_cast<bool>(m_sample_parameters["useEnvelope"]);
    if (m_sample_parameters.count("envelopeLayers"))
        config.envelope_layers = boost::any_cast<int>(m_sample_parameters["envelopeLayers"]);
    if (m_sample_parameters.count("P_ext"))
        config.P_ext = boost::any_cast<real>(m_sample_parameters["P_ext"]);

    const real gamma = m_param->physics.gamma;
    const real G = m_param->gravity.constant;

    std::cout << "==========================================================" << std::endl;
    std::cout << "  Isothermal Bonnor-Ebert Sphere at T = " << config.T_cloud << " K" << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << std::endl;
    std::cout << "Physical model:" << std::endl;
    std::cout << "  - Pure isothermal EOS: P = rho * c_s^2 (constant T)" << std::endl;
    std::cout << "  - Hydrostatic equilibrium: dP/dr = -rho * G * M(r)/r^2" << std::endl;
    std::cout << "  - Truncation at P(R) = P_ext" << std::endl;
    std::cout << std::endl;

    // ========================================================================
    // COMPUTE SOUND SPEED AND EXTERNAL PRESSURE
    // ========================================================================

    // Isothermal sound speed at constant T
    const real c_s = sample::isothermal_sound_speed(config.T_cloud, mu_mol);
    const real c_s_sq = c_s * c_s;

    std::cout << "Isothermal parameters:" << std::endl;
    std::cout << "  T = " << config.T_cloud << " K (CONSTANT)" << std::endl;
    std::cout << "  mu = " << mu_mol << " (molecular H2)" << std::endl;
    std::cout << "  c_s = " << c_s << " km/s" << std::endl;
    std::cout << "  c_s^2 = " << c_s_sq << " (km/s)^2" << std::endl;
    std::cout << std::endl;

    // External pressure from K&I 2000 at T=10K
    // P_ext/k_B [K cm^-3] -> P_ext [code] = (P_ext/k_B) * k_B / (code pressure)
    // Code pressure: [M_sun/pc^3] * [km/s]^2 = [M_sun km^2 / pc^3 s^2]
    // P_ext_cgs = P_ext_K_cm3 * k_B_cgs  [dyn/cm^2 = erg/cm^3]
    // P_ext_code = P_ext_cgs / (code_density * code_velocity^2)
    //            = P_ext_cgs / (Msun/pc^3 * (km/s)^2)
    //            = P_ext_cgs / (Msun_cgs/pc_cgs^3 * kms_cgs^2)

    const real P_ext_cgs = config.P_ext * k_B_cgs;  // erg/cm^3
    const real pressure_code_to_cgs = density_code_to_cgs * kms_cgs * kms_cgs;
    const real P_ext_code = P_ext_cgs / pressure_code_to_cgs;

    std::cout << "External pressure (K&I 2000 at T=" << config.T_cloud << "K):" << std::endl;
    std::cout << "  P_ext/k_B = " << config.P_ext << " K cm^-3" << std::endl;
    std::cout << "  P_ext (code) = " << P_ext_code << std::endl;
    std::cout << std::endl;

    // ========================================================================
    // CREATE MODIFIED ISOTHERMAL SPHERE (BE-LIKE) PROFILE
    // ========================================================================
    //
    // Given n_center and M_cloud, we find R_cloud and r_core such that:
    //   ρ(r) = ρ_c / (1 + (r/r_c)²)
    //   M(<R) = M_cloud
    //
    // For a modified isothermal sphere:
    //   M(<R) = 4π ρ_c r_c³ [R/r_c - arctan(R/r_c)]
    // ========================================================================

    // Use n_center directly as specified by user
    real rho_center = sample::n_to_code_density(config.n_center);
    const real M_target = config.M_cloud;

    std::cout << "Target cloud parameters:" << std::endl;
    std::cout << "  M_cloud = " << M_target << " M_sun" << std::endl;
    std::cout << "  n_center = " << config.n_center << " cm^-3" << std::endl;
    std::cout << "  rho_center = " << rho_center << " M_sun/pc^3" << std::endl;
    std::cout << std::endl;

    // For modified isothermal sphere, choose r_c such that the edge density
    // matches K&I 2000 external pressure for perfect equilibrium:
    //   ρ(R) = ρ_c / (1 + (R/r_c)²) = ρ_ext
    //   => (R/r_c)² = ρ_c/ρ_ext - 1
    //   => R/r_c = sqrt(n_center/n_ext - 1)

    std::cout << "Creating modified isothermal sphere (BE-like) profile..." << std::endl;

    // K&I 2000 external density for equilibrium
    real n_ext = config.P_ext / config.T_cloud;  // cm^-3 (= 71 for P=710, T=10)
    real rho_ext = sample::n_to_code_density(n_ext);
    real density_contrast = config.n_center / n_ext;  // = 1000/71 ≈ 14.08

    std::cout << "  K&I 2000 equilibrium:" << std::endl;
    std::cout << "    n_ext = P_ext/T = " << n_ext << " cm^-3" << std::endl;
    std::cout << "    Density contrast n_c/n_ext = " << density_contrast << std::endl;

    // Compute R/r_c ratio from density contrast
    real R_over_rc = std::sqrt(density_contrast - 1.0);  // ≈ 3.62

    std::cout << "    R/r_c = sqrt(" << density_contrast << " - 1) = " << R_over_rc << std::endl;

    // Mass function: M(<R) = 4π ρ_c r_c³ [x - arctan(x)] where x = R/r_c
    auto compute_mass = [&](real R, real rc) -> real {
        real x = R / rc;
        return 4.0 * M_PI * rho_center * rc * rc * rc * (x - std::atan(x));
    };

    // With fixed R/r_c ratio, find R that gives M_target
    // r_c = R / R_over_rc
    // M = 4π ρ_c (R/R_over_rc)³ [R_over_rc - arctan(R_over_rc)]
    real x = R_over_rc;
    real mass_factor = 4.0 * M_PI * (x - std::atan(x)) / (x * x * x);  // ≈ 0.69 for x=3.62
    real R_estimate = std::pow(M_target / (mass_factor * rho_center), 1.0/3.0);
    real r_c = R_estimate / R_over_rc;

    // Iterate to refine (keeping R/r_c ratio fixed for K&I equilibrium)
    real R_cloud = R_estimate;
    for (int iter = 0; iter < 20; ++iter) {
        real M_current = compute_mass(R_cloud, r_c);
        real error = (M_current - M_target) / M_target;

        if (iter < 5) {
            std::cout << "  iter " << iter << ": R = " << R_cloud << " pc, r_c = " << r_c
                      << " pc, M = " << M_current << " M_sun, error = " << error * 100 << "%" << std::endl;
        }

        if (std::abs(error) < 0.001) {
            std::cout << "  Converged at iter " << iter << std::endl;
            break;
        }

        // Scale R to match mass (M ~ R³ approximately), keeping R/r_c fixed
        R_cloud *= std::pow(M_target / M_current, 1.0/3.0);
        r_c = R_cloud / R_over_rc;  // Maintain fixed ratio for K&I equilibrium
    }

    const real R_target = R_cloud;

    std::cout << std::endl;
    std::cout << "Modified isothermal sphere parameters:" << std::endl;
    std::cout << "  R_cloud = " << R_cloud << " pc" << std::endl;
    std::cout << "  r_core = " << r_c << " pc" << std::endl;
    std::cout << "  rho_center = " << rho_center << " M_sun/pc^3" << std::endl;
    std::cout << "  n_center = " << config.n_center << " cm^-3" << std::endl;

    // Build the profile
    BEProfile profile;
    const real dr = R_cloud / 500.0;  // 500 steps across cloud
    int n_steps = 500;
    profile.r.reserve(n_steps);
    profile.rho.reserve(n_steps);
    profile.M_enc.reserve(n_steps);

    real M_enc = 0.0;
    for (int i = 0; i < n_steps; ++i) {
        real r = (i + 1) * dr;
        real rho = rho_center / (1.0 + (r/r_c)*(r/r_c));
        real dM = 4.0 * M_PI * r * r * rho * dr;
        M_enc += dM;

        profile.r.push_back(r);
        profile.rho.push_back(rho);
        profile.M_enc.push_back(M_enc);
    }
    profile.R_cloud = R_target;
    profile.M_cloud = M_enc;

    real rho_edge = rho_center / (1.0 + (R_cloud/r_c)*(R_cloud/r_c));
    std::cout << "  rho_edge = " << rho_edge << " M_sun/pc^3" << std::endl;
    std::cout << "  n_edge = " << sample::code_density_to_n(rho_edge) << " cm^-3" << std::endl;
    std::cout << "  Density contrast = " << rho_center / rho_edge << std::endl;
    std::cout << "  M_enclosed at R = " << M_enc << " M_sun" << std::endl;

    std::cout << "Profile integration complete:" << std::endl;
    std::cout << "  R_BE = " << profile.R_cloud << " pc" << std::endl;
    std::cout << "  M_BE = " << profile.M_cloud << " M_sun" << std::endl;
    std::cout << "  rho_center = " << profile.rho.front() << " M_sun/pc^3" << std::endl;
    std::cout << "  rho_edge = " << profile.rho.back() << " M_sun/pc^3" << std::endl;
    std::cout << "  n_center = " << sample::code_density_to_n(profile.rho.front()) << " cm^-3" << std::endl;
    std::cout << "  n_edge = " << sample::code_density_to_n(profile.rho.back()) << " cm^-3" << std::endl;
    std::cout << "  Profile points: " << profile.r.size() << std::endl;
    std::cout << std::endl;

    // Use actual cloud parameters from integration
    const real R_actual = profile.R_cloud;
    const real M_actual = profile.M_cloud;

    // ========================================================================
    // CREATE PARTICLES FOLLOWING BE DENSITY PROFILE
    // ========================================================================

    const int N_particles = config.N_particles;

    std::cout << "Creating particles:" << std::endl;
    std::cout << "  N_target = " << N_particles << std::endl;

    // Internal energy: u = c_s^2 / (gamma - 1)
    const real gamma_eff = std::max(gamma - 1.0, 1.0e-4);
    const real u_cloud = c_s_sq / gamma_eff;

    // ========================================================================
    // HCP LATTICE-BASED PARTICLE PLACEMENT
    // ========================================================================
    // 1. Generate HCP lattice in unit sphere
    // 2. Stretch radially to match cumulative mass profile
    // This gives much better neighbor uniformity than random placement!

    std::cout << "Generating HCP lattice for " << N_particles << " particles..." << std::endl;

    // Generate HCP lattice in unit sphere
    std::vector<vec_t> lattice_points = generate_hcp_sphere(N_particles);

    // Stretch lattice to match BE density profile (radial mass mapping)
    stretch_to_density_profile(lattice_points, profile);

    // IMPORTANT: Calculate particle mass from ACTUAL number of lattice points
    // to ensure total mass matches M_cloud exactly
    const int N_actual = static_cast<int>(lattice_points.size());
    const real particle_mass = M_actual / N_actual;

    std::cout << "  Lattice points after stretching: " << N_actual << std::endl;
    std::cout << "  Particle mass = " << particle_mass << " M_sun (adjusted for actual N)" << std::endl;

    std::vector<SPHParticle> particles;
    particles.reserve(N_actual);

    int particle_id = 0;

    // Convert to SPH particles
    const int N_neighbor = m_param->physics.neighbor_number;
    constexpr real A_vol = 4.0 * M_PI / 3.0;

    for (const auto& pos : lattice_points) {
        real r = std::abs(pos);
        if (r > R_actual) continue;  // Skip any points outside cloud

        // Get local density from BE profile
        real rho_local = interpolate_density(profile, r);

        // Smoothing length from local density
        real h_local = std::pow(N_neighbor * particle_mass / (rho_local * A_vol), 1.0 / 3.0);

        SPHParticle p;
        p.pos = pos;
#if DIM == 3
        p.vel = vec_t{-config.v_rel, 0.0, 0.0};
#elif DIM == 2
        p.vel = vec_t{-config.v_rel, 0.0};
#else
        p.vel = vec_t{-config.v_rel};
#endif
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
    // CREATE GHOST ENVELOPE FOR PRESSURE CONFINEMENT
    // ========================================================================
    //
    // IMPORTANT: Envelope density must match K&I 2000 external pressure!
    // P_ext = rho_ext * c_s^2  =>  rho_ext = P_ext / c_s^2
    // At T=10K: n_ext = P_ext/(k_B*T) = 710/10 = 71 cm^-3
    // ========================================================================

    if (config.use_envelope) {
        // Envelope density from K&I 2000 external pressure (NOT cloud edge!)
        real n_ext = config.P_ext / config.T_cloud;  // cm^-3
        real rho_ext = sample::n_to_code_density(n_ext);

        std::cout << "Envelope (K&I 2000 external pressure):" << std::endl;
        std::cout << "  P_ext/k_B = " << config.P_ext << " K cm^-3" << std::endl;
        std::cout << "  T = " << config.T_cloud << " K" << std::endl;
        std::cout << "  n_ext = " << n_ext << " cm^-3 (from P_ext/T)" << std::endl;
        std::cout << "  rho_ext = " << rho_ext << " M_sun/pc^3" << std::endl;
        std::cout << "  Cloud edge: n = " << sample::code_density_to_n(profile.rho.back()) << " cm^-3" << std::endl;

        GhostEnvelopeConfig env_config;
        env_config.R_cloud = R_actual;
        env_config.rho_edge = rho_ext;  // Use K&I external density!
        env_config.u_envelope = u_cloud;
        env_config.particle_mass = particle_mass;
        env_config.N_neighbor = m_param->physics.neighbor_number;
        env_config.num_layers = config.envelope_layers;
        env_config.gamma = m_param->physics.gamma;  // CRITICAL: must match simulation gamma!

        auto envelope = GhostEnvelopeGenerator::generate(env_config);

        for (auto& p : envelope) {
            p.id = particle_id++;
#if DIM == 3
            p.vel = vec_t{-config.v_rel, 0.0, 0.0};
#elif DIM == 2
            p.vel = vec_t{-config.v_rel, 0.0};
#else
            p.vel = vec_t{-config.v_rel};
#endif
        }

        particles.insert(particles.end(), envelope.begin(), envelope.end());

        GhostEnvelopeGenerator::print_summary(env_config, envelope.size());
    }

    // ========================================================================
    // STORE PARAMETERS
    // ========================================================================

    m_sample_parameters["BH_position_x"] = config.b_impact;
    m_sample_parameters["BH_position_y"] = real(0.0);
    m_sample_parameters["BH_position_z"] = real(0.0);
    m_sample_parameters["BH_mass"] = config.M_BH;
    m_sample_parameters["BH_softening"] = config.epsilon_BH;
    m_sample_parameters["T_cloud_K"] = config.T_cloud;
    m_sample_parameters["c_s_kms"] = c_s;
    m_sample_parameters["R_cloud"] = R_actual;
    m_sample_parameters["M_cloud"] = M_actual;
    m_sample_parameters["P_ext_K_cm3"] = config.P_ext;
    m_sample_parameters["rho_center_code"] = profile.rho.front();
    m_sample_parameters["density_to_n"] = density_to_n_factor;

    // Parameters for isothermal relaxation (uses modified isothermal sphere profile)
    m_sample_parameters["T_cloud"] = config.T_cloud;
    m_sample_parameters["r_c_code"] = r_c;
    m_sample_parameters["R_cloud_code"] = R_actual;
    m_sample_parameters["P_ext"] = config.P_ext;

    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());

    // ========================================================================
    // PRINT SUMMARY
    // ========================================================================

    real t_sound = R_actual / c_s;
    real t_ff = std::sqrt(3.0 * M_PI / (32.0 * G * profile.rho.front()));
    real r_tidal = R_actual * std::pow(config.M_BH / M_actual, 1.0 / 3.0);

    std::cout << "==========================================================" << std::endl;
    std::cout << "  Isothermal Bonnor-Ebert IC Complete" << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << "Cloud:" << std::endl;
    std::cout << "  R = " << R_actual << " pc" << std::endl;
    std::cout << "  M = " << M_actual << " M_sun" << std::endl;
    std::cout << "  T = " << config.T_cloud << " K (isothermal)" << std::endl;
    std::cout << "  c_s = " << c_s << " km/s" << std::endl;
    std::cout << "  n_center = " << sample::code_density_to_n(profile.rho.front()) << " cm^-3" << std::endl;
    std::cout << "  n_edge = " << sample::code_density_to_n(profile.rho.back()) << " cm^-3" << std::endl;
    std::cout << std::endl;
    std::cout << "Particles: " << particles.size() << std::endl;
    std::cout << std::endl;
    if (config.M_BH > 0) {
        std::cout << "IMBH:" << std::endl;
        std::cout << "  M_BH = " << config.M_BH << " M_sun" << std::endl;
        std::cout << "  Position: (" << config.b_impact << ", 0, 0) pc" << std::endl;
        std::cout << "  r_tidal ~ " << r_tidal << " pc" << std::endl;
        std::cout << std::endl;
    }
    std::cout << "Timescales:" << std::endl;
    std::cout << "  t_sound = " << t_sound << " Myr" << std::endl;
    std::cout << "  t_ff = " << t_ff << " Myr" << std::endl;
    std::cout << "==========================================================" << std::endl;
#endif // DIM == 3
}

} // namespace sph
