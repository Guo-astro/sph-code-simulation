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
    // CREATE PARTICLES - SELECT METHOD
    // ========================================================================
    // Method 0: HEALPIX (default) - HEALPix equal-area + density weighting
    //           Best IC, minimal relaxation needed (~500 steps)
    // Method 1: SHELL - Concentric shells with Fibonacci spiral
    //           Good IC, minimal relaxation needed (~1000 steps)
    // Method 2: RANDOM - Random + cumulative mass mapping
    //           Requires longer relaxation (~10000+ steps)
    // ========================================================================

    // IC method: "healpix" (default), "shell", or "random"
    std::string ic_method = "healpix";  // Default to best method
    if (m_sample_parameters.count("ic_method")) {
        ic_method = boost::any_cast<std::string>(m_sample_parameters["ic_method"]);
    } else if (m_sample_parameters.count("use_random_ic")) {
        // Backward compatibility
        bool use_random = boost::any_cast<bool>(m_sample_parameters["use_random_ic"]);
        ic_method = use_random ? "random" : "healpix";
    }
    
    bool use_healpix_method = (ic_method == "healpix");
    bool use_shell_method = (ic_method == "shell");

    const int N_target = N * N * N;
    std::vector<vec_t> mapped_positions;
    std::vector<real> local_densities;

    if (use_healpix_method) {
        // ====================================================================
        // METHOD 0: HEALPIX + DENSITY WEIGHTING (BEST)
        // ====================================================================
        // HEALPix provides Hierarchical Equal Area isoLatitude Pixelization
        // Key advantage: equal-area pixels give uniform angular sampling
        // Combined with density-weighted radial shells for optimal IC
        // ====================================================================
        
        std::cout << "\nCreating particles (HEALPix + Density Weighting):" << std::endl;
        std::cout << "  N_target = " << N_target << std::endl;
        std::cout << "  This method gives fastest relaxation convergence" << std::endl;
        
        // HEALPix implementation: base resolution N_side
        // Total pixels = 12 * N_side^2
        // Start with N_side = 2^k for hierarchical structure
        
        // Lambda to compute HEALPix pixel center (ring scheme)
        // For pixel p in [0, 12*nside^2 - 1]:
        auto healpix_pix2vec = [](int nside, int pix) -> vec_t {
            int npix = 12 * nside * nside;
            int ncap = 2 * nside * (nside - 1);  // Pixels in north polar cap
            
            real z, phi;
            
            if (pix < ncap) {
                // North polar cap
                int ph = (pix + 1) / 2;
                int ring = static_cast<int>(0.5 * (1 + std::sqrt(1 + 2*ph)));
                int iphi = (pix + 1) - 2 * ring * (ring - 1);
                z = 1.0 - (ring * ring) / (3.0 * nside * nside);
                phi = (iphi - 0.5) * M_PI / (2.0 * ring);
            } else if (pix < npix - ncap) {
                // Equatorial region
                int ip = pix - ncap;
                int ring = ip / (4 * nside) + nside;
                int iphi = ip % (4 * nside) + 1;
                int fodd = ((ring + nside) % 2 == 0) ? 1 : 0;
                z = (2 * nside - ring) / (1.5 * nside);
                phi = (iphi - 0.5 * (1 + fodd)) * M_PI / (2.0 * nside);
            } else {
                // South polar cap
                int ip = npix - pix;
                int ph = (ip + 1) / 2;
                int ring = static_cast<int>(0.5 * (1 + std::sqrt(2*ph - 1)));
                int iphi = 4 * ring + 1 - (ip - 2 * ring * (ring - 1));
                z = -1.0 + (ring * ring) / (3.0 * nside * nside);
                phi = (iphi - 0.5) * M_PI / (2.0 * ring);
            }
            
            real sin_theta = std::sqrt(std::max(0.0, 1.0 - z*z));
            return vec_t(sin_theta * std::cos(phi), sin_theta * std::sin(phi), z);
        };
        
        // Determine N_side to approximate N_target total particles
        // We'll use multiple shells, each with 12*nside^2 particles
        // Total ~ N_shells * 12 * nside^2 ~ N_target
        // Balance: more shells = better radial resolution, higher nside = better angular
        
        // Estimate: N_shells ~ N_target^(1/3), N_pix_per_shell ~ N_target^(2/3)
        // N_side ~ sqrt(N_pix_per_shell / 12) ~ (N_target^(2/3) / 12)^(1/2)
        
        int N_shells = std::max(10, static_cast<int>(std::pow(N_target, 1.0/3.0)));
        int nside = std::max(2, static_cast<int>(std::sqrt(N_target / (12.0 * N_shells))));
        // Round nside to power of 2 for proper HEALPix structure
        int nside_pow2 = 1;
        while (nside_pow2 * 2 <= nside) nside_pow2 *= 2;
        nside = nside_pow2;
        
        int npix = 12 * nside * nside;  // Pixels per shell
        
        // Recompute N_shells to hit N_target
        N_shells = std::max(1, N_target / npix);
        
        std::cout << "  HEALPix N_side = " << nside << " (" << npix << " pixels/shell)" << std::endl;
        std::cout << "  N_shells = " << N_shells << std::endl;
        std::cout << "  Expected particles = " << N_shells * npix << std::endl;
        
        mapped_positions.reserve(N_shells * npix);
        local_densities.reserve(N_shells * npix);
        
        // Use cumulative mass profile to determine shell radii
        // Density weighting: inner dense shells have more particles per unit volume
        // We achieve this by using equal-mass shells (same as Fibonacci method)
        
        real dM_shell = M_total / N_shells;
        
        for (int s = 0; s < N_shells; ++s) {
            // Target mass for this shell
            real M_target = (s + 0.5) * dM_shell;
            
            // Find radius where M_enc = M_target
            real r_shell = 0.0;
            for (int i = 1; i <= i_s; ++i) {
                if (M_enc_arr[i] >= M_target) {
                    real f = (M_target - M_enc_arr[i-1]) / (M_enc_arr[i] - M_enc_arr[i-1] + 1e-20);
                    r_shell = r_arr[i-1] + f * (r_arr[i] - r_arr[i-1]);
                    break;
                }
            }
            if (r_shell < 1e-10) r_shell = r_arr[1] * (s + 1.0) / N_shells;
            if (r_shell > R_cloud * 0.99) r_shell = R_cloud * 0.99;
            
            // Get local density at this radius
            real rho_local = rho_c;
            for (int i = 1; i <= i_s; ++i) {
                if (r_shell <= r_arr[i]) {
                    real f = (r_shell - r_arr[i-1]) / (r_arr[i] - r_arr[i-1] + 1e-20);
                    rho_local = rho_arr[i-1] + f * (rho_arr[i] - rho_arr[i-1]);
                    break;
                }
            }
            if (rho_local < rho_edge) rho_local = rho_edge;
            
            // Rotate each shell by golden angle to avoid radial alignment
            real shell_rotation = s * 2.399963229728653;  // Golden angle in radians
            
            // Place HEALPix pixels on this shell
            for (int p = 0; p < npix; ++p) {
                vec_t dir = healpix_pix2vec(nside, p);
                
                // Apply shell rotation around z-axis to avoid radial structure
                real cos_rot = std::cos(shell_rotation);
                real sin_rot = std::sin(shell_rotation);
                real x_rot = dir[0] * cos_rot - dir[1] * sin_rot;
                real y_rot = dir[0] * sin_rot + dir[1] * cos_rot;
                
                vec_t pos(x_rot * r_shell, y_rot * r_shell, dir[2] * r_shell);
                mapped_positions.push_back(pos);
                local_densities.push_back(rho_local);
            }
        }
        
        std::cout << "  Created " << N_shells << " shells × " << npix << " pixels" << std::endl;
        std::cout << "  Total HEALPix particles: " << mapped_positions.size() << std::endl;
        std::cout << "  (HEALPix gives uniform angular distribution + optimal radial spacing)" << std::endl;
        
    } else if (use_shell_method) {
        // ====================================================================
        // METHOD 1: CONCENTRIC SHELLS WITH FIBONACCI SPIRAL (BEST)
        // ====================================================================
        // - Place particles on spherical shells at radii determined by M(r)
        // - Use Fibonacci spiral for uniform angular distribution on each shell
        // - Number of shells and particles per shell computed to hit N_target
        // ====================================================================
        
        std::cout << "\nCreating particles (Concentric Shells + Fibonacci):" << std::endl;
        std::cout << "  N_target = " << N_target << std::endl;

        // Fibonacci spiral constants
        const real golden_ratio = (1.0 + std::sqrt(5.0)) / 2.0;
        const real golden_angle = 2.0 * M_PI / (golden_ratio * golden_ratio);

        // For N particles in 3D sphere: N ~ (4/3)*pi*(R/h)^3
        // So h ~ R * (4*pi/(3*N))^(1/3)
        // Number of shells ~ R / (spacing * h) ~ N^(1/3)
        
        const int N_neighbor = m_param->physics.neighbor_number;
        constexpr real A_vol = 4.0 * M_PI / 3.0;
        
        // Estimate number of shells: roughly cube root of N for 3D
        int N_shells = static_cast<int>(std::pow(N_target, 1.0/3.0) * 1.5);
        N_shells = std::max(N_shells, 10);  // Minimum shells
        
        // Build shells from inside out using cumulative mass
        mapped_positions.reserve(N_target);
        local_densities.reserve(N_target);
        
        // Use cumulative mass profile to determine shell radii
        // Each shell contains equal mass dM = M_total / N_shells
        real dM_shell = M_total / N_shells;
        
        int total_placed = 0;
        int shell_num = 0;
        real M_placed = 0.0;
        
        for (int s = 0; s < N_shells && total_placed < N_target; ++s) {
            // Target mass for this shell
            real M_target = (s + 0.5) * dM_shell;
            
            // Find radius where M_enc = M_target
            real r_shell = 0.0;
            for (int i = 1; i <= i_s; ++i) {
                if (M_enc_arr[i] >= M_target) {
                    real f = (M_target - M_enc_arr[i-1]) / (M_enc_arr[i] - M_enc_arr[i-1] + 1e-20);
                    r_shell = r_arr[i-1] + f * (r_arr[i] - r_arr[i-1]);
                    break;
                }
            }
            if (r_shell < 1e-10) r_shell = r_arr[1] * (s + 1.0) / N_shells;  // Fallback for center
            if (r_shell > R_cloud * 0.99) r_shell = R_cloud * 0.99;
            
            // Get local density at this radius
            real rho_local = rho_c;
            for (int i = 1; i <= i_s; ++i) {
                if (r_shell <= r_arr[i]) {
                    real f = (r_shell - r_arr[i-1]) / (r_arr[i] - r_arr[i-1] + 1e-20);
                    rho_local = rho_arr[i-1] + f * (rho_arr[i] - rho_arr[i-1]);
                    break;
                }
            }
            if (rho_local < rho_edge) rho_local = rho_edge;
            
            // Particles on this shell: proportional to shell area * local density
            // Total particles = N_target, distribute by: N_shell / N_total ~ 4*pi*r^2 * rho * dr
            // Simplified: N_shell ~ (N_target / N_shells) * (rho_local / mean_rho) for equal mass shells
            int N_shell = (N_target - total_placed) / (N_shells - s);
            N_shell = std::max(N_shell, 4);  // Minimum for a shell
            
            // Don't exceed remaining
            if (total_placed + N_shell > N_target) {
                N_shell = N_target - total_placed;
            }
            
            // Place particles on shell using Fibonacci spiral
            for (int i = 0; i < N_shell; ++i) {
                // Fibonacci spiral: uniform distribution on sphere
                real y = 1.0 - (2.0 * i + 1.0) / N_shell;  // y from +1 to -1
                real radius_at_y = std::sqrt(std::max(0.0, 1.0 - y * y));
                real theta = golden_angle * i;
                
                real x = radius_at_y * std::cos(theta);
                real z = radius_at_y * std::sin(theta);
                
                // Add small random perturbation to avoid exact symmetry
                std::mt19937 rng(42 + s * 1000 + i);
                std::uniform_real_distribution<real> perturb(-0.01, 0.01);
                x += perturb(rng) * r_shell / R_cloud;
                y += perturb(rng) * r_shell / R_cloud;
                z += perturb(rng) * r_shell / R_cloud;
                
                // Renormalize to shell radius
                real norm = std::sqrt(x*x + y*y + z*z);
                if (norm > 1e-10) {
                    x *= r_shell / norm;
                    y *= r_shell / norm;
                    z *= r_shell / norm;
                }
                
                vec_t pos(x, y, z);
                mapped_positions.push_back(pos);
                local_densities.push_back(rho_local);
            }
            
            total_placed += N_shell;
            shell_num++;
        }
        
        std::cout << "  Created " << shell_num << " shells" << std::endl;
        std::cout << "  Placed " << mapped_positions.size() << " particles" << std::endl;
        std::cout << "  (Shell method gives ~10x faster relaxation)" << std::endl;
        
    } else {
        // ====================================================================
        // METHOD 2: RANDOM + CUMULATIVE MASS MAPPING (Original)
        // ====================================================================
        
        std::cout << "\nCreating particles (Random + mass mapping):" << std::endl;
        std::cout << "  N_target = " << N_target << std::endl;
        std::cout << "  WARNING: This method requires longer relaxation" << std::endl;

        // Generate random points uniformly in a unit sphere
        std::vector<vec_t> random_points;
        random_points.reserve(N_target);

        std::mt19937 rng(42);
        std::uniform_real_distribution<real> uniform(-1.0, 1.0);

        while (static_cast<int>(random_points.size()) < N_target) {
            real x = uniform(rng);
            real y = uniform(rng);
            real z = uniform(rng);
            real r_sq = x*x + y*y + z*z;

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

        // Map random points to BE density profile
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
        env_config.gamma = gamma;  // CRITICAL: must match simulation gamma!

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
