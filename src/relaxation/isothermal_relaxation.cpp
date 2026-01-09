/**
 * @file isothermal_relaxation.cpp
 * @brief Implementation of analytical relaxation for TRUE Bonnor-Ebert spheres
 *
 * Uses TRUE LANE-EMDEN relaxation:
 * - Solves isothermal Lane-Emden equation: d²ψ/dξ² + (2/ξ)dψ/dξ = exp(-ψ)
 * - Density profile: ρ(ξ) = ρ_c × exp(-ψ(ξ))
 * - Equilibrium force: a_eq = (c_s²/r_0) × (dψ/dξ) × r̂
 * - Velocities zeroed every step in the solver loop (quasi-static)
 * - Position update: Δx = ½at² (no velocity term)
 */

#include "relaxation/isothermal_relaxation.hpp"
#include "bhtree.hpp"
#include "kernel/kernel_function.hpp"
#include "periodic.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

namespace sph {

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

IsothermalRelaxation::IsothermalRelaxation()
    : m_initialized(false)
    , m_c_s_sq(0.0)
    , m_r_0(0.0)
    , m_rho_edge(0.0)
    , m_n_profile(1000)
{
}

void IsothermalRelaxation::initialize(const IsothermalRelaxationParams& params)
{
    m_params = params;

    std::cout << "=== Initializing TRUE Bonnor-Ebert Relaxation ===" << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "  T_cloud = " << params.T_cloud << " K" << std::endl;
    std::cout << "  rho_center = " << params.rho_center << " [code]" << std::endl;
    std::cout << "  r_0 = " << params.r_0 << " [code]" << std::endl;
    std::cout << "  xi_s = " << params.xi_s << std::endl;
    std::cout << "  R_cloud = " << params.R_cloud << " [code]" << std::endl;
    std::cout << "  P_ext = " << params.P_ext << " [code]" << std::endl;
    std::cout << "  G = " << params.G << " [code]" << std::endl;
    std::cout << "  mu = " << params.mu << std::endl;

    // Compute sound speed squared
    // c_s^2 = k_B T / (mu m_H) [CGS: cm^2/s^2]
    // Convert to code units: [km/s]^2
    const real m_n_cgs = params.mu * m_proton_cgs;
    const real c_s_sq_cgs = k_B_cgs * params.T_cloud / m_n_cgs;
    m_c_s_sq = c_s_sq_cgs / (kms_cgs * kms_cgs);

    std::cout << "  c_s = " << std::sqrt(m_c_s_sq) << " km/s" << std::endl;
    std::cout << "  c_s^2 = " << m_c_s_sq << " (km/s)^2" << std::endl;

    // Store scale length
    m_r_0 = params.r_0;

    // Solve Lane-Emden equation
    std::cout << "\nSolving isothermal Lane-Emden equation..." << std::endl;
    m_n_profile = 1000;
    solve_isothermal_lane_emden(params.xi_s * 1.1, m_n_profile, m_xi_arr, m_psi_arr, m_dpsi_arr);

    // Find psi at truncation radius
    real psi_s = interpolate_psi(params.xi_s);
    real dpsi_s = interpolate_dpsi(params.xi_s);

    std::cout << "  At xi_s = " << params.xi_s << ":" << std::endl;
    std::cout << "    psi = " << psi_s << std::endl;
    std::cout << "    dpsi/dxi = " << dpsi_s << std::endl;
    std::cout << "    rho_edge/rho_center = " << std::exp(-psi_s) << std::endl;

    // Compute edge density from Lane-Emden solution
    m_rho_edge = params.rho_center * std::exp(-psi_s);

    real n_center = params.rho_center * params.density_to_n;
    real n_edge = m_rho_edge * params.density_to_n;

    std::cout << "\nTRUE Bonnor-Ebert Profile:" << std::endl;
    std::cout << "  rho(xi) = rho_c * exp(-psi(xi))" << std::endl;
    std::cout << "  n_center = " << n_center << " cm^-3" << std::endl;
    std::cout << "  n_edge = " << n_edge << " cm^-3" << std::endl;
    std::cout << "  rho_edge/rho_center = " << m_rho_edge / params.rho_center << std::endl;

    // Check K&I equilibrium: n_edge should equal P_ext/T
    real n_ext_KI = params.P_ext / params.T_cloud;
    std::cout << "\nK&I 2000 equilibrium check:" << std::endl;
    std::cout << "  n_ext (K&I) = P_ext/T = " << n_ext_KI << " cm^-3" << std::endl;
    std::cout << "  n_edge (profile) = " << n_edge << " cm^-3" << std::endl;

    if (std::abs(n_edge - n_ext_KI) / n_ext_KI > 0.1) {
        std::cout << "  WARNING: Cloud edge density does not match K&I external!" << std::endl;
        std::cout << "  This may indicate pressure mismatch at cloud boundary." << std::endl;
    } else {
        std::cout << "  Cloud boundary matches K&I equilibrium (good)" << std::endl;
    }

    std::cout << "\nRelaxation method: TRUE LANE-EMDEN (subtract analytical force from SPH)" << std::endl;
    std::cout << "Equilibrium force: a_eq = (c_s^2/r_0) * (dpsi/dxi) * r_hat" << std::endl;
    std::cout << "Net acceleration: a_net = a_SPH - a_eq (drives to equilibrium)" << std::endl;
    std::cout << "Ghost particle treatment:" << std::endl;
    std::cout << "  For r > R_cloud: rho_eq = constant = rho_edge = " << m_rho_edge << std::endl;
    std::cout << "  This provides proper external pressure confinement" << std::endl;

    m_initialized = true;
    std::cout << "=== TRUE Bonnor-Ebert Relaxation Initialized ===" << std::endl;
}

real IsothermalRelaxation::interpolate_psi(real xi) const
{
    if (xi <= m_xi_arr[0]) {
        // Near origin: psi ~ xi^2/6
        return xi * xi / 6.0;
    }
    if (xi >= m_xi_arr[m_n_profile - 1]) {
        return m_psi_arr[m_n_profile - 1];
    }

    // Binary search for interval
    int lo = 0, hi = m_n_profile - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (m_xi_arr[mid] <= xi) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    // Linear interpolation
    real f = (xi - m_xi_arr[lo]) / (m_xi_arr[hi] - m_xi_arr[lo]);
    return m_psi_arr[lo] + f * (m_psi_arr[hi] - m_psi_arr[lo]);
}

real IsothermalRelaxation::interpolate_dpsi(real xi) const
{
    if (xi <= m_xi_arr[0]) {
        // Near origin: dpsi/dxi ~ xi/3
        return xi / 3.0;
    }
    if (xi >= m_xi_arr[m_n_profile - 1]) {
        return m_dpsi_arr[m_n_profile - 1];
    }

    // Binary search for interval
    int lo = 0, hi = m_n_profile - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (m_xi_arr[mid] <= xi) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    // Linear interpolation
    real f = (xi - m_xi_arr[lo]) / (m_xi_arr[hi] - m_xi_arr[lo]);
    return m_dpsi_arr[lo] + f * (m_dpsi_arr[hi] - m_dpsi_arr[lo]);
}

real IsothermalRelaxation::get_rho_eq(real r) const
{
    // For r > R_cloud: return constant rho_edge (external medium density)
    // This ensures ghost/envelope particles provide proper external pressure
    if (r > m_params.R_cloud) {
        return m_rho_edge;
    }

    // TRUE Bonnor-Ebert profile: rho(xi) = rho_c * exp(-psi(xi))
    // where xi = r / r_0
    real xi = r / m_r_0;
    real psi = interpolate_psi(xi);
    return m_params.rho_center * std::exp(-psi);
}

real IsothermalRelaxation::get_P_eq(real r) const
{
    // P = rho * c_s^2
    return get_rho_eq(r) * m_c_s_sq;
}

vec_t IsothermalRelaxation::compute_relaxation_force(const SPHParticle& p) const
{
    // TRUE Bonnor-Ebert equilibrium force
    // From hydrostatic equilibrium: dP/dr = -rho * g
    // For isothermal: P = rho * c_s^2, so c_s^2 * d(rho)/dr = -rho * g
    // Therefore: a_pressure = -(1/rho) * dP/dr = -c_s^2 * d(ln rho)/dr
    //
    // Since rho = rho_c * exp(-psi), we have ln(rho) = ln(rho_c) - psi
    // So d(ln rho)/dr = -dpsi/dr = -(dpsi/dxi) * (dxi/dr) = -(dpsi/dxi) / r_0
    //
    // Therefore: a_pressure = c_s^2 * (dpsi/dxi) / r_0   (pointing outward)
    //
    // This is the equilibrium pressure gradient force that balances gravity

    if (!m_initialized) {
        return vec_t(0.0);
    }

    real r = std::abs(p.pos);
    if (r < 1.0e-10) {
        return vec_t(0.0);
    }

    if (r > m_params.R_cloud * 1.1) {
        return vec_t(0.0);
    }

    // Compute xi = r / r_0
    real xi = r / m_r_0;

    // Get dpsi/dxi from Lane-Emden solution
    real dpsi_dxi = interpolate_dpsi(xi);

    // Equilibrium acceleration magnitude: a_eq = c_s^2 * (dpsi/dxi) / r_0
    // This points OUTWARD (pressure support against gravity)
    real a_r_mag = m_c_s_sq * dpsi_dxi / m_r_0;

    vec_t r_hat = p.pos / r;
    return r_hat * a_r_mag;
}

void IsothermalRelaxation::apply_relaxation(std::shared_ptr<Simulation> sim, real damping_factor)
{
    if (!m_initialized) {
        return;
    }

    auto& particles = sim->get_particles();
    const int num = sim->get_particle_num();

    // ================================================================
    // TRUE LANE-EMDEN RELAXATION (same approach as lane_emden_relaxation.cpp):
    // - SUBTRACT analytical pressure gradient from SPH acceleration
    // - Net acceleration = SPH forces - Lane-Emden equilibrium pressure gradient
    // - This drives particles toward equilibrium positions
    // - Velocities are zeroed in the solver loop (quasi-static)
    // ================================================================

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < num; ++i) {
        auto& p_i = particles[i];
        if (p_i.is_ghost) continue;

        // Compute relaxation acceleration (analytical pressure gradient from Lane-Emden)
        vec_t relax_acc = compute_relaxation_force(p_i);

        // SUBTRACT analytical pressure gradient from SPH acceleration
        // Net acceleration = SPH forces - Lane-Emden equilibrium pressure gradient
        // This drives particles toward equilibrium positions
        p_i.acc[0] -= relax_acc[0];
        p_i.acc[1] -= relax_acc[1];
#if DIM == 3
        p_i.acc[2] -= relax_acc[2];
#endif

        // NOTE: Velocities are zeroed in the solver loop, not here
        // No damping needed - quasi-static relaxation zeros velocities each step
    }
}

bool IsothermalRelaxation::update_profile_from_sph(std::shared_ptr<Simulation> sim)
{
    // For TRUE Lane-Emden relaxation, we do NOT update the analytical profile
    // The Lane-Emden solution is fixed and defines the target equilibrium
    // The particles should relax toward this fixed profile
    (void)sim;  // Unused
    return false;
}

} // namespace sph
