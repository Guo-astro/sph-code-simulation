/**
 * @file koyama_inutsuka_relaxation.cpp
 * @brief Implementation of K&I 2000 Bonnor-Ebert analytical relaxation
 */

#include "relaxation/koyama_inutsuka_relaxation.hpp"
#include "exception.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sph {

KoyamaInutsukaRelaxation::KoyamaInutsukaRelaxation()
    : m_initialized(false)
    , m_R_actual(0.0)
    , m_M_actual(0.0)
    , m_ki_cooling()
{
}

void KoyamaInutsukaRelaxation::initialize(const KIRelaxationParams& params)
{
    m_params = params;

    std::cout << "KoyamaInutsukaRelaxation: Initializing with:" << std::endl;
    std::cout << "  Target R_cloud = " << params.R_cloud << " [code]" << std::endl;
    std::cout << "  Target M_cloud = " << params.M_cloud << " [code]" << std::endl;
    std::cout << "  P_ext = " << params.P_ext << " [code]" << std::endl;
    std::cout << "  rho_center = " << params.rho_center << " [code]" << std::endl;
    std::cout << "  N_H = " << params.N_H << " [cm^-2]" << std::endl;
    std::cout << "  G = " << params.G << std::endl;
    std::cout << "  density_to_n = " << params.density_to_n << std::endl;

    // Compute equilibrium profile by integrating hydrostatic ODE
    compute_equilibrium_profile();

    m_initialized = true;

    std::cout << "KoyamaInutsukaRelaxation: Profile computed:" << std::endl;
    std::cout << "  Actual R_cloud = " << m_R_actual << " [code]" << std::endl;
    std::cout << "  Actual M_cloud = " << m_M_actual << " [code]" << std::endl;
    std::cout << "  n_H(center) = " << m_params.rho_center * m_params.density_to_n << " [cm^-3]" << std::endl;
    std::cout << "  T_eq(center) = " << m_T_table[0] << " [K]" << std::endl;
    std::cout << "  Profile points: " << m_r_table.size() << std::endl;
}

real KoyamaInutsukaRelaxation::get_T_eq_from_rho(real rho) const
{
    // Convert code density to n_H [cm^-3]
    real n_H = rho * m_params.density_to_n;
    n_H = std::max(n_H, 0.01);  // Clamp to table range
    n_H = std::min(n_H, 1.0e4);

    // Get equilibrium temperature from K&I 2000 tables
    return m_ki_cooling.equilibrium_temperature(n_H, m_params.N_H);
}

real KoyamaInutsukaRelaxation::get_c_eff_squared(real rho) const
{
    // c_eff² = dP/dρ = (k_B T / μ m_H) × (1 + d ln T / d ln n)
    //
    // For ideal gas: P = n k_B T = (ρ / m_n) k_B T
    // dP/dρ = (k_B / m_n) × [T + n dT/dn]
    //       = (k_B T / m_n) × [1 + (n/T) dT/dn]
    //       = (k_B T / m_n) × [1 + d ln T / d ln n]

    real n_H = rho * m_params.density_to_n;
    n_H = std::max(n_H, 0.01);
    n_H = std::min(n_H, 1.0e4);

    real T_eq = m_ki_cooling.equilibrium_temperature(n_H, m_params.N_H);

    // Numerical derivative d ln T / d ln n using finite difference
    const real eps = 0.02;  // 2% for numerical stability
    real n_lo = n_H * (1.0 - eps);
    real n_hi = n_H * (1.0 + eps);
    n_lo = std::max(n_lo, 0.01);
    n_hi = std::min(n_hi, 1.0e4);

    real T_lo = m_ki_cooling.equilibrium_temperature(n_lo, m_params.N_H);
    real T_hi = m_ki_cooling.equilibrium_temperature(n_hi, m_params.N_H);

    real dlnT_dlnn = 0.0;
    if (T_lo > 0.0 && T_hi > 0.0) {
        dlnT_dlnn = (std::log(T_hi) - std::log(T_lo)) / (std::log(n_hi) - std::log(n_lo));
    }

    // c_eff² in CGS
    real c_eff_sq_cgs = (k_B * T_eq / m_n) * (1.0 + dlnT_dlnn);

    // Convert to code units: [cm/s]² → [code velocity]²
    // Code velocity = 1 km/s = 1e5 cm/s
    // So c_eff²[code] = c_eff²[CGS] / (1e5)² = c_eff²[CGS] / 1e10
    real c_eff_sq_code = c_eff_sq_cgs / 1.0e10;

    return c_eff_sq_code;
}

void KoyamaInutsukaRelaxation::compute_equilibrium_profile()
{
    // Integrate the hydrostatic equilibrium ODE:
    //   dρ/dr = -ρ G M(r) / (r² c_eff²)
    //   dM/dr = 4 π r² ρ
    //
    // Stopping condition: P(ρ) ≤ P_ext (pressure truncation)

    const int n_points_max = 2000;
    const real dr_factor = 0.001;  // dr = dr_factor × R_target

    m_r_table.clear();
    m_rho_table.clear();
    m_P_table.clear();
    m_M_table.clear();
    m_T_table.clear();
    m_drho_dr_table.clear();

    m_r_table.reserve(n_points_max);
    m_rho_table.reserve(n_points_max);
    m_P_table.reserve(n_points_max);
    m_M_table.reserve(n_points_max);
    m_T_table.reserve(n_points_max);
    m_drho_dr_table.reserve(n_points_max);

    // Initial conditions at center
    real r = 1.0e-6 * m_params.R_cloud;  // Start at small r, not zero
    real rho = m_params.rho_center;
    real M = (4.0 / 3.0) * M_PI * rho * r * r * r;

    // Get initial pressure
    real n_H_center = rho * m_params.density_to_n;
    real T_center = get_T_eq_from_rho(rho);
    real P_center_cgs = n_H_center * k_B * T_center;  // [dyn/cm²]

    // Convert P_ext to CGS for comparison
    // P_ext is in code units; assume P_ext = (P/k_B) in K cm^-3
    // P[dyn/cm²] = (P/k_B) × k_B
    real P_ext_cgs = m_params.P_ext * k_B;

    std::cout << "  Integrating profile:" << std::endl;
    std::cout << "    P_center = " << P_center_cgs / k_B << " K cm^-3" << std::endl;
    std::cout << "    P_ext = " << P_ext_cgs / k_B << " K cm^-3" << std::endl;

    real dr = dr_factor * m_params.R_cloud;
    real r_max = 5.0 * m_params.R_cloud;  // Safety limit

    bool truncated = false;

    for (int i = 0; i < n_points_max && r < r_max; ++i) {
        // Store current state
        real n_H = rho * m_params.density_to_n;
        real T = get_T_eq_from_rho(rho);
        real P_cgs = n_H * k_B * T;  // [dyn/cm²]

        m_r_table.push_back(r);
        m_rho_table.push_back(rho);
        m_P_table.push_back(P_cgs / k_B);  // Store as P/k_B [K cm^-3]
        m_M_table.push_back(M);
        m_T_table.push_back(T);

        // Check truncation condition
        if (P_cgs <= P_ext_cgs) {
            truncated = true;
            m_R_actual = r;
            m_M_actual = M;
            std::cout << "    Truncated at r = " << r << ", P/k = " << P_cgs/k_B << std::endl;
            break;
        }

        // Get effective sound speed
        real c_eff_sq = get_c_eff_squared(rho);

        // Ensure c_eff² is positive and physical
        if (c_eff_sq <= 0.0) {
            c_eff_sq = 1.0e-6;  // Minimum value
        }

        // Hydrostatic equilibrium: dρ/dr = -ρ G M / (r² c_eff²)
        real drho_dr = -rho * m_params.G * M / (r * r * c_eff_sq);
        m_drho_dr_table.push_back(drho_dr);

        // RK4 integration step
        // k1
        real k1_rho = drho_dr;
        real k1_M = 4.0 * M_PI * r * r * rho;

        // k2
        real r_mid = r + 0.5 * dr;
        real rho_mid = rho + 0.5 * dr * k1_rho;
        real M_mid = M + 0.5 * dr * k1_M;
        rho_mid = std::max(rho_mid, 1.0e-30);
        real c_eff_sq_mid = get_c_eff_squared(rho_mid);
        c_eff_sq_mid = std::max(c_eff_sq_mid, 1.0e-6);
        real k2_rho = -rho_mid * m_params.G * M_mid / (r_mid * r_mid * c_eff_sq_mid);
        real k2_M = 4.0 * M_PI * r_mid * r_mid * rho_mid;

        // k3
        rho_mid = rho + 0.5 * dr * k2_rho;
        M_mid = M + 0.5 * dr * k2_M;
        rho_mid = std::max(rho_mid, 1.0e-30);
        c_eff_sq_mid = get_c_eff_squared(rho_mid);
        c_eff_sq_mid = std::max(c_eff_sq_mid, 1.0e-6);
        real k3_rho = -rho_mid * m_params.G * M_mid / (r_mid * r_mid * c_eff_sq_mid);
        real k3_M = 4.0 * M_PI * r_mid * r_mid * rho_mid;

        // k4
        real r_end = r + dr;
        real rho_end = rho + dr * k3_rho;
        real M_end = M + dr * k3_M;
        rho_end = std::max(rho_end, 1.0e-30);
        real c_eff_sq_end = get_c_eff_squared(rho_end);
        c_eff_sq_end = std::max(c_eff_sq_end, 1.0e-6);
        real k4_rho = -rho_end * m_params.G * M_end / (r_end * r_end * c_eff_sq_end);
        real k4_M = 4.0 * M_PI * r_end * r_end * rho_end;

        // Update
        rho = rho + (dr / 6.0) * (k1_rho + 2.0*k2_rho + 2.0*k3_rho + k4_rho);
        M = M + (dr / 6.0) * (k1_M + 2.0*k2_M + 2.0*k3_M + k4_M);
        r = r_end;

        // Ensure physical values
        rho = std::max(rho, 1.0e-30);
    }

    // Add final gradient entry if we truncated
    if (truncated && m_drho_dr_table.size() < m_rho_table.size()) {
        m_drho_dr_table.push_back(m_drho_dr_table.back());
    }

    if (!truncated) {
        std::cout << "  WARNING: Profile did not truncate at P_ext!" << std::endl;
        std::cout << "           Final r = " << r << ", final P/k = " << m_P_table.back() << std::endl;
        m_R_actual = m_r_table.back();
        m_M_actual = m_M_table.back();
    }

    std::cout << "  Profile integration complete:" << std::endl;
    std::cout << "    Points: " << m_r_table.size() << std::endl;
    std::cout << "    R_cloud: " << m_R_actual << " [code]" << std::endl;
    std::cout << "    M_cloud: " << m_M_actual << " [code]" << std::endl;
    std::cout << "    rho_center: " << m_rho_table[0] << " [code]" << std::endl;
    std::cout << "    rho_edge: " << m_rho_table.back() << " [code]" << std::endl;
}

real KoyamaInutsukaRelaxation::interpolate(const std::vector<real>& table, real r) const
{
    if (m_r_table.empty() || table.empty()) {
        return 0.0;
    }

    // Clamp r to table range
    if (r <= m_r_table.front()) {
        return table.front();
    }
    if (r >= m_r_table.back()) {
        return table.back();
    }

    // Binary search for interval
    auto it = std::lower_bound(m_r_table.begin(), m_r_table.end(), r);
    size_t i = std::distance(m_r_table.begin(), it);
    if (i == 0) i = 1;

    // Linear interpolation
    real r0 = m_r_table[i-1];
    real r1 = m_r_table[i];
    real t = (r - r0) / (r1 - r0);

    return table[i-1] + t * (table[i] - table[i-1]);
}

real KoyamaInutsukaRelaxation::get_rho_eq(real r) const
{
    return interpolate(m_rho_table, r);
}

real KoyamaInutsukaRelaxation::get_P_eq(real r) const
{
    return interpolate(m_P_table, r);
}

real KoyamaInutsukaRelaxation::get_T_eq(real r) const
{
    return interpolate(m_T_table, r);
}

real KoyamaInutsukaRelaxation::get_drho_dr(real r) const
{
    return interpolate(m_drho_dr_table, r);
}

vec_t KoyamaInutsukaRelaxation::compute_relaxation_force(const SPHParticle& p) const
{
    if (!m_initialized) {
        THROW_ERROR("KoyamaInutsukaRelaxation not initialized");
    }

    vec_t force(0.0);

    // Spherical radius
    const real r = std::abs(p.pos);

    if (r < 1e-12) {
        return force;  // No force at center
    }

    if (r > m_R_actual) {
        return force;  // Outside cloud
    }

    // Get equilibrium values at this radius
    real rho_eq = get_rho_eq(r);
    real drho_dr = get_drho_dr(r);
    real c_eff_sq = get_c_eff_squared(rho_eq);

    if (rho_eq <= 1.0e-30 || std::abs(drho_dr) < 1.0e-30) {
        return force;
    }

    // Analytical pressure gradient acceleration:
    // a_r = -(1/ρ) dP/dr
    //
    // For barotropic EOS: dP/dr = (dP/dρ) (dρ/dr) = c_eff² dρ/dr
    // So: a_r = -c_eff² (1/ρ) dρ/dr
    //
    // This is the pressure support that balances gravity in equilibrium.
    // By subtracting this from SPH acceleration, we drive particles
    // toward positions where SPH pressure gradient matches analytical.

    real a_r = -c_eff_sq * drho_dr / rho_eq;

    // Convert to Cartesian direction (radial outward)
    const real r_inv = 1.0 / r;
    force[0] = a_r * p.pos[0] * r_inv;
    force[1] = a_r * p.pos[1] * r_inv;
#if DIM == 3
    force[2] = a_r * p.pos[2] * r_inv;
#endif

    return force;
}

void KoyamaInutsukaRelaxation::apply_relaxation(std::shared_ptr<Simulation> sim, real damping_factor)
{
    if (!m_initialized) {
        THROW_ERROR("KoyamaInutsukaRelaxation not initialized");
    }

    auto& particles = sim->get_particles();
    const int num_p = sim->get_particle_num();

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < num_p; ++i) {
        // Compute analytical pressure gradient acceleration
        vec_t relax_acc = compute_relaxation_force(particles[i]);

        // SUBTRACT analytical pressure gradient from SPH acceleration
        // Net acceleration = SPH forces - analytical equilibrium forces
        // This drives particles toward equilibrium positions
        particles[i].acc[0] -= relax_acc[0];
        particles[i].acc[1] -= relax_acc[1];
#if DIM == 3
        particles[i].acc[2] -= relax_acc[2];
#endif
    }
}

} // namespace sph
