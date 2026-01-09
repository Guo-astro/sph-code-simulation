/**
 * SRMHD Fluid Force Calculation
 *
 * Special Relativistic Magnetohydrodynamics using Godunov SPH.
 * Based on SR-GSPH formulation (Kitajima et al.) + MHD extension.
 *
 * KEY: Uses volume-based formulation like SR-GSPH, NOT mass-based like GSPMHD!
 * - Kernel: sqrt(2) * h
 * - Force: V² * Ω * ∇W where V = ν/N
 * - Energy: power = dot(v_star, force)
 * - Normalize: divide by ν
 */

#include "srmhd/srmhd_fluid_force.hpp"
#include "srmhd/srmhd_riemann.hpp"
#include "srmhd/srmhd_primitive_recovery.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "periodic.hpp"
#include "openmp.hpp"
#include "kernel/kernel_function.hpp"
#include "bhtree.hpp"
#include "logger.hpp"
#include <array>
#include <cmath>
#include <algorithm>

#ifdef EXHAUSTIVE_SEARCH
#include "exhaustive_search.hpp"
#endif

namespace sph
{
namespace srmhd
{

void FluidForce::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::FluidForce::initialize(param);

    m_gamma = param->physics.gamma;
    m_c_speed = param->srgsph.c_speed;
    m_c_shock = (param->srgsph.c_shock > 0.0) ? param->srgsph.c_shock : 3.0;
    m_c_cd = (param->srgsph.c_cd > 0.0) ? param->srgsph.c_cd : 0.2;
    m_use_muscl = param->srgsph.is_2nd_order;
    m_use_powell = param->mhd.use_powell_correction;
    m_use_mhd = param->mhd.use_mhd;
    m_riemann_solver = param->srgsph.riemann_solver;

    std::cout << "[SRMHD] Special Relativistic MHD initialized (SR-GSPH style)" << std::endl;
    if (m_use_mhd) {
        std::cout << "  Mode: Full MHD" << std::endl;
        std::cout << "  Riemann solver: "
                  << (m_riemann_solver == RiemannSolverType::HLLC ? "HLLC" : "HLL")
                  << std::endl;
        std::cout << "  Powell div-B correction: " << (m_use_powell ? "yes" : "no") << std::endl;
    } else {
        std::cout << "  Mode: HYDRO-ONLY (using SR-GSPH exact Riemann solver)" << std::endl;
    }
    std::cout << "  2nd order MUSCL: " << (m_use_muscl ? "yes" : "no") << std::endl;
    std::cout << "  gamma: " << m_gamma << ", c: " << m_c_speed << std::endl;
}

real FluidForce::limiter(real dq1, real dq2) const
{
    const real dq1dq2 = dq1 * dq2;
    if (dq1dq2 <= 0) {
        return 0.0;
    } else {
        return 2.0 * dq1dq2 / (dq1 + dq2);
    }
}

void FluidForce::srmhd_riemann_solver(
    real rho_L, real Pt_L, real v_L, real c_L, real H_L, real gamma_L,
    real rho_R, real Pt_R, real v_R, real c_R, real H_R, real gamma_R,
    real& Pt_star, real& v_star) const
{
    riemann::SRMHDState L, R;
    L.rho = rho_L;
    L.P_t = Pt_L;
    L.v = v_L;
    L.c_f = c_L;
    L.H = H_L;
    L.gamma_lor = gamma_L;

    R.rho = rho_R;
    R.P_t = Pt_R;
    R.v = v_R;
    R.c_f = c_R;
    R.H = H_R;
    R.gamma_lor = gamma_R;

    if (m_riemann_solver == RiemannSolverType::HLLC) {
        riemann::hllc_solver(L, R, Pt_star, v_star);
    } else {
        riemann::hll_solver(L, R, Pt_star, v_star);
    }
}

void FluidForce::moc_alfven_relativistic(
    real vy_L, real vz_L, real By_L, real Bz_L,
    real rho_L, real H_L, real gamma_L,
    real vy_R, real vz_R, real By_R, real Bz_R,
    real rho_R, real H_R, real gamma_R,
    real B_parallel,
    real& vy_star, real& vz_star, real& By_star, real& Bz_star) const
{
    // Relativistic MOC: impedance Z = sqrt(rho * H * gamma^2)
    const real Z_L = std::sqrt(std::max(rho_L * H_L * gamma_L * gamma_L, 1e-30));
    const real Z_R = std::sqrt(std::max(rho_R * H_R * gamma_R * gamma_R, 1e-30));
    const real sign_B = (B_parallel >= 0) ? 1.0 : -1.0;

    const real inv_Z_L = 1.0 / Z_L;
    const real inv_Z_R = 1.0 / Z_R;
    const real inv_Z_sum = 1.0 / (inv_Z_L + inv_Z_R);

    By_star = (By_L * inv_Z_L + By_R * inv_Z_R + sign_B * (vy_R - vy_L)) * inv_Z_sum;
    Bz_star = (Bz_L * inv_Z_L + Bz_R * inv_Z_R + sign_B * (vz_R - vz_L)) * inv_Z_sum;

    const real Z_sum = Z_L + Z_R;
    vy_star = (vy_L * Z_L + vy_R * Z_R + sign_B * (By_R - By_L)) / Z_sum;
    vz_star = (vz_L * Z_L + vz_R * Z_R + sign_B * (Bz_R - Bz_L)) / Z_sum;
}

void FluidForce::calculation(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();

    // Reset derivatives (stored in dS, de like SR-GSPH)
#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        particles[i].dS = vec_t(0.0);
        particles[i].de = 0.0;
        particles[i].dS_t = 0.0;
        particles[i].dB = vec3_t(0.0, 0.0, 0.0);
        // MHD transverse accelerations (must reset before accumulation!)
        particles[i].acc_y_mhd = 0.0;
        particles[i].acc_z_mhd = 0.0;
    }

    const real sqrt_two = std::sqrt(2.0);

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < num; ++i) {
        auto & p_i = particles[i];

        if (p_i.is_ghost) continue;

        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

#ifdef EXHAUSTIVE_SEARCH
        const int n_neighbor = exhaustive_search(p_i, p_i.sml * 6.0, particles, num,
                                                  neighbor_list, m_neighbor_number * neighbor_list_size,
                                                  periodic, false);
#else
        const int n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, false);
#endif

        // Particle i state
        const real rho_i = p_i.dens;
        const real P_i = p_i.pres;
        const real H_i = p_i.enthalpy;
        const real gamma_i = p_i.gamma_lor;
        const real h_i = p_i.sml;
        const vec3_t& B_i = p_i.B;
        const real vx_i = p_i.vel[0];
#if DIM >= 2
        const real vy_i = p_i.vel[1];
#else
        const real vy_i = p_i.vy_mhd;
#endif
#if DIM == 3
        const real vz_i = p_i.vel[2];
#else
        const real vz_i = p_i.vz_mhd;
#endif

        // Volume V = ν / N (SR-GSPH formulation)
        const real V_i = p_i.nu / p_i.N;
        const real omega_i = p_i.gradh;
        const real B_sq_i = inner_product(B_i, B_i);
        const real c_fast_i = p_i.sound;

        for (int n = 0; n < n_neighbor; ++n) {
            const int j = neighbor_list[n];
            if (j == i) continue;

            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(p_i.pos, p_j.pos);
            const real r = std::abs(r_ij);

            // SR-GSPH uses sqrt(2)*h cutoff
            if (r > 3.0 * sqrt_two * std::max(h_i, p_j.sml)) continue;

            const real h_j = p_j.sml;

            // Unit vector from j to i
#if DIM == 1
            const real e_x = (r_ij[0] >= 0) ? 1.0 : -1.0;
            const real e_y = 0.0;
            const real e_z = 0.0;
#elif DIM == 2
            const real r_inv = 1.0 / r;
            const real e_x = r_ij[0] * r_inv;
            const real e_y = r_ij[1] * r_inv;
            const real e_z = 0.0;
#else
            const real r_inv = 1.0 / r;
            const real e_x = r_ij[0] * r_inv;
            const real e_y = r_ij[1] * r_inv;
            const real e_z = r_ij[2] * r_inv;
#endif

            // Particle j state
            const real rho_j = p_j.dens;
            const real P_j = p_j.pres;
            const real H_j = p_j.enthalpy;
            const real gamma_j = p_j.gamma_lor;
            const vec3_t& B_j = p_j.B;
            const real vx_j = p_j.vel[0];
#if DIM >= 2
            const real vy_j = p_j.vel[1];
#else
            const real vy_j = p_j.vy_mhd;
#endif
#if DIM == 3
            const real vz_j = p_j.vel[2];
#else
            const real vz_j = p_j.vz_mhd;
#endif

            const real V_j = p_j.nu / p_j.N;
            const real omega_j = p_j.gradh;
            const real B_sq_j = inner_product(B_j, B_j);
            const real c_fast_j = p_j.sound;

            // Decompose velocities into normal and transverse
            const real v_n_i = vx_i * e_x + vy_i * e_y + vz_i * e_z;
            const real v_n_j = vx_j * e_x + vy_j * e_y + vz_j * e_z;

            const real vt_y_i = vy_i - e_y * v_n_i;
            const real vt_z_i = vz_i - e_z * v_n_i;
            const real vt_y_j = vy_j - e_y * v_n_j;
            const real vt_z_j = vz_j - e_z * v_n_j;

            // Decompose magnetic fields
            const real B_n_i = B_i[0] * e_x + B_i[1] * e_y + B_i[2] * e_z;
            const real B_n_j = B_j[0] * e_x + B_j[1] * e_y + B_j[2] * e_z;
            const real B_n_avg = 0.5 * (B_n_i + B_n_j);

            const real Bt_y_i = B_i[1] - e_y * B_n_i;
            const real Bt_z_i = B_i[2] - e_z * B_n_i;
            const real Bt_y_j = B_j[1] - e_y * B_n_j;
            const real Bt_z_j = B_j[2] - e_z * B_n_j;

            const real B_t_sq_i = Bt_y_i * Bt_y_i + Bt_z_i * Bt_z_i;
            const real B_t_sq_j = Bt_y_j * Bt_y_j + Bt_z_j * Bt_z_j;

            // Total pressure: gas + transverse magnetic
            const real Pt_i = P_i + 0.5 * B_t_sq_i;
            const real Pt_j = P_j + 0.5 * B_t_sq_j;

            // ================================================================
            // Riemann solver for compressive waves
            // ================================================================
            real Pt_star, v_n_star;
            real vt_y_star = 0.0, vt_z_star = 0.0;
            real Bt_y_star = 0.0, Bt_z_star = 0.0;

            if (p_j.is_ghost) {
                // Wall boundary: P* = P_i, v* = 0
                Pt_star = P_i;  // Gas pressure only for wall
                v_n_star = 0.0;
                vt_y_star = vt_y_i;
                vt_z_star = vt_z_i;
                Bt_y_star = Bt_y_i;
                Bt_z_star = Bt_z_i;
            } else if (!m_use_mhd) {
                // ============================================================
                // HYDRO-ONLY: Use SR-GSPH exact Riemann solver
                // This matches what the working SR-GSPH uses
                // ============================================================
                // Sound speed for hydro (not fast magnetosonic)
                const real cs_i = std::sqrt((m_gamma - 1.0) * (H_i - 1.0) / H_i) * m_c_speed;
                const real cs_j = std::sqrt((m_gamma - 1.0) * (H_j - 1.0) / H_j) * m_c_speed;

                // Create Riemann states (matching SR-GSPH format)
                srgsph::riemann::RiemannState left{v_n_i, rho_i, P_i, cs_i};
                srgsph::riemann::RiemannState right{v_n_j, rho_j, P_j, cs_j};

                // Use the same tangent velocity as particle i/j (no tangent velocity in 1D hydro)
                real v_t_i = 0.0, v_t_j = 0.0;
#if DIM == 1
                v_t_i = p_i.vel_t;
                v_t_j = p_j.vel_t;
#endif
                real v_t_star_temp = 0.0;

                // Call SR-GSPH exact Riemann solver
                const bool converged = srgsph::riemann::exact_riemann_solver(
                    left, right, v_t_i, v_t_j, m_gamma, m_c_speed,
                    Pt_star, v_n_star, v_t_star_temp, 100, 1e-10);

                // If exact solver fails, use HLLC acoustic fallback (like SR-GSPH)
                if (!converged) {
                    // HLLC acoustic fallback
                    const real Z_i = rho_i * H_i * gamma_i * gamma_i * cs_i;
                    const real Z_j = rho_j * H_j * gamma_j * gamma_j * cs_j;
                    const real Z_sum = Z_i + Z_j;

                    if (Z_sum > 1e-15) {
                        Pt_star = (Z_i * P_j + Z_j * P_i + Z_i * Z_j * (v_n_i - v_n_j)) / Z_sum;
                        v_n_star = (Z_i * v_n_i + Z_j * v_n_j + P_i - P_j) / Z_sum;
                    } else {
                        Pt_star = 0.5 * (P_i + P_j);
                        v_n_star = 0.5 * (v_n_i + v_n_j);
                    }

                    // Upwind tangent velocity
                    v_t_star_temp = (v_n_star > 0.0) ? v_t_i : v_t_j;
                }

                // Floor pressure
                Pt_star = std::max(Pt_star, 1e-6);

                // Clamp velocity to subluminal
                if (std::abs(v_n_star) > 0.99) {
                    v_n_star = std::copysign(0.99, v_n_star);
                }

                // Tangent velocity: upwind based on contact wave direction
                if (v_n_star > 0.0) {
                    vt_y_star = vt_y_i;
                    vt_z_star = vt_z_i;
                } else {
                    vt_y_star = vt_y_j;
                    vt_z_star = vt_z_j;
                }
            } else {
                // ============================================================
                // FULL MHD: Use SRMHD Riemann solver + MOC for Alfven
                // ============================================================
                // SR-GSPH convention: LEFT = particle i, RIGHT = particle j
                srmhd_riemann_solver(
                    rho_i, Pt_i, v_n_i, c_fast_i, H_i, gamma_i,
                    rho_j, Pt_j, v_n_j, c_fast_j, H_j, gamma_j,
                    Pt_star, v_n_star);

                // MOC for Alfven waves
                moc_alfven_relativistic(
                    vt_y_i, vt_z_i, Bt_y_i, Bt_z_i,
                    rho_i, H_i, gamma_i,
                    vt_y_j, vt_z_j, Bt_y_j, Bt_z_j,
                    rho_j, H_j, gamma_j,
                    B_n_avg,
                    vt_y_star, vt_z_star, Bt_y_star, Bt_z_star);
            }

            // Construct full interface velocity vector
            vec_t v_star_vec(0.0);
#if DIM == 1
            v_star_vec[0] = e_x * v_n_star;
#elif DIM == 2
            v_star_vec[0] = e_x * v_n_star;
            v_star_vec[1] = e_y * v_n_star + vt_y_star;
#else
            v_star_vec[0] = e_x * v_n_star;
            v_star_vec[1] = e_y * v_n_star + vt_y_star;
            v_star_vec[2] = e_z * v_n_star + vt_z_star;
#endif

            // ================================================================
            // Kernel gradients: sqrt(2)*h (SR-GSPH style)
            // ================================================================
            const vec_t grad_W_i = kernel->dw(r_ij, r, sqrt_two * h_i);
            const vec_t grad_W_j = kernel->dw(r_ij, r, sqrt_two * h_j);

            // ================================================================
            // Momentum equation: SR-GSPH volume-based formulation
            // F = -P* * (V_i² * Ω_i * ∇W_i + V_j² * Ω_j * ∇W_j)
            // ================================================================
            // For MHD: effective pressure P* - P_B_parallel
            // For hydro-only: just use the gas pressure P*
            real P_eff;
            if (m_use_mhd) {
                const real P_B_par = 0.25 * (B_n_i * B_n_i + B_n_j * B_n_j);
                P_eff = Pt_star - P_B_par;
            } else {
                P_eff = Pt_star;  // Pure gas pressure for hydro-only
            }

            const vec_t force = grad_W_i * (-P_eff * V_i * V_i * omega_i)
                              + grad_W_j * (-P_eff * V_j * V_j * omega_j);

            // Energy: power = v* · F (SR-GSPH style)
            const real power = inner_product(v_star_vec, force);

            // Accumulate into dS and de (will be normalized later)
            p_i.dS[0] += force[0];
#if DIM >= 2
            p_i.dS[1] += force[1];
#endif
#if DIM == 3
            p_i.dS[2] += force[2];
#endif
            p_i.de += power;

            // MHD-specific terms for 1D (only when MHD is enabled)
#if DIM == 1
            if (m_use_mhd) {
                // Use GSPMHD-style formulas with baryon number/density
                // GSPMHD uses: m_j/ρ² for magnetic tension, m_j/ρ for induction
                // In SR-GSPH: m → ν (baryon number), ρ → N (baryon density)
                const real nu_j = p_j.nu;
                const real N_i = p_i.N;
                const real N_j = p_j.N;
                const real N_i_inv = 1.0 / std::max(N_i, 1e-30);
                const real N_j_inv = 1.0 / std::max(N_j, 1e-30);

                // Magnetic tension force in transverse direction
                // GSPMHD formula: acc_y += B_n * ΔB_t * (m_j/ρ² * Ω * ∇W)
                // SR-GSPH: acc_y += B_n * ΔB_t * (ν_j/N² * Ω * ∇W)
                const real Bt_y_diff = Bt_y_star - Bt_y_i;
                const real Bt_z_diff = Bt_z_star - Bt_z_i;
                const real t_coeff_i = nu_j * N_i_inv * N_i_inv * omega_i;  // ν_j / N_i² * Ω_i
                const real t_coeff_j = nu_j * N_j_inv * N_j_inv * omega_j;  // ν_j / N_j² * Ω_j
                const real dw_sum = grad_W_i[0] * t_coeff_i + grad_W_j[0] * t_coeff_j;
                p_i.acc_y_mhd += B_n_avg * Bt_y_diff * dw_sum;
                p_i.acc_z_mhd += B_n_avg * Bt_z_diff * dw_sum;

                // Induction equation (1D) - Iwasaki & Inutsuka (2011) Eq. (53)
                // TEMPORARILY DISABLED: induction causes exponential B field growth
                // TODO: Debug why induction is unstable (possible feedback loop)
                // Keep B field frozen for now to test magnetic tension alone
#if 0  // DISABLED FOR DEBUGGING
                {
                const real dv_n_i = v_n_star - v_n_i;
                const real dv_n_j = v_n_star - v_n_j;
                const real ind_coeff_i = nu_j * N_i_inv * omega_i;  // ν_j / N_i * Ω_i
                const real ind_coeff_j = nu_j * N_j_inv * omega_j;  // ν_j / N_j * Ω_j
                p_i.dB[1] += Bt_y_i * dv_n_i * ind_coeff_i * grad_W_i[0]
                           + Bt_y_j * dv_n_j * ind_coeff_j * grad_W_j[0];
                p_i.dB[2] += Bt_z_i * dv_n_i * ind_coeff_i * grad_W_i[0]
                           + Bt_z_j * dv_n_j * ind_coeff_j * grad_W_j[0];
                }
#endif

                // Suppress unused variable warnings for star values (used in higher-D)
                (void)vt_y_star; (void)vt_z_star; (void)Bt_y_star; (void)Bt_z_star;
            }
#endif
        }
    }

    // Normalize by baryon number ν (SR-GSPH style)
#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        if (particles[i].is_ghost) continue;
        normalize_srmhd_derivatives(particles[i]);
    }
}

}
}
