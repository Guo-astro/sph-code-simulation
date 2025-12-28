/**
 * @file mhd_fluid_force.cpp
 * @brief Godunov SPMHD fluid force calculation
 *
 * Implementation of Iwasaki & Inutsuka (2011) GSPMHD method.
 */

#include "defines.hpp"
#include "particle.hpp"
#include "periodic.hpp"
#include "simulation.hpp"
#include "bhtree.hpp"
#include "kernel/kernel_function.hpp"
#include "gspmhd/mhd_fluid_force.hpp"

#include <cmath>
#include <algorithm>
#include <iostream>

#ifdef EXHAUSTIVE_SEARCH
#include "exhaustive_search.hpp"
#endif

namespace sph
{
namespace gspmhd
{

void FluidForce::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::FluidForce::initialize(param);
    m_is_2nd_order = param->mhd.is_2nd_order;
    m_gamma = param->physics.gamma;
    m_use_powell = param->mhd.use_powell_correction;

    std::cout << "[GSPMHD] Godunov SPMHD initialized (Iwasaki & Inutsuka 2011)" << std::endl;
    std::cout << "  2nd order: " << (m_is_2nd_order ? "yes" : "no") << std::endl;
    std::cout << "  Powell div-B correction: " << (m_use_powell ? "yes" : "no") << std::endl;
    std::cout << "  gamma: " << m_gamma << std::endl;
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

/**
 * @brief HLL-type MHD Riemann solver (same formula as GSPH)
 */
void FluidForce::mhd_riemann_solver(
    real rho_L, real Pt_L, real v_L, real c_L,
    real rho_R, real Pt_R, real v_R, real c_R,
    real& Pt_star, real& v_star) const
{
    const real roe_l = std::sqrt(std::max(rho_L, 1e-30));
    const real roe_r = std::sqrt(std::max(rho_R, 1e-30));
    const real roe_inv = 1.0 / (roe_l + roe_r);

    const real u_t = (roe_l * v_L + roe_r * v_R) * roe_inv;
    const real c_t = (roe_l * c_L + roe_r * c_R) * roe_inv;
    const real s_l = std::min(v_L - c_L, u_t - c_t);
    const real s_r = std::max(v_R + c_R, u_t + c_t);

    const real c1 = rho_L * (s_l - v_L);
    const real c2 = rho_R * (s_r - v_R);
    const real c3 = 1.0 / (c1 - c2 + 1e-30);
    const real c4 = Pt_L - v_L * c1;
    const real c5 = Pt_R - v_R * c2;

    v_star = (c5 - c4) * c3;
    Pt_star = (c1 * c5 - c2 * c4) * c3;
    Pt_star = std::max(Pt_star, 1e-30);
}

/**
 * @brief Method of Characteristics for Alfven waves (Eq. 55-56)
 */
void FluidForce::moc_alfven_3d(
    real vy_L, real vz_L, real By_L, real Bz_L, real rho_L,
    real vy_R, real vz_R, real By_R, real Bz_R, real rho_R,
    real B_parallel,
    real& vy_star, real& vz_star, real& By_star, real& Bz_star) const
{
    const real sqrt_rho_L = std::sqrt(std::max(rho_L, 1e-30));
    const real sqrt_rho_R = std::sqrt(std::max(rho_R, 1e-30));
    const real sign_B = (B_parallel >= 0) ? 1.0 : -1.0;

    const real inv_sqrt_L = 1.0 / sqrt_rho_L;
    const real inv_sqrt_R = 1.0 / sqrt_rho_R;
    const real inv_sqrt_sum = 1.0 / (inv_sqrt_L + inv_sqrt_R);

    By_star = (By_L * inv_sqrt_L + By_R * inv_sqrt_R + (vy_R - vy_L) * sign_B) * inv_sqrt_sum;
    Bz_star = (Bz_L * inv_sqrt_L + Bz_R * inv_sqrt_R + (vz_R - vz_L) * sign_B) * inv_sqrt_sum;

    const real sqrt_sum = sqrt_rho_L + sqrt_rho_R;
    vy_star = (vy_L * sqrt_rho_L + vy_R * sqrt_rho_R + (By_R - By_L) * sign_B) / sqrt_sum;
    vz_star = (vz_L * sqrt_rho_L + vz_R * sqrt_rho_R + (Bz_R - Bz_L) * sign_B) / sqrt_sum;
}

void FluidForce::calculation(std::shared_ptr<Simulation> sim)
{
    auto& particles = sim->get_particles();
    auto* periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto* kernel = sim->get_kernel().get();
    auto* tree = sim->get_tree().get();

#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        auto& p_i = particles[i];
        if (p_i.is_ghost) continue;

        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

#ifdef EXHAUSTIVE_SEARCH
        int const n_neighbor = exhaustive_search(p_i, p_i.sml, particles, num,
            neighbor_list, m_neighbor_number * neighbor_list_size, periodic, true);
#else
        int const n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, true);
#endif

        const vec_t& r_i = p_i.pos;
        const real h_i = p_i.sml;
        const real rho_i = p_i.dens;
        const real P_i = p_i.pres;
        const real omega_i = p_i.gradh;

        const real vx_i = p_i.vel[0];
        const real vy_i = p_i.vy_mhd;
        const real vz_i = p_i.vz_mhd;
        const vec3_t& B_i = p_i.B;

        const real rho_i_inv = 1.0 / std::max(rho_i, 1e-30);
        const real rho2_inv_i = rho_i_inv * rho_i_inv;
        const real B_sq_i = inner_product(B_i, B_i);
        const real c_fast_i = std::sqrt((m_gamma * P_i + B_sq_i) * rho_i_inv);

        real acc_x = 0.0, acc_y = 0.0, acc_z = 0.0;
        real dene = 0.0;
        real dBx = 0.0, dBy = 0.0, dBz = 0.0;

        for (int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            if (i == j) continue;
            auto& p_j = particles[j];

            const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
            const real r = std::abs(r_ij);

            // Use h cutoff like GSPH
            if (r >= std::max(h_i, p_j.sml) || r == 0.0) continue;

            const real r_inv = 1.0 / r;
            const real vx_j = p_j.vel[0];
            const real vy_j = p_j.vy_mhd;
            const real vz_j = p_j.vz_mhd;
            const vec3_t& B_j = p_j.B;
            const real rho_j = p_j.dens;
            const real P_j = p_j.pres;
            const real omega_j = p_j.gradh;

            const real rho_j_inv = 1.0 / std::max(rho_j, 1e-30);
            const real rho2_inv_j = rho_j_inv * rho_j_inv;
            const real B_sq_j = inner_product(B_j, B_j);
            const real c_fast_j = std::sqrt((m_gamma * P_j + B_sq_j) * rho_j_inv);

            // Unit vector
#if DIM == 1
            const real e_x = (r_ij[0] >= 0) ? 1.0 : -1.0;
            const real e_y = 0.0, e_z = 0.0;
#elif DIM == 2
            const real e_x = r_ij[0] * r_inv;
            const real e_y = r_ij[1] * r_inv;
            const real e_z = 0.0;
#else
            const real e_x = r_ij[0] * r_inv;
            const real e_y = r_ij[1] * r_inv;
            const real e_z = r_ij[2] * r_inv;
#endif

            // Normal velocities
            const real v_n_i = vx_i * e_x + vy_i * e_y + vz_i * e_z;
            const real v_n_j = vx_j * e_x + vy_j * e_y + vz_j * e_z;

            // Normal B fields
            const real B_n_i = B_i[0] * e_x + B_i[1] * e_y + B_i[2] * e_z;
            const real B_n_j = B_j[0] * e_x + B_j[1] * e_y + B_j[2] * e_z;
            const real B_n_avg = 0.5 * (B_n_i + B_n_j);

            // Transverse velocities
            const real vt_y_i = vy_i - e_y * v_n_i;
            const real vt_z_i = vz_i - e_z * v_n_i;
            const real vt_y_j = vy_j - e_y * v_n_j;
            const real vt_z_j = vz_j - e_z * v_n_j;

            // Transverse B fields
            const real Bt_y_i = B_i[1] - e_y * B_n_i;
            const real Bt_z_i = B_i[2] - e_z * B_n_i;
            const real Bt_y_j = B_j[1] - e_y * B_n_j;
            const real Bt_z_j = B_j[2] - e_z * B_n_j;

            const real B_t_sq_i = Bt_y_i * Bt_y_i + Bt_z_i * Bt_z_i;
            const real B_t_sq_j = Bt_y_j * Bt_y_j + Bt_z_j * Bt_z_j;

            // Total pressure (gas + magnetic)
            const real Pt_i = P_i + 0.5 * B_t_sq_i;
            const real Pt_j = P_j + 0.5 * B_t_sq_j;

            // Riemann solver using fast magnetosonic speed
            real Pt_star, v_n_star;
            mhd_riemann_solver(rho_j, Pt_j, v_n_j, c_fast_j,
                               rho_i, Pt_i, v_n_i, c_fast_i,
                               Pt_star, v_n_star);

            // MOC for transverse
            real vt_y_star, vt_z_star, Bt_y_star, Bt_z_star;
            moc_alfven_3d(vt_y_j, vt_z_j, Bt_y_j, Bt_z_j, rho_j,
                          vt_y_i, vt_z_i, Bt_y_i, Bt_z_i, rho_i,
                          B_n_avg, vt_y_star, vt_z_star, Bt_y_star, Bt_z_star);

            // Interface velocity vector (like GSPH: v_ij = e_ij * vstar)
            const real v_ij_x = e_x * v_n_star;
            const real v_ij_y = e_y * v_n_star;
            const real v_ij_z = e_z * v_n_star;

            // Kernel gradients
            const vec_t dw_i = kernel->dw(r_ij, r, h_i);
            const vec_t dw_j = kernel->dw(r_ij, r, p_j.sml);

            // Momentum equation (Eq. 52)
            const real P_B_par = 0.25 * (B_n_i * B_n_i + B_n_j * B_n_j);
            const real P_eff = Pt_star - P_B_par;

            // Force like GSPH: f = dw_i * coeff_i + dw_j * coeff_j
            const real coeff_i = p_j.mass * P_eff * rho2_inv_i * omega_i;
            const real coeff_j = p_j.mass * P_eff * rho2_inv_j * omega_j;

#if DIM == 1
            const real f_x = dw_i[0] * coeff_i + dw_j[0] * coeff_j;
            acc_x -= f_x;
            // Energy like GSPH: dene -= inner_product(f, v_ij - v_i)
            dene -= f_x * (v_ij_x - vx_i);
#else
            const real f_x = dw_i[0] * coeff_i + dw_j[0] * coeff_j;
            const real f_y = dw_i[1] * coeff_i + dw_j[1] * coeff_j;
#if DIM == 3
            const real f_z = dw_i[2] * coeff_i + dw_j[2] * coeff_j;
#else
            const real f_z = 0.0;
#endif
            acc_x -= f_x;
            acc_y -= f_y;
            acc_z -= f_z;
            dene -= (f_x * (v_ij_x - vx_i) + f_y * (v_ij_y - vy_i) + f_z * (v_ij_z - vz_i));
#endif

            // Magnetic tension (1D)
#if DIM == 1
            const real Bt_y_diff = Bt_y_star - Bt_y_i;
            const real Bt_z_diff = Bt_z_star - Bt_z_i;
            const real t_coeff_i = p_j.mass * rho2_inv_i * omega_i;
            const real t_coeff_j = p_j.mass * rho2_inv_j * omega_j;
            const real dw_sum = dw_i[0] * t_coeff_i + dw_j[0] * t_coeff_j;
            acc_y += B_n_avg * Bt_y_diff * dw_sum;
            acc_z += B_n_avg * Bt_z_diff * dw_sum;
#endif

            // Induction equation (1D)
#if DIM == 1
            const real dvy = vt_y_star - vt_y_i;
            const real dvz = vt_z_star - vt_z_i;
            const real ind_coeff_i = p_j.mass * rho_i_inv * omega_i;
            const real ind_coeff_j = p_j.mass * rho_j_inv * omega_j;
            const real f_ind_x = dw_i[0] * ind_coeff_i + dw_j[0] * ind_coeff_j;
            const real Bn_gradW = B_n_avg * f_ind_x * e_x;
            dBy += dvy * Bn_gradW;
            dBz += dvz * Bn_gradW;
#endif
        }

        p_i.acc[0] = acc_x;
#if DIM >= 2
        p_i.acc[1] = acc_y;
#endif
#if DIM >= 3
        p_i.acc[2] = acc_z;
#endif
        p_i.acc_y_mhd = acc_y;
        p_i.acc_z_mhd = acc_z;
        p_i.dene = dene;

#if DIM == 1
        p_i.dB = vec3_t(0.0, dBy, dBz);
#else
        p_i.dB = vec3_t(dBx, dBy, dBz);
#endif
    }
}

} // namespace gspmhd
} // namespace sph
