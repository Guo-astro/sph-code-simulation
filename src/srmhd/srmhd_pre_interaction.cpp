/**
 * SRMHD Pre-Interaction Module
 *
 * Special Relativistic MHD pre-interaction combining:
 * - Kitajima et al. (2025) volume-based density and smoothing length
 * - Iwasaki & Inutsuka (2011) MHD div(B) calculation
 *
 * Performs per-timestep setup:
 * 1. Compute smoothing length h (volume-based iteration)
 * 2. Compute particle volume V_p and lab-frame density N
 * 3. Recover primitive variables from conserved (with MHD wave speeds)
 * 4. Compute divergence of B for Powell correction
 * 5. Compute gradients for MUSCL reconstruction (if enabled)
 */

#include "srmhd/srmhd_pre_interaction.hpp"
#include "srmhd/srmhd_primitive_recovery.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "periodic.hpp"
#include "openmp.hpp"
#include "kernel/kernel_function.hpp"
#include "bhtree.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>


namespace sph
{
namespace srmhd
{

void PreInteraction::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::PreInteraction::initialize(param);

    m_gamma = param->physics.gamma;
    m_c_speed = param->srgsph.c_speed;
    m_is_2nd_order = param->srgsph.is_2nd_order;
    m_iteration = param->iterative_sml;
    m_eta = param->srgsph.eta;
    m_c_smooth = 2.0;
    m_use_powell = param->mhd.use_powell_correction;
    m_first_calculation = true;

    std::cout << "[SRMHD PreInteraction] Initialized" << std::endl;
    std::cout << "  eta: " << m_eta << ", c: " << m_c_speed << std::endl;
    std::cout << "  Powell div-B correction: " << (m_use_powell ? "yes" : "no") << std::endl;
}

real PreInteraction::compute_volume(
    const SPHParticle & p_i,
    const std::vector<SPHParticle> & particles,
    const std::vector<int> & neighbor_list,
    const int n_neighbor,
    const Periodic * periodic,
    const KernelFunction * kernel,
    const real h,
    real & gradh_out
)
{
    real sum_W = 0.0;
    real sum_dW_dh = 0.0;

    for (int n = 0; n < n_neighbor; ++n) {
        const int j = neighbor_list[n];
        const auto & p_j = particles[j];

        const vec_t r_ij = periodic->calc_r_ij(p_i.pos, p_j.pos);
        const real r = std::abs(r_ij);

        sum_W += kernel->w(r, h);
        sum_dW_dh += kernel->dhw(r, h);
    }

    if (sum_W < 1e-15) {
        gradh_out = 1.0;
        return 1.0;
    }

    // Grad-h correction factor: Omega = 1 / (1 + h * sum dW/dh / (D * sum W))
    const real dh_term = h * sum_dW_dh / (DIM * sum_W);

    // Safety check: if dh_term approaches -1, gradh becomes infinite
    // Limit dh_term to prevent this singularity
    const real dh_term_safe = std::max(dh_term, -0.9);  // Ensure 1 + dh_term >= 0.1
    gradh_out = 1.0 / (1.0 + dh_term_safe);

    // Also limit gradh to reasonable range [0.1, 10]
    gradh_out = std::max(0.1, std::min(10.0, gradh_out));

    return 1.0 / sum_W;
}

real PreInteraction::compute_smoothing_length(
    const SPHParticle & p_i,
    const std::vector<SPHParticle> & particles,
    const std::vector<int> & neighbor_list,
    const int n_neighbor,
    const Periodic * periodic,
    const KernelFunction * kernel,
    const real search_radius
)
{
    real h = p_i.sml;

    constexpr int max_iter = 50;
    for (int iter = 0; iter < max_iter; ++iter) {
        real sum_W_star = 0.0;

        for (int n = 0; n < n_neighbor; ++n) {
            const int j = neighbor_list[n];
            const auto & p_j = particles[j];

            const vec_t r_ij = periodic->calc_r_ij(p_i.pos, p_j.pos);
            const real r = std::abs(r_ij);

            sum_W_star += kernel->w(r, m_c_smooth * h);
        }

        if (sum_W_star < 1e-15) {
            break;
        }

        const real Vp_star = 1.0 / sum_W_star;
        const real h_new = m_eta * std::pow(Vp_star, 1.0 / DIM);

        if (std::abs(h_new - h) / h < 1e-12) {
            h = h_new;
            break;
        }

        h = h_new;
    }

    return h;
}

real PreInteraction::compute_div_B(
    const SPHParticle & p_i,
    const std::vector<SPHParticle> & particles,
    const std::vector<int> & neighbor_list,
    const int n_neighbor,
    const Periodic * periodic,
    const KernelFunction * kernel
)
{
    const vec_t& r_i = p_i.pos;
    const vec3_t& B_i = p_i.B;
    const real h_i = p_i.sml;
    const real V_i = p_i.nu / p_i.N;

    real div_B = 0.0;

    for (int n = 0; n < n_neighbor; ++n) {
        const int j = neighbor_list[n];
        const auto & p_j = particles[j];

        const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
        const real r = std::abs(r_ij);

        if (r >= std::max(h_i, p_j.sml) || r == 0.0) {
            continue;
        }

        const real r_inv = 1.0 / r;

        // Unit vector
#if DIM == 1
        const real n_x = (r_ij[0] >= 0) ? 1.0 : -1.0;
        const real n_y = 0.0;
        const real n_z = 0.0;
#elif DIM == 2
        const real n_x = r_ij[0] * r_inv;
        const real n_y = r_ij[1] * r_inv;
        const real n_z = 0.0;
#else
        const real n_x = r_ij[0] * r_inv;
        const real n_y = r_ij[1] * r_inv;
        const real n_z = r_ij[2] * r_inv;
#endif

        // B_parallel* = (B_i + B_j) . n_ij / 2
        const vec3_t& B_j = p_j.B;
        const real B_par_i = B_i[0] * n_x + B_i[1] * n_y + B_i[2] * n_z;
        const real B_par_j = B_j[0] * n_x + B_j[1] * n_y + B_j[2] * n_z;
        const real B_par_star = 0.5 * (B_par_i + B_par_j);

        // Kernel gradient
        const vec_t dw_i = kernel->dw(r_ij, r, h_i);
        const vec_t dw_j = kernel->dw(r_ij, r, p_j.sml);

        // Volume-based formulation
        const real V_j = p_j.nu / p_j.N;
        const real omega_i = p_i.gradh;
        const real omega_j = p_j.gradh;

        // F_ij for div B calculation (antisymmetric)
#if DIM == 1
        const real dw_dot_n_i = dw_i[0] * n_x;
        const real dw_dot_n_j = dw_j[0] * n_x;
#elif DIM == 2
        const real dw_dot_n_i = dw_i[0] * n_x + dw_i[1] * n_y;
        const real dw_dot_n_j = dw_j[0] * n_x + dw_j[1] * n_y;
#else
        const real dw_dot_n_i = dw_i[0] * n_x + dw_i[1] * n_y + dw_i[2] * n_z;
        const real dw_dot_n_j = dw_j[0] * n_x + dw_j[1] * n_y + dw_j[2] * n_z;
#endif
        const real F_ij = dw_dot_n_i * V_i * V_i * omega_i
                        + dw_dot_n_j * V_j * V_j * omega_j;

        // Divergence of B
        div_B += B_par_star * F_ij;
    }

    return div_B;
}

void PreInteraction::initial_smoothing(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    const int num = sim->get_particle_num();

#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        auto & p_i = particles[i];

        // If sml is already set, don't overwrite
        if (p_i.sml > 1.0e-10) continue;

        // Initial guess
        constexpr real A = DIM == 1 ? 2.0 :
                           DIM == 2 ? M_PI :
                                      4.0 * M_PI / 3.0;

        p_i.sml = std::pow(m_neighbor_number * p_i.mass / (p_i.dens * A), 1.0 / DIM) * m_kernel_ratio;
    }
}

void PreInteraction::calculation(std::shared_ptr<Simulation> sim)
{
    if (m_first) {
        initial_smoothing(sim);
        m_first = false;
    }

    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();
    auto * tree = sim->get_tree().get();

#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        auto & p_i = particles[i];

        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

        // Extended search radius for volume-based h iteration
        const real search_r = p_i.sml * 6.0;

        const int n_neighbor_tmp = tree->neighbor_search(p_i, neighbor_list, particles, false);

        // Update smoothing length (skip on first timestep)
        if (m_iteration && !m_first_calculation) {
            p_i.sml = compute_smoothing_length(p_i, particles, neighbor_list,
                                               n_neighbor_tmp, periodic, kernel, p_i.sml);
        }

        // Ghost particles: only compute h, skip density/primitive updates
        if (p_i.is_ghost) {
            p_i.gradh = 1.0;
            // Set fast wave speed for ghost particles using RELATIVISTIC formula
            // Non-relativistic formula can exceed c, causing timestep collapse!
            const real rho = std::max(p_i.dens, 1e-30);
            const real P = p_i.pres;
            const real u_int = P / ((m_gamma - 1.0) * rho);
            const real H = 1.0 + u_int / (m_c_speed * m_c_speed)
                         + P / (rho * m_c_speed * m_c_speed);

            // Relativistic sound speed squared: c_s^2 = (gamma-1)(H-1)/H * c^2
            const real c_s2 = (m_gamma - 1.0) * (H - 1.0) / H * m_c_speed * m_c_speed;

            // Relativistic Alfven speed: c_A^2 = B^2 / (rho*H*gamma^2 + B^2)
            // For ghost at rest, gamma_lor = 1
            const real B_sq = inner_product(p_i.B, p_i.B);
            const real rho_H = rho * H;
            const real c_A2 = B_sq / (rho_H + B_sq) * m_c_speed * m_c_speed;

            // Relativistic fast magnetosonic: c_f^2 = c_s^2 + c_A^2 - c_s^2*c_A^2/c^2
            // This formula guarantees c_f <= c
            const real c_f2 = c_s2 + c_A2 - c_s2 * c_A2 / (m_c_speed * m_c_speed);
            p_i.sound = std::sqrt(std::max(c_f2, 0.0));
            continue;
        }

        // Compute particle volume and grad-h correction
        real gradh_i;
        const real Vp = compute_volume(p_i, particles, neighbor_list,
                                       n_neighbor_tmp, periodic, kernel, p_i.sml, gradh_i);
        p_i.gradh = gradh_i;

        // Number density: N = nu / V_p
        // Add floor to prevent N from becoming too small, which causes
        // V = nu/N to become large and destabilizes the force calculation
        const real N_floor = p_i.nu / (100.0 * p_i.sml);  // Corresponds to particle spacing = 100*h
        const real N_new = std::max(p_i.nu / Vp, N_floor);

        // First timestep: recompute conserved variables from primitives
        if (m_first_calculation) {
            // Compute total velocity squared including MHD components
            real v2 = inner_product(p_i.vel, p_i.vel);
#if DIM == 1
            const real vy = p_i.vy_mhd;
            const real vz = p_i.vz_mhd;
            v2 += vy * vy + vz * vz;
#endif
            const real gamma_lor = 1.0 / std::sqrt(1.0 - v2 / (m_c_speed * m_c_speed));
            p_i.gamma_lor = gamma_lor;

            // Rest-frame density from lab-frame
            const real rest_density = std::max(N_new / gamma_lor, 1e-20);
            const real u_internal = p_i.pres / ((m_gamma - 1.0) * rest_density);
            const real H = 1.0 + u_internal / (m_c_speed * m_c_speed)
                         + p_i.pres / (rest_density * m_c_speed * m_c_speed);
            p_i.enthalpy = H;

            // Conserved variables
            p_i.N = N_new;
            p_i.S = p_i.vel * (gamma_lor * H);
#if DIM == 1
            p_i.S_t = std::sqrt(vy * vy + vz * vz) * (gamma_lor * H);
#else
            p_i.S_t = 0.0;
#endif
            const real X = m_gamma / (m_gamma - 1.0);
            p_i.e = (H * (X * gamma_lor * gamma_lor - 1.0) + 1.0) / (X * gamma_lor);

            // Update primitive storage
            p_i.dens = rest_density;
            p_i.ene = u_internal;

            // MHD wave speeds
            const real B_sq = inner_product(p_i.B, p_i.B);
            const real c_s = std::sqrt((m_gamma - 1.0) * (H - 1.0) / H) * m_c_speed;
            const real rho_H_gamma2 = rest_density * H * gamma_lor * gamma_lor;
            const real c_A = std::sqrt(B_sq / (rho_H_gamma2 + B_sq));
            const real c_f2 = c_s * c_s + c_A * c_A - c_s * c_s * c_A * c_A;
            p_i.sound = std::sqrt(std::max(c_f2, 0.0));
        } else {
            // Normal timestep: recover primitives from conserved with MHD
            p_i.N = N_new;

            const auto prim = PrimitiveRecovery::conserved_to_primitive_mhd(
                p_i.S, p_i.S_t, p_i.e, p_i.N, p_i.B, m_gamma, m_c_speed);

            // Store primitive variables
            p_i.vel = prim.vel;
            p_i.vel_t = prim.vel_t;
            p_i.dens = prim.density;
            p_i.pres = prim.pressure;
            p_i.gamma_lor = prim.gamma_lor;
            p_i.enthalpy = prim.enthalpy;

            // Store MHD wave speeds
            p_i.sound = prim.fast_speed;  // Fast magnetosonic for timestep
        }

        // Compute divergence of B for Powell correction
        if (m_use_powell) {
            p_i.div_B = compute_div_B(p_i, particles, neighbor_list,
                                       n_neighbor_tmp, periodic, kernel);
        }

        p_i.neighbor = n_neighbor_tmp;
    }

    // After first calculation, switch to normal mode
    m_first_calculation = false;
}

} // namespace srmhd
} // namespace sph
