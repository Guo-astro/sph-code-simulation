#include "defines.hpp"
#include "particle.hpp"
#include "periodic.hpp"
#include "simulation.hpp"
#include "bhtree.hpp"
#include "kernel/kernel_function.hpp"
#include "gspmhd/mhd_pre_interaction.hpp"

#include <cmath>
#include <iostream>

#ifdef EXHAUSTIVE_SEARCH
#include "exhaustive_search.hpp"
#endif

namespace sph
{
namespace gspmhd
{

void PreInteraction::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::PreInteraction::initialize(param);
    m_gamma = param->physics.gamma;
    m_use_powell = param->mhd.use_powell_correction;

    std::cout << "[GSPMHD PreInteraction] Initialized" << std::endl;
}

void PreInteraction::calculation(std::shared_ptr<Simulation> sim)
{
    // First, call the base class to compute density, smoothing length, and grad-h
    sph::PreInteraction::calculation(sim);

    // Now compute MHD-specific quantities: div B and sound speed with magnetic pressure
    auto& particles = sim->get_particles();
    auto* periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto* kernel = sim->get_kernel().get();
    auto* tree = sim->get_tree().get();

#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        auto& p_i = particles[i];

        // Ghost particles: skip MHD-specific computation, keep initial values
        // Their properties are set by initialization and preserved
        if (p_i.is_ghost) {
            // Set sound speed including magnetic contribution for ghost particles
            const real B_sq = inner_product(p_i.B, p_i.B);
            const real c_sound = std::sqrt(m_gamma * p_i.pres / p_i.dens);
            const real c_alfven = std::sqrt(B_sq / p_i.dens);
            p_i.sound = std::sqrt(c_sound * c_sound + c_alfven * c_alfven);
            continue;
        }

        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

        // Neighbor search
#ifdef EXHAUSTIVE_SEARCH
        int const n_neighbor = exhaustive_search(p_i, p_i.sml, particles, num, neighbor_list, m_neighbor_number * neighbor_list_size, periodic, true);
#else
        int const n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, true);
#endif

        const vec_t& r_i = p_i.pos;
        const vec3_t& B_i = p_i.B;
        const real h_i = p_i.sml;
        const real rho_i = p_i.dens;

        real div_B = 0.0;

        for (int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto& p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
            const real r = std::abs(r_ij);

            // Use same cutoff as GSPH - smoothing length, not extended support
            if (r >= std::max(h_i, p_j.sml) || r == 0.0) {
                continue;
            }

            const real r_inv = 1.0 / r;

            // Unit vector components (for 1D, only x matters)
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

            // Compute B_parallel* = (B_i + B_j) . n_ij / 2
            const vec3_t& B_j = p_j.B;
            const real B_par_i = B_i[0] * n_x + B_i[1] * n_y + B_i[2] * n_z;
            const real B_par_j = B_j[0] * n_x + B_j[1] * n_y + B_j[2] * n_z;
            const real B_par_star = 0.5 * (B_par_i + B_par_j);

            // Kernel gradient
            const vec_t dw_i = kernel->dw(r_ij, r, h_i);
            const vec_t dw_j = kernel->dw(r_ij, r, p_j.sml);

            const real rho2_inv_i = 1.0 / sqr(rho_i);
            const real rho2_inv_j = 1.0 / sqr(p_j.dens);

            // F_ij for div B calculation
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
            const real F_ij = dw_dot_n_i * rho2_inv_i * p_j.mass
                            + dw_dot_n_j * rho2_inv_j * p_j.mass;

            // Divergence of B (Eq. 31)
            div_B += B_par_star * F_ij;
        }

        // Store div B / rho for monitoring and Powell correction
        p_i.div_B = div_B * rho_i;

        // Update sound speed to include magnetic contribution (fast wave speed)
        const real B_sq = inner_product(B_i, B_i);
        const real c_sound = std::sqrt(m_gamma * p_i.pres / rho_i);
        const real c_alfven = std::sqrt(B_sq / rho_i);
        p_i.sound = std::sqrt(c_sound * c_sound + c_alfven * c_alfven);  // Fast wave speed
    }
}

}
}
