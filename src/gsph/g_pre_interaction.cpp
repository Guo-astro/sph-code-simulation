#include <algorithm>
#include <cmath>

#include "parameters.hpp"
#include "gsph/g_pre_interaction.hpp"
#include "simulation.hpp"
#include "periodic.hpp"
#include "openmp.hpp"
#include "kernel/kernel_function.hpp"
#include "exception.hpp"
#include "bhtree.hpp"


namespace sph
{
namespace gsph
{

void PreInteraction::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::PreInteraction::initialize(param);
    m_is_2nd_order = param->gsph.is_2nd_order;
    m_use_gradh = param->gsph.use_gradh;
    m_use_volume_based = param->gsph.use_volume_based;
    m_eta = param->gsph.eta;
    m_c_smooth = param->gsph.c_smooth;
    m_first_calculation = true;
}

/**
 * Compute particle volume V_p(x_i) = [Σ_j W(x_i - x_j, h)]^(-1)
 * and grad-h correction factor Ω_i (volume-based formulation)
 *
 * Based on Kitajima et al.: V_p = [Σ_j W]^(-1)
 * Grad-h correction: Ω = 1 / (1 + h * Σ dW/dh / (D * Σ W))
 */
real PreInteraction::compute_volume(
    const SPHParticle & p_i,
    const std::vector<SPHParticle> & particles,
    const std::vector<int> & neighbor_list,
    const int n_neighbor,
    const Periodic * periodic,
    const KernelFunction * kernel,
    const real h,
    real & gradh_out)
{
    real sum_W = 0.0;
    real sum_dW_dh = 0.0;

    for (int n = 0; n < n_neighbor; ++n) {
        const int j = neighbor_list[n];
        const auto & p_j = particles[j];

        const vec_t r_ij = periodic->calc_r_ij(p_i.pos, p_j.pos);
        const real r = std::abs(r_ij);

        if (r >= h) continue;

        sum_W += kernel->w(r, h);
        sum_dW_dh += kernel->dhw(r, h);
    }

    if (sum_W < 1e-15) {
        gradh_out = 1.0;
        return 1.0;
    }

    // Grad-h correction factor: Ω = 1 / (1 + h * Σ dW/dh / (D * Σ W))
    const real dh_term = h * sum_dW_dh / (DIM * sum_W);
    gradh_out = 1.0 / (1.0 + dh_term);

    return 1.0 / sum_W;
}

/**
 * Compute variable smoothing length using volume-based approach
 * Iteratively solves: h = η * [V_p*]^(1/d) where V_p* = [Σ_j W(r, c_smooth*h)]^(-1)
 */
real PreInteraction::compute_smoothing_length_volume(
    const SPHParticle & p_i,
    const std::vector<SPHParticle> & particles,
    const std::vector<int> & neighbor_list,
    const int n_neighbor,
    const Periodic * periodic,
    const KernelFunction * kernel)
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

        if (sum_W_star < 1e-15) break;

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

void PreInteraction::calculation(std::shared_ptr<Simulation> sim)
{
    if(m_first) {
        initial_smoothing(sim);
        m_first = false;
    }

    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();
    auto * tree = sim->get_tree().get();

    omp_real h_per_v_sig(std::numeric_limits<real>::max());

    // for MUSCL
    auto & grad_d = sim->get_vector_array("grad_density");
    auto & grad_p = sim->get_vector_array("grad_pressure");
    vec_t * grad_v[DIM] = {
        sim->get_vector_array("grad_velocity_0").data(),
#if DIM == 2
        sim->get_vector_array("grad_velocity_1").data(),
#elif DIM == 3
        sim->get_vector_array("grad_velocity_1").data(),
        sim->get_vector_array("grad_velocity_2").data(),
#endif
    };

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];

        // Ghost particles: skip density computation, keep initial values
        // Their properties are set by initialization and maintained by update_ghost_particles()
        if (p_i.is_ghost) {
            p_i.gradh = 1.0;  // Set gradh=1.0 so forces don't vanish at boundaries
            continue;
        }

        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

        // guess smoothing length
        constexpr real A = DIM == 1 ? 2.0 :
                           DIM == 2 ? M_PI :
                                      4.0 * M_PI / 3.0;
        p_i.sml = std::pow(m_neighbor_number * p_i.mass / (p_i.dens * A), 1.0 / DIM) * m_kernel_ratio;

        // neighbor search - use larger search radius for volume-based approach
        const real search_radius = m_use_volume_based ? p_i.sml * 6.0 : p_i.sml;
        const int n_neighbor_tmp = tree->neighbor_search(p_i, neighbor_list, particles, false);

        // smoothing length update
        if(m_iteration) {
            if(m_use_volume_based && !m_first_calculation) {
                // Volume-based smoothing length: h = η * V_p*^(1/d)
                p_i.sml = compute_smoothing_length_volume(p_i, particles, neighbor_list,
                                                          n_neighbor_tmp, periodic, kernel);
            } else {
                // Standard mass-based smoothing length
                p_i.sml = newton_raphson(p_i, particles, neighbor_list, n_neighbor_tmp, periodic, kernel);
            }
        }

        // Compute density and grad-h based on approach
        real dens_i = 0.0;
        real v_sig_max = p_i.sound * 2.0;
        const vec_t & pos_i = p_i.pos;
        int n_neighbor = 0;

        if(m_use_volume_based) {
            // ===== VOLUME-BASED APPROACH (Kitajima et al.) =====
            // Compute particle volume V_p and grad-h correction
            real gradh_i;
            const real Vp = compute_volume(p_i, particles, neighbor_list,
                                           n_neighbor_tmp, periodic, kernel, p_i.sml, gradh_i);

            // Density from volume: ρ = m / V_p
            dens_i = p_i.mass / Vp;
            p_i.nu = Vp;

            if(m_use_gradh) {
                p_i.gradh = gradh_i;
            } else {
                p_i.gradh = 1.0;
            }

            // Count neighbors and compute signal velocity
            for(int n = 0; n < n_neighbor_tmp; ++n) {
                int const j = neighbor_list[n];
                auto & p_j = particles[j];
                const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
                const real r = std::abs(r_ij);

                if(r >= p_i.sml) continue;

                ++n_neighbor;

                if(i != j) {
                    const real v_sig = p_i.sound + p_j.sound - 3.0 * inner_product(r_ij, p_i.vel - p_j.vel) / r;
                    if(v_sig > v_sig_max) {
                        v_sig_max = v_sig;
                    }
                }
            }
        } else {
            // ===== MASS-BASED APPROACH (Standard GSPH) =====
            real dh_dens_i = 0.0;

            for(int n = 0; n < n_neighbor_tmp; ++n) {
                int const j = neighbor_list[n];
                auto & p_j = particles[j];
                const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
                const real r = std::abs(r_ij);

                if(r >= p_i.sml) {
                    break;
                }

                ++n_neighbor;
                dens_i += p_j.mass * kernel->w(r, p_i.sml);
                dh_dens_i += p_j.mass * kernel->dhw(r, p_i.sml);

                if(i != j) {
                    const real v_sig = p_i.sound + p_j.sound - 3.0 * inner_product(r_ij, p_i.vel - p_j.vel) / r;
                    if(v_sig > v_sig_max) {
                        v_sig_max = v_sig;
                    }
                }
            }

            // Grad-h correction for mass-based: Ω = 1 / (1 + (h/Dρ) * dρ/dh)
            if(m_use_gradh) {
                p_i.gradh = 1.0 / (1.0 + p_i.sml / (DIM * dens_i) * dh_dens_i);
            } else {
                p_i.gradh = 1.0;
            }
        }

        p_i.dens = dens_i;
        p_i.pres = (m_gamma - 1.0) * dens_i * p_i.ene;
        p_i.neighbor = n_neighbor;

        const real h_per_v_sig_i = p_i.sml / v_sig_max;
        if(h_per_v_sig.get() > h_per_v_sig_i) {
            h_per_v_sig.get() = h_per_v_sig_i;
        }

        // Gradient calculation for MUSCL
        if(!m_is_2nd_order) {
            continue;
        }

        vec_t dd, du;
        vec_t dv[DIM];

        if(m_use_volume_based) {
            // Volume-based gradient: Σ V_j (A_j - A_i) ∇W
            for(int n = 0; n < n_neighbor_tmp; ++n) {
                int const j = neighbor_list[n];
                if(i == j) continue;

                auto & p_j = particles[j];
                const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
                const real r = std::abs(r_ij);

                if(r >= p_i.sml * 2.0) continue;

                const vec_t dw_ij = kernel->dw(r_ij, r, p_i.sml);
                const real Vp_j = p_j.nu;

                dd += dw_ij * Vp_j;
                du += dw_ij * (Vp_j * (p_j.ene - p_i.ene));
                for(int k = 0; k < DIM; ++k) {
                    dv[k] += dw_ij * (Vp_j * (p_j.vel[k] - p_i.vel[k]));
                }
            }
        } else {
            // Mass-based gradient: Σ m_j (A_j - A_i) ∇W / ρ
            for(int n = 0; n < n_neighbor; ++n) {
                int const j = neighbor_list[n];
                auto & p_j = particles[j];
                const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
                const real r = std::abs(r_ij);
                const vec_t dw_ij = kernel->dw(r_ij, r, p_i.sml);
                dd += dw_ij * p_j.mass;
                du += dw_ij * (p_j.mass * (p_j.ene - p_i.ene));
                for(int k = 0; k < DIM; ++k) {
                    dv[k] += dw_ij * (p_j.mass * (p_j.vel[k] - p_i.vel[k]));
                }
            }
        }

        grad_d[i] = dd;
        grad_p[i] = (dd * p_i.ene + du) * (m_gamma - 1.0);
        const real rho_inv = 1.0 / p_i.dens;
        for(int k = 0; k < DIM; ++k) {
            grad_v[k][i] = dv[k] * rho_inv;
        }
    }

    // After first calculation, enable volume-based h iteration
    m_first_calculation = false;

    sim->set_h_per_v_sig(h_per_v_sig.min());

    tree->set_kernel();
}

}
}
