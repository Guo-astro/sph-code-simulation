#include "srgsph/sr_fluid_force.hpp"
#include "srgsph/sr_exact_riemann.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "periodic.hpp"
#include "openmp.hpp"
#include "kernel/kernel_function.hpp"
#include "bhtree.hpp"
#include <array>
#include <cmath>
#include <algorithm>

#ifdef EXHAUSTIVE_SEARCH
#include "exhaustive_search.hpp"
#endif

namespace sph
{
namespace srgsph
{

void FluidForce::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::FluidForce::initialize(param);

    m_gamma = param->physics.gamma;
    m_c_speed = param->srgsph.c_speed;
    m_c_shock = (param->srgsph.c_shock > 0.0) ? param->srgsph.c_shock : 3.0;
    m_c_cd = (param->srgsph.c_cd > 0.0) ? param->srgsph.c_cd : 0.2;
    m_use_muscl = param->srgsph.is_2nd_order;
}

void FluidForce::solve_interface_state(
    const real left_state[5],
    const real right_state[5],
    real & pstar,
    real & vstar,
    real & vt_star) const
{
    riemann::RiemannState left{left_state[0], left_state[1], left_state[2], left_state[3]};
    riemann::RiemannState right{right_state[0], right_state[1], right_state[2], right_state[3]};
    const real vt_left = left_state[4];
    const real vt_right = right_state[4];

    const bool converged = riemann::exact_riemann_solver(
        left, right, vt_left, vt_right, m_gamma, m_c_speed,
        pstar, vstar, vt_star, 100, 1e-10);

    if (!converged) {
        riemann::hllc_riemann_solver(
            left, right, vt_left, vt_right, m_gamma, m_c_speed,
            pstar, vstar, vt_star);
    }
}

/**
 * Main calculation loop for fluid forces
 * Based on Python compute_forces() (lines 284-597)
 */
void FluidForce::calculation(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();
    auto * tree = sim->get_tree().get();

    std::vector<vec_t> * grad_p = nullptr;
    std::vector<vec_t> * grad_rho = nullptr;
    std::array<std::vector<vec_t> *, DIM> grad_v{};
    grad_v.fill(nullptr);
    bool gradients_available = false;

    if (m_use_muscl) {
        try {
            grad_p = &sim->get_vector_array("grad_pressure");
            grad_rho = &sim->get_vector_array("grad_density");
            grad_v[0] = &sim->get_vector_array("grad_velocity_0");
#if DIM >= 2
            grad_v[1] = &sim->get_vector_array("grad_velocity_1");
#endif
#if DIM == 3
            grad_v[2] = &sim->get_vector_array("grad_velocity_2");
#endif

            gradients_available = (grad_p != nullptr && grad_rho != nullptr);
            for (int dim = 0; gradients_available && dim < DIM; ++dim) {
                gradients_available = grad_v[dim] != nullptr;
            }
        } catch (...) {
            grad_p = nullptr;
            grad_rho = nullptr;
            grad_v.fill(nullptr);
            gradients_available = false;
        }
    }

    // Reset derivatives
#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        particles[i].dS = vec_t(0.0);
        particles[i].de = 0.0;
    }

    const real sqrt_two = std::sqrt(2.0);

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

#ifdef EXHAUSTIVE_SEARCH
        const int n_neighbor = exhaustive_search(p_i, p_i.sml * 6.0, particles, num,
                                                  neighbor_list, m_neighbor_number * neighbor_list_size,
                                                  periodic, false);
#else
        const int n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, false);
#endif

        for (int n = 0; n < n_neighbor; ++n) {
            const int j = neighbor_list[n];
            if (j == i) {
                continue;
            }

            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(p_i.pos, p_j.pos);
            const real r = std::abs(r_ij);

            if (r > 3.0 * std::max(p_i.sml, p_j.sml)) {
                continue;
            }

            const vec_t n_ij = r_ij * (-1.0 / r);

            const vec_t midpoint = (p_i.pos + p_j.pos) * 0.5;
            const vec_t dr_i = midpoint - p_i.pos;
            const vec_t dr_j = midpoint - p_j.pos;

            bool use_second_order = gradients_available;
            if (use_second_order) {
                const vec_t e_ij = n_ij * (-1.0);
                const vec_t dv = p_i.vel - p_j.vel;
                const real shock_indicator = m_c_shock * inner_product(e_ij, dv);
                const real min_cs = std::min(p_i.sound, p_j.sound);
                const bool is_shock = (-shock_indicator) > min_cs;

                bool is_cd = false;
                if (p_i.pres > 0.0 && p_j.pres > 0.0) {
                    const real log_p = std::abs(std::log10(p_i.pres / p_j.pres));
                    is_cd = log_p > m_c_cd;
                } else {
                    is_cd = true;
                }

                if (is_shock || is_cd) {
                    use_second_order = false;
                }
            }

            real p_L = p_i.pres;
            real rho_L = p_i.dens;
            real vx_L = p_i.vel[0];
#if DIM >= 2
            real vy_L = p_i.vel[1];
#endif
#if DIM == 3
            real vz_L = p_i.vel[2];
#endif

            real p_R = p_j.pres;
            real rho_R = p_j.dens;
            real vx_R = p_j.vel[0];
#if DIM >= 2
            real vy_R = p_j.vel[1];
#endif
#if DIM == 3
            real vz_R = p_j.vel[2];
#endif

            if (use_second_order) {
                const real p_L_trial = p_i.pres + inner_product((*grad_p)[i], dr_i);
                const real rho_L_trial = p_i.dens + inner_product((*grad_rho)[i], dr_i);
                const real vx_L_trial = p_i.vel[0] + inner_product((*grad_v[0])[i], dr_i);
#if DIM >= 2
                const real vy_L_trial = p_i.vel[1] + inner_product((*grad_v[1])[i], dr_i);
#endif
#if DIM == 3
                const real vz_L_trial = p_i.vel[2] + inner_product((*grad_v[2])[i], dr_i);
#endif

                const real p_R_trial = p_j.pres + inner_product((*grad_p)[j], dr_j);
                const real rho_R_trial = p_j.dens + inner_product((*grad_rho)[j], dr_j);
                const real vx_R_trial = p_j.vel[0] + inner_product((*grad_v[0])[j], dr_j);
#if DIM >= 2
                const real vy_R_trial = p_j.vel[1] + inner_product((*grad_v[1])[j], dr_j);
#endif
#if DIM == 3
                const real vz_R_trial = p_j.vel[2] + inner_product((*grad_v[2])[j], dr_j);
#endif

                const bool valid_left = (p_L_trial > 0.0 && rho_L_trial > 0.0);
                const bool valid_right = (p_R_trial > 0.0 && rho_R_trial > 0.0);

                if (valid_left && valid_right) {
                    p_L = p_L_trial;
                    rho_L = rho_L_trial;
                    vx_L = vx_L_trial;
#if DIM >= 2
                    vy_L = vy_L_trial;
#endif
#if DIM == 3
                    vz_L = vz_L_trial;
#endif
                    p_R = p_R_trial;
                    rho_R = rho_R_trial;
                    vx_R = vx_R_trial;
#if DIM >= 2
                    vy_R = vy_R_trial;
#endif
#if DIM == 3
                    vz_R = vz_R_trial;
#endif
                } else {
                    use_second_order = false;
                }
            }

            if (!use_second_order) {
                p_L = p_i.pres;
                rho_L = p_i.dens;
                vx_L = p_i.vel[0];
#if DIM >= 2
                vy_L = p_i.vel[1];
#endif
#if DIM == 3
                vz_L = p_i.vel[2];
#endif
                p_R = p_j.pres;
                rho_R = p_j.dens;
                vx_R = p_j.vel[0];
#if DIM >= 2
                vy_R = p_j.vel[1];
#endif
#if DIM == 3
                vz_R = p_j.vel[2];
#endif
            }

#if DIM == 1
            const vec_t v_L(vx_L);
            const vec_t v_R(vx_R);
#elif DIM == 2
            const vec_t v_L(vx_L, vy_L);
            const vec_t v_R(vx_R, vy_R);
#else
            const vec_t v_L(vx_L, vy_L, vz_L);
            const vec_t v_R(vx_R, vy_R, vz_R);
#endif

            const real vn_L = inner_product(v_L, n_ij);
            const real vn_R = inner_product(v_R, n_ij);
            const vec_t vt_vec_L = v_L - n_ij * vn_L;
            const vec_t vt_vec_R = v_R - n_ij * vn_R;
            const real vt_L = std::abs(vt_vec_L);
            const real vt_R = std::abs(vt_vec_R);

            real left_state[5] = {vn_L, rho_L, p_L, p_i.sound, vt_L};
            real right_state[5] = {vn_R, rho_R, p_R, p_j.sound, vt_R};

            real P_star = 0.0;
            real v_star_n = 0.0;
            real vt_star = 0.0;
            solve_interface_state(left_state, right_state, P_star, v_star_n, vt_star);

            vec_t v_star_vec = n_ij * v_star_n;
            if (vt_star > 1e-10) {
                if (v_star_n > 0.0 && vt_L > 1e-10) {
                    v_star_vec = v_star_vec + vt_vec_L * (vt_star / vt_L);
                } else if (v_star_n <= 0.0 && vt_R > 1e-10) {
                    v_star_vec = v_star_vec + vt_vec_R * (vt_star / vt_R);
                }
            }

            const real V_i = p_i.nu / p_i.N;
            const real V_j = p_j.nu / p_j.N;
            const real Vij2 = 0.5 * (V_i * V_i + V_j * V_j);

            const vec_t grad_W_i = kernel->dw(r_ij, r, sqrt_two * p_i.sml);
            const vec_t grad_W_j = kernel->dw(r_ij, r, sqrt_two * p_j.sml);
            const vec_t term_grad = grad_W_i - (grad_W_j * (-1.0));

            const vec_t force = term_grad * (-P_star * Vij2);
            const real power = inner_product(v_star_vec, force);

            p_i.dS[0] += force[0];
#if DIM >= 2
            p_i.dS[1] += force[1];
#endif
#if DIM == 3
            p_i.dS[2] += force[2];
#endif
            p_i.de += power;
        }
    }

#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        FluidForce::normalize_sr_derivatives(particles[i]);
    }
}

} // namespace srgsph
} // namespace sph
