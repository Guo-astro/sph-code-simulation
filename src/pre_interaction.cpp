#include <algorithm>
#include <iostream>
#include <cmath>

#include "parameters.hpp"
#include "pre_interaction.hpp"
#include "simulation.hpp"
#include "periodic.hpp"
#include "openmp.hpp"
#include "kernel/kernel_function.hpp"
#include "exception.hpp"
#include "bhtree.hpp"

#ifdef EXHAUSTIVE_SEARCH
#include "exhaustive_search.hpp"
#endif

namespace sph
{

void PreInteraction::initialize(std::shared_ptr<SPHParameters> param)
{
    m_use_time_dependent_av = param->av.use_time_dependent_av;
    if(m_use_time_dependent_av) {
        m_alpha_max = param->av.alpha_max;
        m_alpha_min = param->av.alpha_min;
        m_epsilon = param->av.epsilon;
    }
    m_use_balsara_switch = param->av.use_balsara_switch;
    m_gamma = param->physics.gamma;
    m_neighbor_number = param->physics.neighbor_number;
    m_iteration = param->iterative_sml;
    if(m_iteration) {
        m_kernel_ratio = 1.2;
    } else {
        m_kernel_ratio = 1.0;
    }
    m_first = true;
    
    // Jeans length resolution check parameters
    m_gravity_enabled = param->gravity.is_valid;
    m_G = param->gravity.constant;
    m_jeans_warning_issued = false;
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
    const real dt = sim->get_dt();
    auto * tree = sim->get_tree().get();

    omp_real h_per_v_sig(std::numeric_limits<real>::max());

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

        // guess smoothing length
        constexpr real A = DIM == 1 ? 2.0 :
                           DIM == 2 ? M_PI :
                                      4.0 * M_PI / 3.0;
        p_i.sml = std::pow(m_neighbor_number * p_i.mass / (p_i.dens * A), 1.0 / DIM) * m_kernel_ratio;
        
        // neighbor search
#ifdef EXHAUSTIVE_SEARCH
        const int n_neighbor_tmp = exhaustive_search(p_i, p_i.sml, particles, num, neighbor_list, m_neighbor_number * neighbor_list_size, periodic, false);
#else
        const int n_neighbor_tmp = tree->neighbor_search(p_i, neighbor_list, particles, false);
#endif
        // smoothing length
        if(m_iteration) {
            p_i.sml = newton_raphson(p_i, particles, neighbor_list, n_neighbor_tmp, periodic, kernel);
        }

        // density etc.
        real dens_i = 0.0;
        real dh_dens_i = 0.0;
        real v_sig_max = p_i.sound * 2.0;
        const vec_t & pos_i = p_i.pos;
        int n_neighbor = 0;
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

        p_i.dens = dens_i;
        p_i.pres = (m_gamma - 1.0) * dens_i * p_i.ene;
        p_i.gradh = 1.0 / (1.0 + p_i.sml / (DIM * dens_i) * dh_dens_i);
        p_i.neighbor = n_neighbor;

        const real h_per_v_sig_i = p_i.sml / v_sig_max;
        if(h_per_v_sig.get() > h_per_v_sig_i) {
            h_per_v_sig.get() = h_per_v_sig_i;
        }

        // Artificial viscosity
        if(m_use_balsara_switch && DIM != 1) {
#if DIM != 1
            // balsara switch
            real div_v = 0.0;
#if DIM == 2
            real rot_v = 0.0;
#else
            vec_t rot_v = 0.0;
#endif
            for(int n = 0; n < n_neighbor; ++n) {
                int const j = neighbor_list[n];
                auto & p_j = particles[j];
                const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
                const real r = std::abs(r_ij);
                const vec_t dw = kernel->dw(r_ij, r, p_i.sml);
                const vec_t v_ij = p_i.vel - p_j.vel;
                div_v -= p_j.mass * inner_product(v_ij, dw);
                rot_v += vector_product(v_ij, dw) * p_j.mass;
            }
            div_v /= p_i.dens;
            rot_v /= p_i.dens;
            p_i.balsara = std::abs(div_v) / (std::abs(div_v) + std::abs(rot_v) + 1e-4 * p_i.sound / p_i.sml);

            // time dependent alpha
            if(m_use_time_dependent_av) {
                const real tau_inv = m_epsilon * p_i.sound / p_i.sml;
                const real dalpha = (-(p_i.alpha - m_alpha_min) * tau_inv + std::max(-div_v, (real)0.0) * (m_alpha_max - p_i.alpha)) * dt;
                p_i.alpha += dalpha;
            }
#endif
        } else if(m_use_time_dependent_av) {
            real div_v = 0.0;
            for(int n = 0; n < n_neighbor; ++n) {
                int const j = neighbor_list[n];
                auto & p_j = particles[j];
                const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
                const real r = std::abs(r_ij);
                const vec_t dw = kernel->dw(r_ij, r, p_i.sml);
                const vec_t v_ij = p_i.vel - p_j.vel;
                div_v -= p_j.mass * inner_product(v_ij, dw);
            }
            div_v /= p_i.dens;
            const real tau_inv = m_epsilon * p_i.sound / p_i.sml;
            const real s_i = std::max(-div_v, (real)0.0);
            p_i.alpha = (p_i.alpha + dt * tau_inv * m_alpha_min + s_i * dt * m_alpha_max) / (1.0 + dt * tau_inv + s_i * dt);
        }
    }

    sim->set_h_per_v_sig(h_per_v_sig.min());

    // ========================================================================
    // JEANS LENGTH RESOLUTION CHECK (Truelove et al. 1997)
    // ========================================================================
    // For self-gravitating systems, the Jeans length λ_J = c_s * sqrt(π / (G ρ))
    // must be resolved to avoid artificial fragmentation.
    // Criterion: h / λ_J < 0.25 (Truelove criterion)
    // ========================================================================
    if(m_gravity_enabled && m_G > 0.0) {
        real max_jeans_ratio = 0.0;
        int worst_particle_id = -1;
        real worst_h = 0.0;
        real worst_lambda_J = 0.0;
        real worst_rho = 0.0;
        real worst_cs = 0.0;
        
        for(int i = 0; i < num; ++i) {
            const auto& p_i = particles[i];
            if(p_i.dens > 0.0 && p_i.sound > 0.0) {
                // Jeans length: λ_J = c_s * sqrt(π / (G ρ))
                const real lambda_J = p_i.sound * std::sqrt(M_PI / (m_G * p_i.dens));
                const real jeans_ratio = p_i.sml / lambda_J;
                
                if(jeans_ratio > max_jeans_ratio) {
                    max_jeans_ratio = jeans_ratio;
                    worst_particle_id = i;
                    worst_h = p_i.sml;
                    worst_lambda_J = lambda_J;
                    worst_rho = p_i.dens;
                    worst_cs = p_i.sound;
                }
            }
        }
        
        // Truelove criterion: h/λ_J < 0.25 to avoid artificial fragmentation
        constexpr real TRUELOVE_THRESHOLD = 0.25;
        if(max_jeans_ratio > TRUELOVE_THRESHOLD && !m_jeans_warning_issued) {
            std::cerr << "\n*** JEANS LENGTH RESOLUTION WARNING ***" << std::endl;
            std::cerr << "Truelove criterion violated: h/λ_J = " << max_jeans_ratio 
                      << " > " << TRUELOVE_THRESHOLD << std::endl;
            std::cerr << "Worst particle: id=" << worst_particle_id << std::endl;
            std::cerr << "  h = " << worst_h << ", λ_J = " << worst_lambda_J << std::endl;
            std::cerr << "  ρ = " << worst_rho << ", c_s = " << worst_cs << std::endl;
            std::cerr << "Risk: Artificial gravitational fragmentation may occur!" << std::endl;
            std::cerr << "Suggestion: Increase resolution (more particles) or use pressure floor." << std::endl;
            std::cerr << "**********************************************\n" << std::endl;
            m_jeans_warning_issued = true;
        } else if(max_jeans_ratio <= TRUELOVE_THRESHOLD && m_jeans_warning_issued) {
            // Reset warning flag if condition is no longer violated
            m_jeans_warning_issued = false;
        }
    }

#ifndef EXHAUSTIVE_SEARCH
    tree->set_kernel();
#endif
}

void PreInteraction::initial_smoothing(std::shared_ptr<Simulation> sim)
{
    std::cout << "\n=== INITIAL_SMOOTHING CALLED ===" << std::endl;
    
    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();
    auto * tree = sim->get_tree().get();
    
    std::cout << "First particle before smoothing: dens=" << particles[0].dens 
             << ", pres=" << particles[0].pres << ", mass=" << particles[0].mass << std::endl;

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        const vec_t & pos_i = p_i.pos;
        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

        // guess smoothing length
        constexpr real A = DIM == 1 ? 2.0 :
                           DIM == 2 ? M_PI :
                                      4.0 * M_PI / 3.0;
        p_i.sml = std::pow(m_neighbor_number * p_i.mass / (p_i.dens * A), 1.0 / DIM);
        
        // DEBUG first particle
        bool is_first = (i == 0);
        if (is_first) {
            #pragma omp critical
            {
                std::cout << "Particle 0: initial sml guess = " << p_i.sml 
                         << " (mass=" << p_i.mass << ", dens=" << p_i.dens << ")" << std::endl;
            }
        }
        
        // neighbor search
#ifdef EXHAUSTIVE_SEARCH
        int const n_neighbor = exhaustive_search(p_i, p_i.sml, particles, num, neighbor_list, m_neighbor_number * neighbor_list_size, periodic, false);
#else
        int const n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, false);
#endif

        if (is_first) {
            #pragma omp critical
            {
                std::cout << "Particle 0: found " << n_neighbor << " neighbors" << std::endl;
            }
        }

        // density
        real dens_i = 0.0;
        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
            const real r = std::abs(r_ij);

            if(r >= p_i.sml) {
                break;
            }

            dens_i += p_j.mass * kernel->w(r, p_i.sml);
            
            if (is_first && n < 5) {
                #pragma omp critical
                {
                    std::cout << "  neighbor " << n << ": j=" << j << ", r=" << r 
                             << ", mass=" << p_j.mass << ", w=" << kernel->w(r, p_i.sml) << std::endl;
                }
            }
        }

        p_i.dens = dens_i;
        
        if (is_first) {
            #pragma omp critical
            {
                std::cout << "Particle 0: final dens = " << dens_i << std::endl;
            }
        }
    }
    
    std::cout << "First particle after smoothing: dens=" << particles[0].dens 
             << ", pres=" << particles[0].pres << ", mass=" << particles[0].mass << std::endl;
    std::cout << "=== INITIAL_SMOOTHING COMPLETE ===" << std::endl;
}

inline real powh_(const real h) {
#if DIM == 1
    return 1;
#elif DIM == 2
    return h;
#elif DIM == 3
    return h * h;
#endif
}

real PreInteraction::newton_raphson(
    const SPHParticle & p_i,
    const std::vector<SPHParticle> & particles,
    const std::vector<int> & neighbor_list,
    const int n_neighbor,
    const Periodic * periodic,
    const KernelFunction * kernel
)
{
    real h_i = p_i.sml / m_kernel_ratio;
    constexpr real A = DIM == 1 ? 2.0 :
                       DIM == 2 ? M_PI :
                                  4.0 * M_PI / 3.0;
    const real b = p_i.mass * m_neighbor_number / A;

    // f = rho h^d - b
    // f' = drho/dh h^d + d rho h^{d-1}

    constexpr real epsilon = 1e-4;
    constexpr int max_iter = 10;
    const auto & r_i = p_i.pos;
    for(int i = 0; i < max_iter; ++i) {
        const real h_b = h_i;

        real dens = 0.0;
        real ddens = 0.0;
        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
            const real r = std::abs(r_ij);

            if(r >= h_i) {
                break;
            }

            dens += p_j.mass * kernel->w(r, h_i);
            ddens += p_j.mass * kernel->dhw(r, h_i);
        }

        const real f = dens * powh(h_i) - b;
        const real df = ddens * powh(h_i) + DIM * dens * powh_(h_i);

        h_i -= f / df;

        if(std::abs(h_i - h_b) < (h_i + h_b) * epsilon) {
            return h_i;
        }
    }

#pragma omp critical
    {
        WRITE_LOG << "Particle id " << p_i.id << " is not convergence";
    }

    return p_i.sml / m_kernel_ratio;
}

}
