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
    m_use_gradh = param->gsph.use_gradh;  // Read grad-h flag (also used for SSPH)
    m_gamma = param->physics.gamma;
    m_neighbor_number = param->physics.neighbor_number;
    m_c_smooth = param->physics.c_smooth;  // C_smooth: smoothing length expansion factor
    m_iteration = param->iterative_sml;
    m_preserve_initial_density = param->preserve_initial_density;
    if(m_iteration) {
        m_kernel_ratio = 1.2;  // Standard value for iterative smoothing
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

        // neighbor search
        const int n_neighbor_tmp = tree->neighbor_search(p_i, neighbor_list, particles, false);
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
        // Grad-h correction factor: Ω_i = 1 / (1 + (h/Dρ) * dρ/dh)
        // This corrects for variable smoothing length in the kernel gradient
        // When disabled (use_gradh=false), set to 1.0 which causes core collapse in hydrostatic tests
        if(m_use_gradh) {
            p_i.gradh = 1.0 / (1.0 + p_i.sml / (DIM * dens_i) * dh_dens_i);
        } else {
            p_i.gradh = 1.0;  // No grad-h correction
        }
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

    tree->set_kernel();
}

void PreInteraction::initial_smoothing(std::shared_ptr<Simulation> sim)
{
    std::cout << "\n=== INITIAL_SMOOTHING CALLED ===" << std::endl;

    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();
    auto * tree = sim->get_tree().get();

    // Count real particles and ghost particles
    int n_real = 0, n_ghost = 0;
    for (int i = 0; i < num; ++i) {
        if (particles[i].is_ghost) n_ghost++;
        else n_real++;
    }
    std::cout << "Total particles: " << num << " (real: " << n_real << ", ghost: " << n_ghost << ")" << std::endl;

    std::cout << "First particle before smoothing: dens=" << particles[0].dens
             << ", pres=" << particles[0].pres << ", mass=" << particles[0].mass << std::endl;

    // Collect debug info for boundary particles
    struct BoundaryDebugInfo {
        int id;
        real pos_x;
        real initial_dens;
        real kernel_sum_dens;
        int n_neighbor;
        int n_ghost_neighbor;
        real sml;
    };
    std::vector<BoundaryDebugInfo> boundary_debug;

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];

        // Ghost particles: skip density computation, keep initial values
        // Their properties are set by initialization and maintained by update_ghost_particles()
        if (p_i.is_ghost) {
            p_i.gradh = 1.0;  // Set gradh=1.0 so forces don't vanish at boundaries
            continue;
        }

        const vec_t & pos_i = p_i.pos;
        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

        // guess smoothing length
        constexpr real A = DIM == 1 ? 2.0 :
                           DIM == 2 ? M_PI :
                                      4.0 * M_PI / 3.0;
        p_i.sml = std::pow(m_neighbor_number * p_i.mass / (p_i.dens * A), 1.0 / DIM);

        // neighbor search
        int const n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, false);

        // Compute density via SPH kernel sum: ρ = Σ m_j * W(r_ij, h)
        real dens_i = 0.0;
        real dh_dens_i = 0.0;
        int n_used = 0;
        int n_ghost_used = 0;
        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
            const real r = std::abs(r_ij);

            if(r >= p_i.sml) {
                break;
            }

            ++n_used;
            if (p_j.is_ghost) ++n_ghost_used;
            dens_i += p_j.mass * kernel->w(r, p_i.sml);
            dh_dens_i += p_j.mass * kernel->dhw(r, p_i.sml);
        }

        // Debug: check if this is a boundary particle (near x=0 or x=0.5)
        real pos_x = pos_i[0];
        bool is_left_boundary = (pos_x < 0.05);
        bool is_right_boundary = (pos_x > 0.45);
        bool is_boundary = is_left_boundary || is_right_boundary;

        if (is_boundary) {
            #pragma omp critical
            {
                BoundaryDebugInfo info;
                info.id = p_i.id;
                info.pos_x = pos_x;
                info.initial_dens = p_i.dens;
                info.kernel_sum_dens = dens_i;
                info.n_neighbor = n_used;
                info.n_ghost_neighbor = n_ghost_used;
                info.sml = p_i.sml;
                boundary_debug.push_back(info);
            }
        }

        // Store computed density (unless preserveInitialDensity is set for shock tubes)
        if(!m_preserve_initial_density) {
            p_i.dens = dens_i;
            p_i.pres = (m_gamma - 1.0) * dens_i * p_i.ene;
        }
        // If preserveInitialDensity is set, keep the initial density/pressure values

        // Grad-h correction factor (use kernel-computed density for gradh even if preserving initial)
        real dens_for_gradh = m_preserve_initial_density ? p_i.dens : dens_i;
        if(m_use_gradh) {
            p_i.gradh = 1.0 / (1.0 + p_i.sml / (DIM * dens_for_gradh) * dh_dens_i);
        } else {
            p_i.gradh = 1.0;
        }
        p_i.neighbor = n_used;
    }

    // Print boundary particle debug info
    std::cout << "\n=== BOUNDARY PARTICLE DENSITY DEBUG ===" << std::endl;
    std::sort(boundary_debug.begin(), boundary_debug.end(),
              [](const BoundaryDebugInfo& a, const BoundaryDebugInfo& b) { return a.pos_x < b.pos_x; });

    int count = 0;
    for (const auto& info : boundary_debug) {
        if (count < 10 || (info.pos_x > 0.45 && count < 20)) {  // First 10 left + first 10 right
            std::cout << "Particle " << info.id << ": x=" << info.pos_x
                     << ", initial_dens=" << info.initial_dens
                     << ", kernel_sum=" << info.kernel_sum_dens
                     << ", ratio=" << (info.kernel_sum_dens / info.initial_dens)
                     << ", n_neighbor=" << info.n_neighbor
                     << ", n_ghost=" << info.n_ghost_neighbor
                     << ", sml=" << info.sml
                     << std::endl;
        }
        count++;
    }
    std::cout << "=== END BOUNDARY DEBUG ===" << std::endl;

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

    // C_smooth approach: use expanded kernel W(r, C_smooth * h) for h-adaptation
    // This makes h vary more smoothly in space, reducing ∇h terms
    // When C_smooth = 1 (default), behavior is unchanged
    // When C_smooth = 2 (typical), h adapts more gradually
    const real cs = m_c_smooth;

    // f = rho h^d - b  where rho is computed with W(r, C_smooth * h)
    // f' = drho/dh h^d + d rho h^{d-1}

    constexpr real epsilon = 1e-4;
    constexpr int max_iter = 10;
    const auto & r_i = p_i.pos;
    for(int i = 0; i < max_iter; ++i) {
        const real h_b = h_i;
        const real h_expanded = cs * h_i;  // Expanded smoothing length

        real dens = 0.0;
        real ddens = 0.0;
        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
            const real r = std::abs(r_ij);

            // Use expanded kernel for h-adaptation (C_smooth effect)
            if(r >= h_expanded) {
                break;
            }

            dens += p_j.mass * kernel->w(r, h_expanded);
            // Chain rule: d/dh W(r, C_smooth*h) = C_smooth * dhW(r, C_smooth*h)
            ddens += p_j.mass * cs * kernel->dhw(r, h_expanded);
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
