#include "defines.hpp"
#include "particle.hpp"
#include "periodic.hpp"
#include "simulation.hpp"
#include "bhtree.hpp"
#include "kernel/kernel_function.hpp"
#include "srgsph/sr_pre_interaction.hpp"
#include "srgsph/sr_primitive_recovery.hpp"
#include <iostream>

#ifdef EXHAUSTIVE_SEARCH
#include "exhaustive_search.hpp"
#endif

namespace sph
{
namespace srgsph
{

void PreInteraction::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::PreInteraction::initialize(param);
    m_eta = param->srgsph.eta;
    m_c_smooth = param->srgsph.c_smooth;
    m_c_speed = param->srgsph.c_speed;
    m_gamma = param->physics.gamma;
    m_iteration = param->iterative_sml;
    m_first = true;
}

void PreInteraction::initial_smoothing(std::shared_ptr<Simulation> sim)
{
    // Initial guess for smoothing length based on particle spacing
    auto & particles = sim->get_particles();
    const int num = sim->get_particle_num();
    
    // Estimate from neighbor number and particle mass
    // Similar to standard SPH initial guess
    constexpr real A = DIM == 1 ? 2.0 :
                       DIM == 2 ? M_PI :
                                  4.0 * M_PI / 3.0;
    const real neighbor = m_neighbor_number;
    
    std::cerr << "\n=== INITIAL SMOOTHING LENGTH CALCULATION ===" << std::endl;
    std::cerr << "Formula: h = η × (ν / N))^(1/d)" << std::endl;
    std::cerr << "  η (eta) = " << m_eta << std::endl;
    std::cerr << "  d (dimension) = " << DIM << std::endl;
    
#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        // Initial h using eta parameter directly
        // h = η × (ν/N)^(1/d)
        p_i.sml = m_eta * std::pow(p_i.nu / p_i.N, 1.0 / DIM);
        
        // Log a few representative particles
        if(i == 0 || i == 1600 || i == 3200 || i == 3599) {
            #pragma omp critical
            {
                std::cerr << "Particle " << i << " (x=" << p_i.pos[0] << "):" << std::endl;
                std::cerr << "  ν = " << p_i.nu << std::endl;
                std::cerr << "  N = " << p_i.N << " (baryon number density from initialization)" << std::endl;
                std::cerr << "  h_initial = η × (ν/N)^(1/d) = " << m_eta << " × (" << p_i.nu << "/" << p_i.N << ")^(1/" << DIM << ")" << std::endl;
                std::cerr << "  h_initial = " << p_i.sml << std::endl;
            }
        }
    }
    std::cerr << "=== END INITIAL SMOOTHING ===" << std::endl << std::endl;
}

real PreInteraction::compute_volume(
    const SPHParticle & p_i,
    const std::vector<SPHParticle> & particles,
    const std::vector<int> & neighbor_list,
    const int n_neighbor,
    const Periodic * periodic,
    const KernelFunction * kernel,
    const real h
)
{
    // Eq. 33: Vp(x) = [Σ_j W(x-x_j, h)]^(-1)
    real sum_w = 0.0;
    const vec_t & r_i = p_i.pos;
    
    for(int n = 0; n < n_neighbor; ++n) {
        int const j = neighbor_list[n];
        auto & p_j = particles[j];
        const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
        const real r = std::abs(r_ij);
        
        if(r < h) {
            sum_w += kernel->w(r, h);
        }
    }
    
    // Self-contribution
    sum_w += kernel->w(0.0, h);
    
    // Volume = 1 / sum_w
    return 1.0 / sum_w;
}

real PreInteraction::compute_smoothing_length(
    const SPHParticle & p_i,
    const std::vector<SPHParticle> & particles,
    const std::vector<int> & neighbor_list,
    const int n_neighbor,
    const Periodic * periodic,
    const KernelFunction * kernel
)
{
    // VOLUME-BASED APPROACH with gather method (Eqs. 35-36 in paper)
    // h = η * Vp*^(1/d)
    // where Vp* uses kernel with C_smooth * h
    //
    // This approach maintains smooth variation of h across contact discontinuities
    // and is numerically more stable than the standard approach.
    //
    // CRITICAL: η must relate to neighbor_number to maintain correct h scale
    // For uniform 1D: Vp* ≈ particle spacing dx
    // We want h ≈ (neighbor_number/2) × dx to contain ~neighbor_number particles
    // Therefore: η = neighbor_number / (2×A) where A is the kernel support factor
    constexpr real A = DIM == 1 ? 2.0 :
                       DIM == 2 ? M_PI :
                                  4.0 * M_PI / 3.0;
    const real eta_corrected = m_eta;  // Use eta from config, not hardcoded neighbor/2A
    
    real h = p_i.sml;  // Initial guess from previous step
    const int max_iter = 20;
    const real tol = 1.0e-6;
    
    // Diagnostic: track a few representative particles
    static int call_count = 0;
    static bool detailed_log = false;
    if(call_count < 3) {
        detailed_log = true;
        call_count++;
    }
    
    const int particle_id = p_i.id;
    const bool is_tracked = detailed_log && (particle_id == 0 || particle_id == 1600 || particle_id == 3200 || particle_id == 3599);
    
    if(is_tracked) {
        std::cerr << "\n=== SMOOTHING LENGTH ITERATION for particle " << particle_id << " ===" << std::endl;
        std::cerr << "Position: x = " << p_i.pos[0] << std::endl;
        std::cerr << "Baryon number: ν = " << p_i.nu << std::endl;
        std::cerr << "Initial h guess: " << h << std::endl;
        std::cerr << "Parameters: η = " << m_eta << " (from JSON)" << std::endl;
        std::cerr << "           C_smooth = " << m_c_smooth << std::endl;
        std::cerr << "Neighbors found: " << n_neighbor << std::endl;
        
        // Theoretical expectation for uniform 1D:
        // For uniform spacing dx, with neighbor_number particles in range 2h,
        // we expect h ≈ neighbor_number × dx / 2
        // Estimate local spacing from nearest neighbor
        real min_dist = 1e10;
        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            if(j == particle_id) continue;
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(p_i.pos, p_j.pos);
            const real r = std::abs(r_ij);
            if(r > 0 && r < min_dist) min_dist = r;
        }
        const real expected_h = m_neighbor_number * min_dist / 2.0;
        std::cerr << "Nearest neighbor distance: " << min_dist << std::endl;
        std::cerr << "Expected h (50 neighbors): " << expected_h << std::endl;
    }
    
    for(int iter = 0; iter < max_iter; ++iter) {
        // Compute Vp* at current h using C_smooth expansion
        // Vp*(x) = [Σ_j W(x-x_j, C_smooth*h)]^(-1)  (Eq. 36)
        const real h_smooth = m_c_smooth * h;
        real sum_w_star = 0.0;
        const vec_t & r_i = p_i.pos;
        
        int neighbors_in_range = 0;
        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
            const real r = std::abs(r_ij);
            
            if(r < h_smooth) {
                sum_w_star += kernel->w(r, h_smooth);
                neighbors_in_range++;
            }
        }
        
        // Self-contribution
        sum_w_star += kernel->w(0.0, h_smooth);
        
        if(sum_w_star < 1.0e-20) {
            // Pathological case - expand search
            if(is_tracked) {
                std::cerr << "  Iter " << iter << ": sum_w_star too small, expanding h" << std::endl;
            }
            h *= 1.5;
            continue;
        }
        
        // Vp* = 1 / sum_w_star
        const real Vp_star = 1.0 / sum_w_star;
        
        // New smoothing length: h = η * Vp*^(1/d)  (Eq. 35)
        // Use corrected η that relates to neighbor_number
        const real h_new = eta_corrected * std::pow(Vp_star, 1.0 / real(DIM));
        
        if(is_tracked && iter < 10) {
            std::cerr << "  Iter " << iter << ":" << std::endl;
            std::cerr << "    h_current = " << h << ", h_smooth = C×h = " << h_smooth << std::endl;
            std::cerr << "    Neighbors in C×h range: " << neighbors_in_range << std::endl;
            std::cerr << "    Σ W(r, C×h) = " << sum_w_star << std::endl;
            std::cerr << "    Vp* = 1/Σ = " << Vp_star << std::endl;
            std::cerr << "    h_new = η_corrected × (Vp*)^(1/d) = " << eta_corrected << " × " << Vp_star << "^" << (1.0/DIM) << " = " << h_new << std::endl;
            std::cerr << "    |h_new - h| / h = " << std::abs(h_new - h) / h << std::endl;
            std::cerr << "    [Physics: Vp* ≈ particle spacing, η×Vp* should give h for " << m_neighbor_number << " neighbors]" << std::endl;
        }
        
        // Damping to avoid oscillations
        const real h_damped = 0.5 * (h + h_new);
        
        // Check convergence
        if(std::abs(h_damped - h) < tol * h) {
            if(is_tracked) {
                std::cerr << "  CONVERGED at iteration " << iter << std::endl;
                std::cerr << "  Final h = " << h_damped << std::endl;
                std::cerr << "  Ratio to expected: " << h_damped / (m_neighbor_number * Vp_star / 2.0) << std::endl;
            }
            return h_damped;
        }
        
        h = h_damped;
        
        // Safety bounds - be more permissive
        if(h < 1.0e-10) {
            if(is_tracked) {
                std::cerr << "  WARNING: h hit lower bound!" << std::endl;
            }
            h = 1.0e-10;
        }
        if(h > 10.0) {
            if(is_tracked) {
                std::cerr << "  WARNING: h hit upper bound!" << std::endl;
            }
            h = 10.0;
        }
    }
    
    if(is_tracked) {
        std::cerr << "  DID NOT CONVERGE after " << max_iter << " iterations" << std::endl;
        std::cerr << "  Final h = " << h << std::endl;
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

    // For MUSCL reconstruction
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
        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

        // Neighbor search with current smoothing length
#ifdef EXHAUSTIVE_SEARCH
        int const n_neighbor = exhaustive_search(p_i, p_i.sml, particles, num,
                                                neighbor_list, m_neighbor_number * neighbor_list_size,
                                                periodic, false);
#else
        int const n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, false);
#endif

        // 1. Update smoothing length
        if(m_iteration) {
            // Iterative smoothing length (volume-based approach)
            p_i.sml = compute_smoothing_length(p_i, particles, neighbor_list, n_neighbor,
                                              periodic, kernel);
        } else {
            // Fixed smoothing length: h = η * (ν/N)^(1/d)
            // Using simple formula without iteration
            constexpr real A = DIM == 1 ? 2.0 :
                               DIM == 2 ? M_PI :
                                          4.0 * M_PI / 3.0;
            // For uniform distribution: N ≈ 1 initially, so h ≈ η * ν^(1/d)
            // Use the volume from previous step or initial guess
            p_i.sml = m_eta * std::pow(p_i.nu / p_i.N, 1.0 / DIM);
        }
        const real h_i = p_i.sml;

        // 2. Compute particle volume Vp using VOLUME-BASED approach (Eq. 33)
        // Vp(x) = [Σ_j W(x - x_j, h)]^(-1)
        const real Vp_i = compute_volume(p_i, particles, neighbor_list, n_neighbor,
                                         periodic, kernel, h_i);
        
        // 3. Compute baryon number density using VOLUME-BASED approach (Eq. 42)
        // N_volume-based(x) = ν(x) / Vp(x)
        p_i.N = p_i.nu / Vp_i;
        
        // Diagnostic for tracked particles
        static int pre_call_count = 0;
        const bool is_first_few_calls = (pre_call_count < 3);
        const bool is_tracked_particle = (i == 0 || i == 1600 || i == 3200 || i == 3599);
        
        if(is_first_few_calls && is_tracked_particle) {
            std::cerr << "\n--- Particle " << i << " at pre-interaction call " << pre_call_count << " ---" << std::endl;
            std::cerr << "  Position: x = " << p_i.pos[0] << std::endl;
            std::cerr << "  h = " << h_i << std::endl;
            std::cerr << "  ν = " << p_i.nu << std::endl;
            std::cerr << "  Vp = " << Vp_i << " (particle volume)" << std::endl;
            std::cerr << "  N = ν/Vp = " << p_i.N << " (baryon number density)" << std::endl;
            std::cerr << "  Expected N for uniform: ν/dx = " << p_i.nu / Vp_i << std::endl;
        }
        
        if(is_first_few_calls && i == num - 1) {
            pre_call_count++;
        }

        // 4. Recover primitive variables from conserved (S, e, N)
        // This gives us velocity, pressure, density for gradients
        auto prim = PrimitiveRecovery::conserved_to_primitive(
            p_i.S, p_i.e, p_i.N, m_gamma, m_c_speed
        );
        
        p_i.vel = prim.vel;
        p_i.pres = prim.pressure;
        p_i.dens = prim.density;  // Rest frame density n
        p_i.sound = prim.sound_speed;
        p_i.gamma_lor = prim.gamma_lor;
        p_i.enthalpy = prim.enthalpy;

        // 4. Compute gradients for MUSCL reconstruction
        // ∇φ_i = (1/Ω_i) Σ_j (m_j/ρ_j) (φ_j - φ_i) ∇W_ij
        vec_t grad_dens(0.0), grad_pres(0.0);
        vec_t grad_vel[DIM];
        for(int k = 0; k < DIM; ++k) {
            grad_vel[k] = vec_t(0.0);
        }

        real omega = 0.0;
        const vec_t & r_i = p_i.pos;

        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
            const real r = std::abs(r_ij);

            if(r >= h_i || r == 0.0) {
                continue;
            }

            // For gradients, need neighbor's primitive variables too
            auto prim_j = PrimitiveRecovery::conserved_to_primitive(
                p_j.S, p_j.e, p_j.N, m_gamma, m_c_speed
            );

            const vec_t dw = kernel->dw(r_ij, r, h_i);
            const real vol_j = p_j.nu / p_j.N;  // Particle volume
            
            omega += vol_j * kernel->w(r, h_i);

            grad_dens += dw * (vol_j * (prim_j.density - p_i.dens));
            grad_pres += dw * (vol_j * (prim_j.pressure - p_i.pres));
            for(int k = 0; k < DIM; ++k) {
                grad_vel[k] += dw * (vol_j * (prim_j.vel[k] - p_i.vel[k]));
            }
        }

        if(omega > 0.0) {
            const real omega_inv = 1.0 / omega;
            grad_d[i] = grad_dens * omega_inv;
            grad_p[i] = grad_pres * omega_inv;
            for(int k = 0; k < DIM; ++k) {
                grad_v[k][i] = grad_vel[k] * omega_inv;
            }
        } else {
            grad_d[i] = vec_t(0.0);
            grad_p[i] = vec_t(0.0);
            for(int k = 0; k < DIM; ++k) {
                grad_v[k][i] = vec_t(0.0);
            }
        }
    }
}

}
}
