#include "defines.hpp"
#include "particle.hpp"
#include "periodic.hpp"
#include "simulation.hpp"
#include "bhtree.hpp"
#include "kernel/kernel_function.hpp"
#include "srgsph/sr_pre_interaction.hpp"
#include "srgsph/sr_primitive_recovery.hpp"

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
    
#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        // Initial h from h^d = neighbor * mass / density
        // For SR, use baryon number instead of mass
        p_i.sml = std::pow(neighbor * p_i.nu / (A * p_i.N), 1.0 / DIM);
    }
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
    // Iterative solution of Eqs. 35-36:
    // h = η Vp*^(1/d)
    // Vp* = [Σ_j W(x_i-x_j, C_smooth*h)]^(-1)
    
    real h = p_i.sml;  // Initial guess from previous step
    const int max_iter = 10;
    const real tol = 1.0e-6;
    
    for(int iter = 0; iter < max_iter; ++iter) {
        // Compute Vp* at current h
        const real h_smooth = m_c_smooth * h;
        const real Vp_star = compute_volume(p_i, particles, neighbor_list, n_neighbor,
                                           periodic, kernel, h_smooth);
        
        // New h = η Vp*^(1/d)
        const real h_new = m_eta * std::pow(Vp_star, 1.0 / DIM);
        
        // Check convergence
        if(std::abs(h_new - h) < tol * h) {
            return h_new;
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

        // 1. Update smoothing length (iterative)
        p_i.sml = compute_smoothing_length(p_i, particles, neighbor_list, n_neighbor,
                                          periodic, kernel);
        const real h_i = p_i.sml;

        // 2. Compute particle volume Vp
        const real Vp = compute_volume(p_i, particles, neighbor_list, n_neighbor,
                                      periodic, kernel, h_i);

        // 3. Compute baryon number density N = ν/Vp  (Eq. 37)
        p_i.N = p_i.nu / Vp;

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

        // 5. Compute gradients for MUSCL reconstruction
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
