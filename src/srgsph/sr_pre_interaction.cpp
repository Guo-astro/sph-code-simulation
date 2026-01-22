#include "srgsph/sr_pre_interaction.hpp"
#include "srgsph/sr_primitive_recovery.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "periodic.hpp"
#include "openmp.hpp"
#include "kernel/kernel_function.hpp"
#include "bhtree.hpp"
#include <iostream>


namespace sph
{
namespace srgsph
{

void PreInteraction::initialize(std::shared_ptr<SPHParameters> param)
{
    sph::PreInteraction::initialize(param);

    m_gamma = param->physics.gamma;
    m_c_speed = param->srgsph.c_speed;
    m_is_2nd_order = param->srgsph.is_2nd_order;
    m_iteration = param->iterative_sml;
    m_eta = param->srgsph.eta;
    m_c_smooth = 2.0; // Default from Python (line 51)
    m_first_calculation = true; // First calculation needs to initialize conserved from primitives
}

/**
 * Compute particle volume V_p(x_i) = [Σ_j W(x_i - x_j, h)]^(-1)
 * and grad-h correction factor Ω_i
 * Based on Python lines 99-118
 *
 * The grad-h correction accounts for variable smoothing length:
 *   Ω_i = 1 / (1 + (h_i / (D * N_i)) * dN/dh)
 *
 * Since N = ν * Σ_j W, we have dN/dh = ν * Σ_j dW/dh, so:
 *   Ω_i = 1 / (1 + h * Σ_j dW/dh / (D * Σ_j W))
 */
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
        return 1.0; // Safety fallback
    }

    // Grad-h correction factor: Ω = 1 / (1 + h * Σ dW/dh / (D * Σ W))
    // This corrects for variable smoothing length in the kernel gradient
    const real dh_term = h * sum_dW_dh / (DIM * sum_W);
    gradh_out = 1.0 / (1.0 + dh_term);

    return 1.0 / sum_W;
}

/**
 * Compute variable smoothing length using Kitajima volume-based approach
 * Based on Python lines 80-131
 * Iteratively solves: h = � * [V_p*]^(1/d) where V_p* = [�_j W(x-x_j, C_smooth*h)]^(-1)
 */
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
    real h = p_i.sml; // Initial guess from previous step

    // Iterate to find self-consistent h
    // Use many iterations and tight tolerance for uniform h in identical-state regions
    constexpr int max_iter = 50;
    for (int iter = 0; iter < max_iter; ++iter) {
        // Calculate V_p* using C_smooth * h
        real sum_W_star = 0.0;

        for (int n = 0; n < n_neighbor; ++n) {
            const int j = neighbor_list[n];
            const auto & p_j = particles[j];

            const vec_t r_ij = periodic->calc_r_ij(p_i.pos, p_j.pos);
            const real r = std::abs(r_ij);

            sum_W_star += kernel->w(r, m_c_smooth * h);
        }

        if (sum_W_star < 1e-15) {
            break; // Safety
        }

        const real Vp_star = 1.0 / sum_W_star;

        // Update h: h = � * (V_p*)^(1/d)
        const real h_new = m_eta * std::pow(Vp_star, 1.0 / DIM);

        // Check convergence - use tight tolerance to ensure uniform h
        // in regions with identical physical states
        if (std::abs(h_new - h) / h < 1e-12) {
            h = h_new;
            break;
        }

        h = h_new;
    }

    return h;
}

/**
 * Initial smoothing length calculation for first timestep
 * Based on standard SPH guess
 */
void PreInteraction::initial_smoothing(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    const int num = sim->get_particle_num();

#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        auto & p_i = particles[i];

        // If sml is already set (e.g. by make_sr_sod), don't overwrite it with a guess!
        if (p_i.sml > 1.0e-10) continue;

        // Initial guess: h ~ (ν/N)^(1/d) * η
        // Use particle's initial properties
        constexpr real A = DIM == 1 ? 2.0 :
                           DIM == 2 ? M_PI :
                                      4.0 * M_PI / 3.0;

        // For SR-SPH, mass is baryon number �, dens is lab-frame N
        p_i.sml = std::pow(m_neighbor_number * p_i.mass / (p_i.dens * A), 1.0 / DIM) * m_kernel_ratio;
    }
}

/**
 * Main calculation loop
 * Based on Python SRGSPH class methods (lines 80-283)
 */
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

    // MUSCL gradient arrays (if 2nd order enabled)
    std::vector<vec_t> * grad_p = nullptr;
    std::vector<vec_t> * grad_rho = nullptr;
    std::vector<vec_t> * grad_v[DIM] = {nullptr};

    if (m_is_2nd_order) {
        grad_p = &sim->get_vector_array("grad_pressure");
        grad_rho = &sim->get_vector_array("grad_density");
        grad_v[0] = &sim->get_vector_array("grad_velocity_0");
#if DIM >= 2
        grad_v[1] = &sim->get_vector_array("grad_velocity_1");
#endif
#if DIM == 3
        grad_v[2] = &sim->get_vector_array("grad_velocity_2");
#endif
    }

#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        
        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

        // Neighbor search
        // We need a large enough radius to cover:
        // 1. Standard kernel W(r, h) -> radius ~ 3h
        // 2. Expanded kernel W(r, C_smooth*h) for Vp* -> radius ~ 3 * C_smooth * h
        // With C_smooth=2.0, we need ~6h. Using 6.0h to be safe.
        const real search_r = p_i.sml * 6.0;
        const int n_neighbor_tmp = tree->neighbor_search(p_i, neighbor_list, particles, false);

        // Update smoothing length if iteration enabled
        // Skip on first timestep to preserve exact h from initial conditions
        // This ensures uniform initial conditions have exactly equal h values
        if (m_iteration && !m_first_calculation) {
            p_i.sml = compute_smoothing_length(p_i, particles, neighbor_list,
                                               n_neighbor_tmp, periodic, kernel, p_i.sml);
        }

        // Ghost particles: only compute h, skip density/primitive updates
        // Their other properties are set by update_ghost_particles() mirroring
        // Set gradh = 1.0 so force terms don't vanish at boundaries
        if (p_i.is_ghost) {
            p_i.gradh = 1.0;
            continue;
        }

        // Compute particle volume V_p (used for number density) and grad-h correction
        // V_p = [Σ_j W(r, h)]^(-1)
        // Ω = 1 / (1 + h * Σ dW/dh / (D * Σ W))
        real gradh_i;
        const real Vp = compute_volume(p_i, particles, neighbor_list,
                                       n_neighbor_tmp, periodic, kernel, p_i.sml, gradh_i);
        p_i.gradh = gradh_i;

        // Number density: N = ν / V_p  (line 204 in Python)
        const real N_new = p_i.nu / Vp;

        // First timestep: recompute conserved variables from primitives using correct N
        // (Initial conditions set primitives, but N was estimated, so S and e are inconsistent)
        if (m_first_calculation) {
            // Use existing primitive values (set by initialization)
            // Compute velocity magnitude including all velocity components
            // In DIM=1: v² = v_x² + v_t² (tangent velocity for 1D tests)
            // In DIM≥2: v² = v_x² + v_y² [+ v_z²] (no tangent velocity needed)
            real v2 = inner_product(p_i.vel, p_i.vel);
#if DIM == 1
            // In 1D, tangent velocity contributes to Lorentz factor
            const real vel_t = p_i.vel_t;
            v2 += vel_t * vel_t;
#endif
            const real gamma_lor = 1.0 / std::sqrt(1.0 - v2 / (m_c_speed * m_c_speed));
            p_i.gamma_lor = gamma_lor;

            // Convert updated lab-frame density N_new into rest-frame density n = N / gamma
            const real rest_density = std::max(N_new / gamma_lor, 1e-20);
            const real u_internal = p_i.pres / ((m_gamma - 1.0) * rest_density);
            const real H = 1.0 + u_internal / (m_c_speed * m_c_speed) + p_i.pres / (rest_density * m_c_speed * m_c_speed);
            p_i.enthalpy = H;

            // Compute conserved variables from primitives using CORRECT N
            // Canonical momentum: S = γHv (vector in DIM dimensions)
            p_i.N = N_new;
            p_i.S = p_i.vel * (gamma_lor * H);  // Vector: S[i] = γH * v[i]
#if DIM == 1
            // In 1D, store tangent momentum separately
            p_i.S_t = vel_t * (gamma_lor * H);
#else
            p_i.S_t = 0.0;  // Not used in 2D/3D
#endif
            const real X = m_gamma / (m_gamma - 1.0);
            p_i.e = (H * (X * gamma_lor * gamma_lor - 1.0) + 1.0) / (X * gamma_lor);

            // Update primitive storage to rest-frame values for snapshots/output
            p_i.dens = rest_density;
            p_i.ene = u_internal;

            // Sound speed
            p_i.sound = std::sqrt((m_gamma - 1.0) * (H - 1.0) / H) * m_c_speed;
        } else {
            // Normal timestep: recover primitives from conserved variables
            p_i.N = N_new;

            // IMPORTANT: In 1D SRGSPH, the tangent MOMENTUM S_t = γHv_t is conserved
            // (dS_t/dt = 0 because there's no force in the tangent direction).
            // The tangent VELOCITY v_t should be recovered from the conserved S_t.
            //
            // Use version with CONSERVED S_t (not fixed v_t) for correct physics.
            // v_t will be computed from: v_t = S_t / (γH)
            const auto prim = PrimitiveRecovery::conserved_to_primitive_with_tangent(
                p_i.S, p_i.S_t, p_i.e, p_i.N, m_gamma, m_c_speed);

            // Store primitive variables in particle
            p_i.vel = prim.vel;
            p_i.vel_t = prim.vel_t;  // v_t is recovered from conserved S_t
            p_i.dens = prim.density;         // Rest-frame density n
            p_i.pres = prim.pressure;
            p_i.sound = prim.sound_speed;
            p_i.gamma_lor = prim.gamma_lor;
            p_i.enthalpy = prim.enthalpy;

            // S_t remains constant (conserved quantity, no force in tangent direction)
            // Don't update S_t here - it's a conserved variable!
        }

        // Compute gradients for MUSCL reconstruction (if 2nd order)
        if (m_is_2nd_order) {
            vec_t grad_p_i(0.0);
            vec_t grad_rho_i(0.0);
            vec_t grad_vx_i(0.0), grad_vy_i(0.0), grad_vz_i(0.0);

            // SPH gradient: A_i H �_j V_j (A_j - A_i) W_ij
            // Based on Python lines 224-283
            for (int n = 0; n < n_neighbor_tmp; ++n) {
                const int j = neighbor_list[n];
                if (i == j) continue;

                const auto & p_j = particles[j];
                const vec_t r_ij = periodic->calc_r_ij(p_i.pos, p_j.pos);
                const real r = std::abs(r_ij);

                if (r >= p_i.sml * 2.0) continue; // Support radius

                // Gradient of kernel
                const vec_t grad_W = kernel->dw(r_ij, r, p_i.sml);

                // Volume of j (approximation: use current Vp)
                const real Vp_j = p_j.nu / p_j.N;

                // Accumulate gradients
                grad_p_i += grad_W * (Vp_j * (p_j.pres - p_i.pres));
                grad_rho_i += grad_W * (Vp_j * (p_j.dens - p_i.dens));

                grad_vx_i += grad_W * (Vp_j * (p_j.vel[0] - p_i.vel[0]));
#if DIM >= 2
                grad_vy_i += grad_W * (Vp_j * (p_j.vel[1] - p_i.vel[1]));
#endif
#if DIM == 3
                grad_vz_i += grad_W * (Vp_j * (p_j.vel[2] - p_i.vel[2]));
#endif
            }

            // Store gradients
            (*grad_p)[i] = grad_p_i;
            (*grad_rho)[i] = grad_rho_i;
            (*grad_v[0])[i] = grad_vx_i;
#if DIM >= 2
            (*grad_v[1])[i] = grad_vy_i;
#endif
#if DIM == 3
            (*grad_v[2])[i] = grad_vz_i;
#endif
        }

        p_i.neighbor = n_neighbor_tmp;
    }

    // After first calculation, switch to normal primitive recovery mode
    m_first_calculation = false;
}

} // namespace srgsph
} // namespace sph
