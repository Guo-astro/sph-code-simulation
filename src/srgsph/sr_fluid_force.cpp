#include "srgsph/sr_fluid_force.hpp"
#include "srgsph/sr_exact_riemann.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "periodic.hpp"
#include "openmp.hpp"
#include "kernel/kernel_function.hpp"
#include "bhtree.hpp"
#include "logger.hpp"
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
    m_riemann_solver = param->srgsph.riemann_solver;
}

/**
 * Sophisticated vacuum state fallback for relativistic Riemann problem
 * 
 * When the exact solver fails (typically in near-vacuum states with high tangent velocity),
 * this function provides a physically-motivated approximation based on:
 * 
 * 1. Acoustic (linearized) Riemann solution for pressure
 * 2. Characteristic-based upwinding for velocity  
 * 3. Proper K = h·W·v_t invariant preservation (Pons et al. 2000)
 * 
 * The key insight is that in near-vacuum regions, the waves are weak
 * and we can use a linearized approximation that respects causality.
 */
static void vacuum_state_fallback(
    const riemann::RiemannState& left,
    const riemann::RiemannState& right,
    real v_t_L, real v_t_R,
    real gamma_c, real c,
    real& P_star, real& v_x_star, real& v_t_star)
{
    // Compute total velocities and Lorentz factors
    const real v2_L = left.v * left.v + v_t_L * v_t_L;
    const real v2_R = right.v * right.v + v_t_R * v_t_R;
    
    // Clamp to subluminal (safety)
    const real v2_L_safe = std::min(v2_L, 0.9999);
    const real v2_R_safe = std::min(v2_R, 0.9999);
    
    const real W_L = 1.0 / std::sqrt(1.0 - v2_L_safe);
    const real W_R = 1.0 / std::sqrt(1.0 - v2_R_safe);
    
    // Enthalpy: h = 1 + γ/(γ-1) * P/n
    const real h_L = 1.0 + gamma_c * left.P / ((gamma_c - 1.0) * std::max(left.n, 1e-10));
    const real h_R = 1.0 + gamma_c * right.P / ((gamma_c - 1.0) * std::max(right.n, 1e-10));
    
    // Sound speeds (ensure positive)
    const real cs_L = std::max(left.c_s, 1e-6);
    const real cs_R = std::max(right.c_s, 1e-6);
    
    // Compute K invariant (Pons et al. Eq. 22, 25): K = h * W * v_t
    const real K_L = h_L * W_L * v_t_L;
    const real K_R = h_R * W_R * v_t_R;
    
    // ========================================================================
    // ACOUSTIC APPROXIMATION for P* (linearized Riemann solution)
    // ========================================================================
    // For weak waves: P* ≈ (Z_L*P_R + Z_R*P_L + Z_L*Z_R*(v_L - v_R)) / (Z_L + Z_R)
    // where Z = ρ*h*W²*c_s is the relativistic acoustic impedance
    
    const real Z_L = left.n * h_L * W_L * W_L * cs_L;
    const real Z_R = right.n * h_R * W_R * W_R * cs_R;
    const real Z_sum = Z_L + Z_R;
    
    if (Z_sum > 1e-20) {
        // Acoustic approximation
        P_star = (Z_L * right.P + Z_R * left.P + Z_L * Z_R * (left.v - right.v)) / Z_sum;
    } else {
        // Fallback to simple average if impedances vanish
        P_star = 0.5 * (left.P + right.P);
    }
    
    // Floor pressure to positive value
    P_star = std::max(P_star, 1e-10);
    
    // ========================================================================
    // CHARACTERISTIC-BASED UPWINDING for v_x*
    // ========================================================================
    // Use wave speed to determine upwind direction
    // Contact wave moves at v_x*, approximate with weighted average
    
    if (Z_sum > 1e-20) {
        // Velocity from acoustic approximation
        v_x_star = (Z_L * left.v + Z_R * right.v + (left.P - right.P)) / Z_sum;
    } else {
        v_x_star = 0.5 * (left.v + right.v);
    }
    
    // ========================================================================
    // K-INVARIANT BASED v_t* (Pons et al. physics)
    // ========================================================================
    // The invariant K = h * W * v_t is conserved across waves
    // We need to compute new h* and W* to find v_t*
    
    // Estimate star region density using isentropic relation: n* ∝ P*^(1/γ)
    // This is approximate but respects thermodynamics
    const real n_star_L = left.n * std::pow(P_star / std::max(left.P, 1e-10), 1.0 / gamma_c);
    const real n_star_R = right.n * std::pow(P_star / std::max(right.P, 1e-10), 1.0 / gamma_c);
    const real n_star = 0.5 * (n_star_L + n_star_R);
    
    // Star region enthalpy
    const real h_star = 1.0 + gamma_c * P_star / ((gamma_c - 1.0) * std::max(n_star, 1e-10));
    
    // Star region Lorentz factor from v_x* (need to account for v_t* too)
    // This requires solving: W* = 1/√(1 - v_x*² - v_t*²)
    // and: K* = h* * W* * v_t*
    // 
    // For upwind state, use the K from the characteristic direction:
    // If v_x* > 0 (flow to right), use K_L; otherwise use K_R
    const real K_star = (v_x_star > 0) ? K_L : K_R;
    
    // Now solve for v_t*:
    // K_star = h_star * W_star * v_t*
    // W_star = 1/√(1 - v_x*² - v_t*²)
    //
    // Let u = v_t*², then:
    // K_star² * (1 - v_x*² - u) = h_star² * u
    // K_star² - K_star² * v_x*² - K_star² * u = h_star² * u
    // K_star² * (1 - v_x*²) = u * (K_star² + h_star²)
    // u = K_star² * (1 - v_x*²) / (K_star² + h_star²)
    
    const real v_x_star2 = v_x_star * v_x_star;
    const real v_x_max2 = 0.9999;  // Leave room for v_t
    
    if (v_x_star2 < v_x_max2 && std::abs(K_star) > 1e-10) {
        const real K2 = K_star * K_star;
        const real h2 = h_star * h_star;
        const real denominator = K2 + h2;
        
        if (denominator > 1e-10) {
            real v_t_star2 = K2 * (1.0 - v_x_star2) / denominator;
            
            // Safety: ensure total velocity is subluminal
            const real v_t_max2 = std::max(0.0, 0.9999 - v_x_star2);
            v_t_star2 = std::min(v_t_star2, v_t_max2);
            
            v_t_star = std::copysign(std::sqrt(v_t_star2), K_star);
        } else {
            // K very small, v_t essentially 0
            v_t_star = 0.0;
        }
    } else if (std::abs(K_star) <= 1e-10) {
        // K is zero, so v_t* = 0
        v_t_star = 0.0;
    } else {
        // v_x* is already near the speed limit, clamp v_t* to what's allowed
        const real v_t_max2 = std::max(0.0, 0.9999 - v_x_star2);
        v_t_star = std::copysign(std::sqrt(v_t_max2), K_star);
    }
    
    // Final safety: ensure v_x* respects the total velocity constraint
    const real v_t_star2 = v_t_star * v_t_star;
    const real v_x_allowed_max = std::sqrt(std::max(0.0, 0.9999 - v_t_star2));
    if (std::abs(v_x_star) > v_x_allowed_max) {
        v_x_star = std::copysign(v_x_allowed_max * 0.99, v_x_star);
    }
}

void FluidForce::solve_interface_state(
    const real left_state[5],
    const real right_state[5],
    real & P_star,
    real & v_x_star,
    real & v_t_star) const
{
    // Build RiemannState structures from input arrays
    // Array format: [v^x, n, P, c_s, v^t]
    riemann::RiemannState left{left_state[0], left_state[1], left_state[2], left_state[3]};
    riemann::RiemannState right{right_state[0], right_state[1], right_state[2], right_state[3]};
    const real v_t_L = left_state[4];
    const real v_t_R = right_state[4];

    // Use exact Riemann solver (Kitajima/Pons et al.)
    const bool converged = riemann::exact_riemann_solver(
        left, right, v_t_L, v_t_R, m_gamma, m_c_speed,
        P_star, v_x_star, v_t_star, 100, 1e-10);

    // If exact solver fails, use sophisticated vacuum state fallback
    if (!converged) {
        // Log detailed failure info (first 10 only)
        static int failure_count = 0;
        if (failure_count < 10) {
            const real v_mag_L = std::sqrt(left.v * left.v + v_t_L * v_t_L);
            const real v_mag_R = std::sqrt(right.v * right.v + v_t_R * v_t_R);
            WRITE_LOG << "[RIEMANN FAILURE #" << (failure_count + 1) << "] Using vacuum fallback."
                      << " L: v_x=" << left.v << " v_t=" << v_t_L << " |v|=" << v_mag_L
                      << " n=" << left.n << " P=" << left.P << " c_s=" << left.c_s
                      << " | R: v_x=" << right.v << " v_t=" << v_t_R << " |v|=" << v_mag_R
                      << " n=" << right.n << " P=" << right.P << " c_s=" << right.c_s;
            ++failure_count;
        }
        
        // Use physics-based vacuum state fallback
        vacuum_state_fallback(left, right, v_t_L, v_t_R, m_gamma, m_c_speed,
                             P_star, v_x_star, v_t_star);
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
        particles[i].dS_t = 0.0;  // Reset tangent momentum derivative for 1D tests
    }

    const real sqrt_two = std::sqrt(2.0);

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        
        // Skip ghost particles - they don't evolve and shouldn't accumulate forces
        // (their properties are set by update_ghost_particles via mirroring)
        if (p_i.is_ghost) continue;
        
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

            if (r > 3.0 * std::sqrt(2.0) * std::max(p_i.sml, p_j.sml)) {
                continue;
            }

            const vec_t n_ij = r_ij * (-1.0 / r);

            const vec_t midpoint = (p_i.pos + p_j.pos) * 0.5;
            const vec_t dr_i = midpoint - p_i.pos;
            const vec_t dr_j = midpoint - p_j.pos;

            bool use_second_order = gradients_available;
            
            // Don't use 2nd order reconstruction for ghost interactions
            // because ghosts don't have computed gradients
            if (p_j.is_ghost) {
                use_second_order = false;
            }
            
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
                
                // Check velocity magnitudes to ensure subluminal (|v|² < c² = 1)
                // Must include tangent velocity v_t for 1D tests with transverse motion
                const real vt_L_sq = p_i.vel_t * p_i.vel_t;
                const real vt_R_sq = p_j.vel_t * p_j.vel_t;
#if DIM == 1
                const real v2_L_trial = vx_L_trial * vx_L_trial + vt_L_sq;
                const real v2_R_trial = vx_R_trial * vx_R_trial + vt_R_sq;
#elif DIM == 2
                const real v2_L_trial = vx_L_trial * vx_L_trial + vy_L_trial * vy_L_trial + vt_L_sq;
                const real v2_R_trial = vx_R_trial * vx_R_trial + vy_R_trial * vy_R_trial + vt_R_sq;
#else
                const real v2_L_trial = vx_L_trial * vx_L_trial + vy_L_trial * vy_L_trial + vz_L_trial * vz_L_trial;
                const real v2_R_trial = vx_R_trial * vx_R_trial + vy_R_trial * vy_R_trial + vz_R_trial * vz_R_trial;
#endif
                const bool velocity_valid = (v2_L_trial < 1.0 && v2_R_trial < 1.0);

                if (valid_left && valid_right && velocity_valid) {
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
            // Get tangent velocities
            real v_t_L = p_i.vel_t;
            real v_t_R = p_j.vel_t;
            
            // Safety: clamp normal velocities to subluminal with tangent velocity
            // v_x^2 + v_t^2 < 1  =>  |v_x| < sqrt(1 - v_t^2)
            const real max_vx_L = std::sqrt(std::max(0.9999 - v_t_L * v_t_L, 0.01));
            const real max_vx_R = std::sqrt(std::max(0.9999 - v_t_R * v_t_R, 0.01));
            
            if (std::abs(vx_L) > max_vx_L) {
                static int clamp_L_count = 0;
                if (clamp_L_count < 5) {
                    WRITE_LOG << "[VEL CLAMP L] particle " << i << " vx_L=" << vx_L 
                              << " clamped to " << std::copysign(max_vx_L, vx_L)
                              << " (v_t=" << v_t_L << ", max_vx=" << max_vx_L << ")";
                    ++clamp_L_count;
                }
                vx_L = std::copysign(max_vx_L, vx_L);
            }
            if (std::abs(vx_R) > max_vx_R) {
                static int clamp_R_count = 0;
                if (clamp_R_count < 5) {
                    WRITE_LOG << "[VEL CLAMP R] particle " << j << " vx_R=" << vx_R 
                              << " clamped to " << std::copysign(max_vx_R, vx_R)
                              << " (v_t=" << v_t_R << ", max_vx=" << max_vx_R << ")";
                    ++clamp_R_count;
                }
                vx_R = std::copysign(max_vx_R, vx_R);
            }
            
            const vec_t v_L(vx_L);
            const vec_t v_R(vx_R);
#elif DIM == 2
            const vec_t v_L(vx_L, vy_L);
            const vec_t v_R(vx_R, vy_R);
#else
            const vec_t v_L(vx_L, vy_L, vz_L);
            const vec_t v_R(vx_R, vy_R, vz_R);
#endif

            // Decompose velocities into normal (v^x) and tangential (v^t) components
            // For DIM=1 with tangent velocity tests, use the separate vel_t field
            real v_x_L = inner_product(v_L, n_ij);
            real v_x_R = inner_product(v_R, n_ij);
            
#if DIM == 1
            // In 1D, tangent velocity is stored in separate vel_t field
            // This is used for Section 3.1.5 tangent velocity tests
            // (v_t_L and v_t_R already set above)
            vec_t v_t_vec_L(0.0);  // Not used in 1D for reconstruction
            vec_t v_t_vec_R(0.0);
#else
            // In 2D/3D, decompose velocity vector into normal and tangent components
            const vec_t v_t_vec_L = v_L - n_ij * v_x_L;
            const vec_t v_t_vec_R = v_R - n_ij * v_x_R;
            real v_t_L = std::abs(v_t_vec_L);
            real v_t_R = std::abs(v_t_vec_R);
#endif

            // Riemann solution variables
            real P_star = 0.0;
            real v_x_star = 0.0;
            real v_t_star = 0.0;
            vec_t v_star_vec(0.0);

            // For ghost-real interactions, use wall boundary condition directly
            // instead of solving Riemann problem which would create spurious pressure
            // Ghost particles represent a stationary reflecting wall
            if (p_j.is_ghost) {
                // Wall boundary condition:
                // - P* = P_L (wall reflects with same pressure, no shock/rarefaction)
                // - v* = 0 (no flow through wall)
                // - v_t* = v_t_L (tangential velocity preserved)
                P_star = p_L;
                v_x_star = 0.0;
                v_t_star = v_t_L;
                
                // Reconstruct star velocity (only tangential component)
#if DIM >= 2
                if (v_t_star > 1e-10 && v_t_L > 1e-10) {
                    v_star_vec = v_t_vec_L * (v_t_star / v_t_L);
                }
#else
                // In 1D, v_star_vec only has normal component (tangent handled separately)
                v_star_vec = vec_t(0.0);
#endif
            } else {
                // Normal real-real interaction: solve Riemann problem
                // State arrays: [v^x, n, P, c_s, v^t]
                real left_state[5] = {v_x_L, rho_L, p_L, p_i.sound, v_t_L};
                real right_state[5] = {v_x_R, rho_R, p_R, p_j.sound, v_t_R};

                solve_interface_state(left_state, right_state, P_star, v_x_star, v_t_star);
                
                // DEBUG: Check Riemann result for discontinuity pairs (large pressure ratio)
                static int riemann_debug = 0;
                static int discontinuity_debug = 0;
                const real p_ratio = (p_L > p_R) ? (p_L / std::max(p_R, 1e-10)) : (p_R / std::max(p_L, 1e-10));
                if (discontinuity_debug < 20 && p_ratio > 1000.0) {
                    WRITE_LOG << "[DISCONTINUITY #" << discontinuity_debug << "] P*=" << P_star
                              << " v*=" << v_x_star << " v_t*=" << v_t_star
                              << " from p_L=" << p_L << " p_R=" << p_R
                              << " v_L=" << v_x_L << " v_R=" << v_x_R
                              << " v_t_L=" << v_t_L << " v_t_R=" << v_t_R
                              << " (expected: P*≈0.89, v*≈0.32)";
                    ++discontinuity_debug;
                }
                if (riemann_debug < 10 && (P_star > 10.0 || P_star < 0.0)) {
                    WRITE_LOG << "[RIEMANN RESULT #" << riemann_debug << "] Suspicious P*=" << P_star
                              << " from p_L=" << p_L << " p_R=" << p_R
                              << " rho_L=" << rho_L << " rho_R=" << rho_R
                              << " cs_L=" << p_i.sound << " cs_R=" << p_j.sound;
                    ++riemann_debug;
                }

                // Reconstruct full star velocity vector
                v_star_vec = n_ij * v_x_star;
#if DIM >= 2
                if (v_t_star > 1e-10) {
                    if (v_x_star > 0.0 && v_t_L > 1e-10) {
                        v_star_vec = v_star_vec + v_t_vec_L * (v_t_star / v_t_L);
                    } else if (v_x_star <= 0.0 && v_t_R > 1e-10) {
                        v_star_vec = v_star_vec + v_t_vec_R * (v_t_star / v_t_R);
                    }
                }
#endif
                // In 1D, tangent momentum is handled separately via dS_t
            }

            const real V_i = p_i.nu / p_i.N;
            const real V_j = p_j.nu / p_j.N;
            const real Vij2 = 0.5 * (V_i * V_i + V_j * V_j);

            // Per Kitajima Eq. 58-59: use h_i for grad_W_i and h_j for grad_W_j
            // This anti-symmetric form guarantees momentum/energy conservation to machine precision
            // grad_W_i = ∇_i W(x_i - x_j, √2 h_i)
            // grad_W_j = ∇_j W(x_i - x_j, √2 h_j)
            const vec_t grad_W_i = kernel->dw(r_ij, r, sqrt_two * p_i.sml);
            const vec_t grad_W_j = kernel->dw(r_ij, r, sqrt_two * p_j.sml);
            const vec_t term_grad = grad_W_i - (grad_W_j * (-1.0));

            const vec_t force = term_grad * (-P_star * Vij2);
            const real power = inner_product(v_star_vec, force);
            
            // DEBUG: Check for extreme forces
            static int force_debug_count = 0;
            if (force_debug_count < 20 && std::abs(force[0]) > 100.0) {
                WRITE_LOG << "[FORCE DEBUG #" << force_debug_count << "] Large force! "
                          << "F=" << force[0]
                          << ", P*=" << P_star 
                          << ", V_i=" << V_i << ", V_j=" << V_j
                          << ", Vij2=" << Vij2
                          << ", grad_W=" << term_grad[0]
                          << ", r=" << r
                          << ", sml_i=" << p_i.sml << ", sml_j=" << p_j.sml;
                ++force_debug_count;
            }

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
        // Skip ghost particles in normalization too
        if (particles[i].is_ghost) continue;
        FluidForce::normalize_sr_derivatives(particles[i]);
    }
}

} // namespace srgsph
} // namespace sph
