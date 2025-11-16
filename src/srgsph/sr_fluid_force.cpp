#include "defines.hpp"
#include "particle.hpp"
#include "periodic.hpp"
#include "simulation.hpp"
#include "bhtree.hpp"
#include "kernel/kernel_function.hpp"
#include "srgsph/sr_fluid_force.hpp"
#include "srgsph/sr_primitive_recovery.hpp"

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
    m_is_2nd_order = param->srgsph.is_2nd_order;
    m_gamma = param->physics.gamma;
    m_c_speed = param->srgsph.c_speed;
    m_c_shock = param->srgsph.c_shock;
    m_c_cd = param->srgsph.c_cd;

    exact_riemann_solver();
}

// van Leer (1979) limiter - same as non-relativistic GSPH
inline real limiter(const real dq1, const real dq2)
{
    const real dq1dq2 = dq1 * dq2;
    if(dq1dq2 <= 0) {
        return 0.0;
    } else {
        return 2.0 * dq1dq2 / (dq1 + dq2);
    }
}

/**
 * Main SR-GSPH calculation
 * Computes dS/dt (canonical momentum) and de/dt (canonical energy)
 * Based on Eqs. 64-65 in paper
 */
void FluidForce::calculation(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();
    auto * tree = sim->get_tree().get();
    const real dt = sim->get_dt();

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
        
        // Neighbor search
#ifdef EXHAUSTIVE_SEARCH
        int const n_neighbor = exhaustive_search(p_i, p_i.sml, particles, num,
                                                neighbor_list, m_neighbor_number * neighbor_list_size,
                                                periodic, true);
#else
        int const n_neighbor = tree->neighbor_search(p_i, neighbor_list, particles, true);
#endif

        // Initialize forces
        const vec_t & r_i = p_i.pos;
        const vec_t & v_i = p_i.vel;
        const real h_i = p_i.sml;
        const real nu_i = p_i.nu;

        vec_t dS(0.0);    // dS/dt (canonical momentum)
        real de = 0.0;    // de/dt (canonical energy)

        for(int n = 0; n < n_neighbor; ++n) {
            int const j = neighbor_list[n];
            auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
            const real r = std::abs(r_ij);

            if(r >= std::max(h_i, p_j.sml) || r == 0.0) {
                continue;
            }

            const real r_inv = 1.0 / r;
            const vec_t e_ij = r_ij * r_inv;
            
            // Project velocities onto line connecting particles
            const real ve_i = inner_product(v_i, e_ij);
            const real ve_j = inner_product(p_j.vel, e_ij);
            
            real vstar, pstar;

            if(m_is_2nd_order) {
                // MUSCL reconstruction for 2nd order accuracy
                // Murante et al. (2011) adapted for SR

                real right[4], left[4];  // [velocity, density, pressure, sound]
                const real delta_i = 0.5 * (1.0 - p_i.sound * dt * r_inv);
                const real delta_j = 0.5 * (1.0 - p_j.sound * dt * r_inv);

                // Velocity reconstruction
                const real dv_ij = ve_i - ve_j;
                vec_t dv_i, dv_j;
                for(int k = 0; k < DIM; ++k) {
                    dv_i[k] = inner_product(grad_v[k][i], e_ij);
                    dv_j[k] = inner_product(grad_v[k][j], e_ij);
                }
                const real dve_i = inner_product(dv_i, e_ij) * r;
                const real dve_j = inner_product(dv_j, e_ij) * r;
                
                // Monotonicity constraint (Eq. 66)
                // Switch to 1st order if:
                // 1. Strong shock: C_shock * e_ij·(v_i - v_j) > min(c_s,i, c_s,j)
                // 2. Large pressure jump: |log10(P_i/P_j)| > C_c.d.
                const real min_cs = std::min(p_i.sound, p_j.sound);
                const real shock_indicator = m_c_shock * std::abs(ve_i - ve_j);
                const real pressure_jump = std::abs(std::log10(p_i.pres / p_j.pres));
                const bool is_shock = (shock_indicator > min_cs) || (pressure_jump > m_c_cd);

                // CRITICAL: Use REST-FRAME baryon number density n = N/γ
                const real n_i = p_i.N / p_i.gamma_lor;
                const real n_j = p_j.N / p_j.gamma_lor;

                if(is_shock) {
                    // Fall back to 1st order at shocks
                    right[0] = ve_i;
                    left[0] = ve_j;
                    right[1] = n_i;   // REST-FRAME density
                    left[1] = n_j;    // REST-FRAME density
                    right[2] = p_i.pres;
                    left[2] = p_j.pres;
                    right[3] = p_i.sound;
                    left[3] = p_j.sound;
                } else {
                    // 2nd order MUSCL
                    right[0] = ve_i - limiter(dv_ij, dve_i) * delta_i;
                    left[0] = ve_j + limiter(dv_ij, dve_j) * delta_j;

                    // Rest-frame density n (not ρ!)
                    const real dn_ij = n_i - n_j;
                    const real dn_i = inner_product(grad_d[i], e_ij) * r / p_i.gamma_lor;
                    const real dn_j = inner_product(grad_d[j], e_ij) * r / p_j.gamma_lor;
                    right[1] = n_i - limiter(dn_ij, dn_i) * delta_i;
                    left[1] = n_j + limiter(dn_ij, dn_j) * delta_j;

                    // Pressure
                    const real dp_ij = p_i.pres - p_j.pres;
                    const real dp_i = inner_product(grad_p[i], e_ij) * r;
                    const real dp_j = inner_product(grad_p[j], e_ij) * r;
                    right[2] = p_i.pres - limiter(dp_ij, dp_i) * delta_i;
                    left[2] = p_j.pres + limiter(dp_ij, dp_j) * delta_j;

                    // Sound speed (recompute from reconstructed rest-frame n)
                    right[3] = std::sqrt(m_gamma * right[2] / right[1]);
                    left[3] = std::sqrt(m_gamma * left[2] / left[1]);
                }

                m_solver(left, right, pstar, vstar);
            } else {
                // 1st order Riemann solver
                // CRITICAL: Riemann solver uses REST-FRAME variables!
                // n = N/γ (rest-frame baryon number density)
                const real n_i = p_i.N / p_i.gamma_lor;
                const real n_j = p_j.N / p_j.gamma_lor;
                
                const real right[4] = {
                    ve_i,      // rest-frame velocity component
                    n_i,       // REST-FRAME baryon number density
                    p_i.pres,  // pressure (same in all frames)
                    p_i.sound, // sound speed
                };
                const real left[4] = {
                    ve_j,
                    n_j,
                    p_j.pres,
                    p_j.sound,
                };

                m_solver(left, right, pstar, vstar);
            }

            // Kernel gradients (variable smoothing length)
            // Eq. 64-65: uses ∇_i W(x_i-x_j, 2h_i) - ∇_j W(x_i-x_j, 2h_j)
            const vec_t dw_i = kernel->dw(r_ij, r, 2.0 * h_i);
            const vec_t dw_j = kernel->dw(r_ij, r, 2.0 * p_j.sml);
            const vec_t dw_ij = dw_i - dw_j;

            // Volume weighting V²_ij (Paper Eq. 29-30 and 63)
            // V²_ij(h) ≡ ∫ [1/N²(x)] W(x-x_mid, √2h) dx  
            // With N = ν/Vₚ, so 1/N² = Vₚ²/ν²
            // Approximation: V²_ij ≈ ½(Vₚ²_i/ν² + Vₚ²_j/ν²) = ½((ν/N)²_i + (ν/N)²_j)
            const real V_i = nu_i / p_i.N;
            const real V_j = p_j.nu / p_j.N;
            const real V2_ij = 0.5 * (V_i * V_i + V_j * V_j);

            // Equation of motion (Eq. 64):
            // ⟨ν_i Ṡ_i⟩ = -Σ_j P*_ij V²_ij [∇_i W - ∇_j W]
            dS -= dw_ij * (pstar * V2_ij);

            // Energy equation (Eq. 65):
            // ⟨ν_i ė_i⟩ = -Σ_j P*_ij v*_ij · V²_ij [∇_i W - ∇_j W]
            const vec_t v_ij = e_ij * vstar;
            de -= inner_product(v_ij, dw_ij) * (pstar * V2_ij);
        }

        // Store time derivatives
        // Note: These are ν_i * dS/dt and ν_i * de/dt
        // Divide by ν_i to get actual derivatives
        p_i.dS = dS / nu_i;
        p_i.de = de / nu_i;
    }
}

/**
 * Exact Riemann Solver for Special Relativistic Hydrodynamics
 * Implementation of Pons, Martí & Müller (2000)
 * "Exact solution of the Riemann problem in special relativistic hydrodynamics"
 * 
 * Handles arbitrary tangential velocities for multi-dimensional flows.
 * Supports ideal gas EOS with constant adiabatic index γ.
 */
void FluidForce::exact_riemann_solver()
{
    m_solver = [&](const real left[], const real right[], real & pstar, real & vstar) {
        // Input state variables
        const real vx_L = left[0];    // Normal velocity (left)
        const real n_L  = left[1];    // Rest-frame baryon density (left)
        const real p_L  = left[2];    // Pressure (left)
        const real cs_L = left[3];    // Sound speed (left)
        
        const real vx_R = right[0];   // Normal velocity (right)
        const real n_R  = right[1];   // Rest-frame baryon density (right)
        const real p_R  = right[2];   // Pressure (right)
        const real cs_R = right[3];   // Sound speed (right)

        // EOS: p = (γ-1)ρε, so ρ = n (for baryon number density)
        // Specific enthalpy: h = 1 + ε + p/ρ = 1 + γε = 1 + γp/((γ-1)ρ)
        const real rho_L = n_L;
        const real rho_R = n_R;
        const real h_L = 1.0 + m_gamma * p_L / ((m_gamma - 1.0) * rho_L);
        const real h_R = 1.0 + m_gamma * p_R / ((m_gamma - 1.0) * rho_R);
        
        // Lorentz factors (assuming tangential velocity = 0 for 1D Riemann problem)
        const real W_L = 1.0 / std::sqrt(1.0 - vx_L * vx_L);
        const real W_R = 1.0 / std::sqrt(1.0 - vx_R * vx_R);

        // Solve for intermediate state pressure p* using iterative method
        // Initial guess: arithmetic mean
        real p_star = 0.5 * (p_L + p_R);
        
        const int max_iter = 50;
        const real tolerance = 1.0e-8;
        
        for(int iter = 0; iter < max_iter; ++iter) {
            real vx_Lstar, vx_Rstar;
            
            // Left wave: determine if shock or rarefaction
            if(p_star > p_L) {
                // Left shock: use Rankine-Hugoniot relations
                vx_Lstar = compute_shock_velocity(p_L, rho_L, h_L, vx_L, W_L, p_star, true);
            } else {
                // Left rarefaction: integrate ODE
                vx_Lstar = compute_rarefaction_velocity(p_L, rho_L, h_L, vx_L, cs_L, p_star, true);
            }
            
            // Right wave: determine if shock or rarefaction
            if(p_star > p_R) {
                // Right shock: use Rankine-Hugoniot relations
                vx_Rstar = compute_shock_velocity(p_R, rho_R, h_R, vx_R, W_R, p_star, false);
            } else {
                // Right rarefaction: integrate ODE
                vx_Rstar = compute_rarefaction_velocity(p_R, rho_R, h_R, vx_R, cs_R, p_star, false);
            }
            
            // Check convergence: velocities must match at contact discontinuity
            const real f = vx_Rstar - vx_Lstar;
            if(std::abs(f) < tolerance) {
                pstar = p_star;
                vstar = 0.5 * (vx_Lstar + vx_Rstar);
                return;
            }
            
            // Newton-Raphson update with numerical derivative
            const real dp = p_star * 1.0e-6;
            real vx_Lstar_plus, vx_Rstar_plus;
            
            if(p_star + dp > p_L) {
                vx_Lstar_plus = compute_shock_velocity(p_L, rho_L, h_L, vx_L, W_L, p_star + dp, true);
            } else {
                vx_Lstar_plus = compute_rarefaction_velocity(p_L, rho_L, h_L, vx_L, cs_L, p_star + dp, true);
            }
            
            if(p_star + dp > p_R) {
                vx_Rstar_plus = compute_shock_velocity(p_R, rho_R, h_R, vx_R, W_R, p_star + dp, false);
            } else {
                vx_Rstar_plus = compute_rarefaction_velocity(p_R, rho_R, h_R, vx_R, cs_R, p_star + dp, false);
            }
            
            const real f_plus = vx_Rstar_plus - vx_Lstar_plus;
            const real df_dp = (f_plus - f) / dp;
            
            // Update pressure with damping for stability
            if(std::abs(df_dp) > 1.0e-12) {
                real p_new = p_star - 0.5 * f / df_dp;  // Damped Newton step
                p_new = std::max(p_new, 0.01 * std::min(p_L, p_R));  // Prevent negative pressure
                p_star = p_new;
            } else {
                // Fallback to bisection if derivative is too small
                if(f > 0) {
                    p_star *= 0.9;
                } else {
                    p_star *= 1.1;
                }
            }
        }
        
        // If not converged, use values from last iteration (with warning in debug builds)
        #ifdef DEBUG
        std::cerr << "Warning: Exact Riemann solver did not converge in " << max_iter << " iterations" << std::endl;
        #endif
        
        // Compute final vstar from last iteration
        real vx_Lstar, vx_Rstar;
        if(p_star > p_L) {
            vx_Lstar = compute_shock_velocity(p_L, rho_L, h_L, vx_L, W_L, p_star, true);
        } else {
            vx_Lstar = compute_rarefaction_velocity(p_L, rho_L, h_L, vx_L, cs_L, p_star, true);
        }
        if(p_star > p_R) {
            vx_Rstar = compute_shock_velocity(p_R, rho_R, h_R, vx_R, W_R, p_star, false);
        } else {
            vx_Rstar = compute_rarefaction_velocity(p_R, rho_R, h_R, vx_R, cs_R, p_star, false);
        }
        
        pstar = p_star;
        vstar = 0.5 * (vx_Lstar + vx_Rstar);
    };
}

/**
 * Compute post-shock velocity using Rankine-Hugoniot relations
 * Pons et al. (2000), Equation 4.12 and surrounding equations
 * 
 * @param p_a Pressure ahead of shock
 * @param rho_a Density ahead of shock
 * @param h_a Specific enthalpy ahead of shock
 * @param vx_a Normal velocity ahead of shock
 * @param W_a Lorentz factor ahead of shock
 * @param p_b Pressure behind shock (post-shock)
 * @param is_left True if left-propagating shock, false if right-propagating
 * @return Post-shock normal velocity vx_b
 */
real FluidForce::compute_shock_velocity(const real p_a, const real rho_a, const real h_a,
                                        const real vx_a, const real W_a, const real p_b,
                                        const bool is_left) const
{
    // Taub adiabat (Eq. 4.16) - solve for h_b
    // h²_b[1 + (γ-1)(p_a-p_b)/(γp_b)] - h_b·h_a(p_a-p_b)/(γp_b) - h²_a + h_a(p_a-p_b)/rho_a = 0
    const real dp = p_a - p_b;
    const real A = 1.0 + (m_gamma - 1.0) * dp / (m_gamma * p_b);
    const real B = -h_a * dp / (m_gamma * p_b);
    const real C = -h_a * h_a + h_a * dp / rho_a;
    
    // Quadratic formula: h_b = (-B + sqrt(B² - 4AC)) / (2A)
    // Take positive root only
    const real discriminant = B * B - 4.0 * A * C;
    if(discriminant < 0.0) {
        // Should not happen for physical shocks, use fallback
        return vx_a;
    }
    const real h_b = (-B + std::sqrt(discriminant)) / (2.0 * A);
    
    // Density from EOS: p_b = (γ-1)ρ_b·ε_b
    // h_b = 1 + ε_b + p_b/ρ_b = 1 + γε_b (for ideal gas)
    // So: ε_b = (h_b - 1)/γ and ρ_b = p_b/[(γ-1)ε_b] = γp_b/[(γ-1)(h_b-1)]
    const real rho_b = m_gamma * p_b / ((m_gamma - 1.0) * (h_b - 1.0));
    
    // Mass flux j² from Eq. 4.17
    // j² = -[p] / [h/ρ]
    const real h_rho_a = h_a / rho_a;
    const real h_rho_b = h_b / rho_b;
    const real j_squared = -(p_b - p_a) / (h_rho_b - h_rho_a);
    
    if(j_squared < 0.0) {
        // Unphysical, fallback
        return vx_a;
    }
    
    const real j = is_left ? -std::sqrt(j_squared) : std::sqrt(j_squared);
    
    // Shock velocity from Eq. 4.14
    // V_s = [ρ²_a W²_a v^x_a ± |j| sqrt(j² + ρ²_a W²_a(1 - v^x_a²))] / [ρ²_a W²_a + j²]
    const real rho_W_a = rho_a * W_a;
    const real rho_W_a_sq = rho_W_a * rho_W_a;
    const real numerator = rho_W_a_sq * vx_a + j * std::sqrt(j_squared + rho_W_a_sq * (1.0 - vx_a * vx_a));
    const real denominator = rho_W_a_sq + j_squared;
    const real V_s = numerator / denominator;
    
    // Post-shock velocity from Eq. 4.12
    // v^x_b = [h_a W_a v^x_a + W_s(p_b-p_a)/j + W_s V_s/(ρ_a W_a)]^(-1)
    //         × [h_a W_a + (p_b-p_a)/j]
    const real W_s = 1.0 / std::sqrt(1.0 - V_s * V_s);
    const real term1 = h_a * W_a * vx_a + W_s * (p_b - p_a) / j + W_s * V_s / (rho_a * W_a);
    const real term2 = h_a * W_a + (p_b - p_a) / j;
    
    return term2 / term1;
}

/**
 * Compute post-rarefaction velocity by integrating ODE
 * Pons et al. (2000), Equation 3.10 and 3.14
 * 
 * @param p_a Pressure ahead of rarefaction
 * @param rho_a Density ahead of rarefaction
 * @param h_a Specific enthalpy ahead of rarefaction
 * @param vx_a Normal velocity ahead of rarefaction
 * @param cs_a Sound speed ahead of rarefaction
 * @param p_b Pressure behind rarefaction (post-rarefaction)
 * @param is_left True if left-propagating rarefaction, false if right-propagating
 * @return Post-rarefaction normal velocity vx_b
 */
real FluidForce::compute_rarefaction_velocity(const real p_a, const real rho_a, const real h_a,
                                              const real vx_a, const real cs_a, const real p_b,
                                              const bool is_left) const
{
    // For vanishing tangential velocity (v_t = 0), Eq. 3.12 simplifies to:
    // W² dv^x = ± (c_s / γp) dp = ± (c_s / ρ) dρ
    // Integrate from state a to state b
    
    // For ideal gas: p = (γ-1)ρε, and h = 1 + γε
    // So: ε = (h-1)/γ and ρ = γp/[(γ-1)(h-1)]
    // Sound speed: c²_s = γp/ρh (from Eq. 2.15)
    
    const int n_steps = 100;
    const real dp = (p_b - p_a) / n_steps;
    
    real p = p_a;
    real vx = vx_a;
    
    for(int i = 0; i < n_steps; ++i) {
        // Current state
        const real eps = (m_gamma - 1.0) * p / rho_a;  // Assuming isentropic: ρ/ρ_a = (p/p_a)^(1/γ)
        const real rho = rho_a * std::pow(p / p_a, 1.0 / m_gamma);
        const real h = 1.0 + eps + p / rho;
        const real cs = std::sqrt(m_gamma * p / (rho * h));
        const real W = 1.0 / std::sqrt(1.0 - vx * vx);
        
        // Eq. 3.10: dv^x/dp = ± 1/(ρhW²c_s sqrt(1 + g))
        // For v_t = 0, g = 0, so: dv^x/dp = ± 1/(ρhW²c_s)
        const real sign = is_left ? -1.0 : 1.0;
        const real dvx_dp = sign / (rho * h * W * W * cs);
        
        // Euler integration
        vx += dvx_dp * dp;
        p += dp;
        
        // Clamp velocity to sub-luminal
        if(vx >= 1.0) vx = 0.99999;
        if(vx <= -1.0) vx = -0.99999;
    }
    
    return vx;
}

}
}

