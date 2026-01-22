#include "relaxation/lane_emden_relaxation.hpp"
#include "exception.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sph
{

LaneEmdenRelaxation::LaneEmdenRelaxation()
    : m_initialized(false)
{
}

void LaneEmdenRelaxation::initialize(const LaneEmdenRelaxationParams& params)
{
    m_params = params;
#if DIM == 2
    m_data.load_from_file("data/lane_emden/n1.5_2d.dat");
#else
    m_data.load_from_file("data/lane_emden/n1.5_3d.dat");
#endif
    m_initialized = true;
    
    std::cout << "LaneEmdenRelaxation: Initialized with:" << std::endl;
    std::cout << "  R = " << params.R << ", M = " << params.M_total << std::endl;
    std::cout << "  ρ_c = " << params.rho_center << ", K = " << params.K << std::endl;
    std::cout << "  α = " << params.alpha_scaling << ", G = " << params.G << std::endl;
}

vec_t LaneEmdenRelaxation::compute_relaxation_force(const SPHParticle& p) const
{
    if(!m_initialized) {
        THROW_ERROR("LaneEmdenRelaxation not initialized");
    }
    
    vec_t force(0.0);
    
    // Spherical radius
    const real r = std::abs(p.pos);
    
    if(r < 1e-12) {
        return force;  // No force at center
    }
    
    // Dimensionless coordinate
    const real xi = r / m_params.alpha_scaling;
    const real xi_1 = m_data.get_xi_1();

    if(xi >= xi_1) {
        // Outside sphere: no relaxation force, particle will be removed
        return force;
    }
    
    // Get Lane-Emden solution
    const real theta = m_data.get_theta(xi);
    const real dtheta = m_data.dtheta_dxi(xi);
    
    if(theta <= 0.0) {
        return force;
    }
    
    // Compute analytical acceleration from Lane-Emden solution
    // For polytrope n=1.5: a_r = -GM(r)/r² where M(r) from Lane-Emden
    // Or equivalently: a_r = -(1/ρ) dP/dr
    
    // Using pressure gradient form:
    // P = K ρ^γ, so dP/dr = K γ ρ^(γ-1) dρ/dr
    // ρ = ρ_c θ^n, so dρ/dr = ρ_c n θ^(n-1) dθ/dr
    // dθ/dr = (1/α) dθ/dξ
    // Therefore: a_r = -(1/ρ)(dP/dr) = -(K γ n / α) θ^(γ-1-n) dθ/dξ
    //                = -(K γ n / α) θ^(γ-n-1) dθ/dξ
    // For n=1.5, γ=5/3: γ-n-1 = 5/3 - 1.5 - 1 = -1.833... 
    // Actually: a_r = -(K γ n / (α ρ)) dθ/dξ where ρ = ρ_c θ^n
    //              = -(K γ n / (α ρ_c)) θ^(-n) dθ/dξ
    //              = -(K γ n / (α ρ_c)) θ^(-1.5) dθ/dξ
    
    // Wait, let me recalculate: a = -(1/ρ) dP/dr
    // P = K ρ^γ, dP/dr = K γ ρ^(γ-1) dρ/dr
    // So a = -(K γ / ρ) ρ^(γ-1) dρ/dr = -K γ ρ^(γ-2) dρ/dr
    // ρ = ρ_c θ^n, dρ/dr = ρ_c n θ^(n-1) (1/α) dθ/dξ
    // a = -K γ (ρ_c θ^n)^(γ-2) * ρ_c n θ^(n-1) (1/α) dθ/dξ
    //   = -K γ ρ_c^(γ-1) θ^(n(γ-2)) n θ^(n-1) (1/α) dθ/dξ
    //   = -K γ ρ_c^(γ-1) n / α * θ^(n(γ-2)+n-1) dθ/dξ
    //   = -K γ ρ_c^(γ-1) n / α * θ^(nγ-2n+n-1) dθ/dξ
    //   = -K γ ρ_c^(γ-1) n / α * θ^(nγ-n-1) dθ/dξ
    // For n=1.5, γ=5/3: nγ-n-1 = 1.5*5/3 - 1.5 - 1 = 2.5 - 2.5 = 0
    // So: a = -K γ ρ_c^(γ-1) n / α * dθ/dξ
    
    // Correct formula derived from hydrostatic equilibrium:
    // a_r = -(1/ρ) dP/dr = -K γ ρ^(γ-2) dρ/dr
    // With ρ = ρ_c θ^n and dρ/dr = ρ_c n θ^(n-1) (1/α) dθ/dξ:
    // a_r = -K γ ρ_c^(γ-1) n / α * θ^(nγ - n - 1) dθ/dξ
    // For n=1.5, γ=5/3: nγ - n - 1 = 1.5*5/3 - 1.5 - 1 = 0
    // Therefore: a_r = -K γ ρ_c^(γ-1) n / α * dθ/dξ
    
    const real n = 1.5;
    const real rho_c_pow_gamma_minus_1 = std::pow(m_params.rho_center, m_params.gamma - 1.0);
    const real prefactor = m_params.K * m_params.gamma * rho_c_pow_gamma_minus_1 * n / m_params.alpha_scaling;
    
    // Equilibrium pressure gradient acceleration: a_r = -K γ ρ_c^(γ-1) n / α * dθ/dξ
    // For interior (dθ/dξ < 0): a_r > 0 (outward pressure support balancing gravity)
    const real a_r = -prefactor * dtheta;
    
    // Convert to Cartesian vector
    // We return this as the "analytical pressure gradient force" to be subtracted
    // from total SPH acceleration: a_net = a_grav + a_pressure_SPH - a_pressure_analytical
    const real r_inv = 1.0 / r;
    force[0] = a_r * p.pos[0] * r_inv;
    force[1] = a_r * p.pos[1] * r_inv;
#if DIM == 3
    force[2] = a_r * p.pos[2] * r_inv;
#endif
    
    return force;
}

void LaneEmdenRelaxation::apply_relaxation(std::shared_ptr<Simulation> sim, real damping_factor)
{
    if(!m_initialized) {
        THROW_ERROR("LaneEmdenRelaxation not initialized");
    }
    
    auto& particles = sim->get_particles();
    const int num_p = sim->get_particle_num();
    
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(int i = 0; i < num_p; ++i) {
        // Compute relaxation acceleration (analytical pressure gradient from Lane-Emden)
        vec_t relax_acc = compute_relaxation_force(particles[i]);
        
        // SUBTRACT analytical pressure gradient from SPH acceleration
        // Net acceleration = SPH forces - Lane-Emden equilibrium pressure gradient
        // This drives particles toward equilibrium positions
        particles[i].acc[0] -= relax_acc[0];
        particles[i].acc[1] -= relax_acc[1];
#if DIM == 3
        particles[i].acc[2] -= relax_acc[2];
#endif
        
        // NOTE: Velocities are zeroed in the solver loop, not here
    }
}

int LaneEmdenRelaxation::remove_escaping_particles(std::shared_ptr<Simulation> sim, real tolerance_factor)
{
    if(!m_initialized) {
        THROW_ERROR("LaneEmdenRelaxation not initialized");
    }

    auto& particles = sim->get_particles();
    const real xi_1 = m_data.get_xi_1();
    const real R_cloud = xi_1 * m_params.alpha_scaling;
    const real R_max = R_cloud * tolerance_factor;

    // Use erase-remove idiom
    auto new_end = std::remove_if(particles.begin(), particles.end(),
        [R_max](const SPHParticle& p) {
            const real r = std::abs(p.pos);
            return r > R_max;
        });

    const int removed = std::distance(new_end, particles.end());
    particles.erase(new_end, particles.end());
    sim->set_particle_num(static_cast<int>(particles.size()));

    if(removed > 0) {
        std::cout << "LaneEmdenRelaxation: Removed " << removed << " particles (r > " << R_max << ")" << std::endl;
    }
    return removed;
}

} // namespace sph
