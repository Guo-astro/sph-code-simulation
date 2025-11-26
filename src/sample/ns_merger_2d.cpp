// 2D Neutron Star Merger
// Simplified 2D model of binary neutron star collision
// Based on: Shibata & Hotokezaka (2019), Baiotti & Rezzolla (2017)
//
// Two neutron stars modeled as n=1 polytropes (Γ=2) collide head-on
// in the xy-plane. Uses SR-GSPH for relativistic hydrodynamics.
//
// Physical scales (code units with c=1):
//   Length: 10 km (typical NS scale)
//   Time: L/c ≈ 0.033 ms
//   Density: 10^14 g/cm³ (nuclear density scale)
//
// Initial conditions:
//   - Two NS with radius ~1.2 code units (12 km)
//   - Separation ~6 code units (60 km)
//   - Collision velocity ~0.1-0.2 c

#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"
#include <cmath>
#include <vector>
#include <iostream>

namespace sph {

// Lane-Emden density profile for n=1 polytrope: rho(r) = rho_c * sin(pi*r/R) / (pi*r/R)
static real lane_emden_n1_density(real r, real R, real rho_c) {
    if (r < 1e-10) return rho_c;  // Central value
    real xi = M_PI * r / R;
    if (xi > M_PI) return 0.0;  // Outside the star
    return rho_c * std::sin(xi) / xi;
}

// Pressure from polytropic EOS: P = K * rho^Gamma
static real polytropic_pressure(real rho, real K, real Gamma) {
    return K * std::pow(rho, Gamma);
}

void Solver::make_ns_merger_2d()
{
#if DIM != 2
    THROW_ERROR("NS Merger 2D requires DIM == 2");
#else
    std::cout << "=== Initializing 2D Neutron Star Merger ===" << std::endl;
    
    // Get parameters from config
    const real gamma = m_param->physics.gamma;  // Should be 2.0 for Γ=2 polytrope
    
    // Default parameters (can be overridden by config)
    real R_star = 1.2;      // NS radius in code units (12 km)
    real rho_c = 2.8;       // Central density (2.8 × 10^14 g/cm³ in physical)
    real separation = 6.0;  // Initial separation (60 km)
    real v_collision = 0.15; // Collision velocity (0.15c)
    int n_radial = 30;      // Particles per radius
    
    // Read from config if present
    try {
        if (m_sample_parameters.count("R_star"))
            R_star = boost::any_cast<real>(m_sample_parameters["R_star"]);
        if (m_sample_parameters.count("rho_c"))
            rho_c = boost::any_cast<real>(m_sample_parameters["rho_c"]);
        if (m_sample_parameters.count("separation"))
            separation = boost::any_cast<real>(m_sample_parameters["separation"]);
        if (m_sample_parameters.count("v_collision"))
            v_collision = boost::any_cast<real>(m_sample_parameters["v_collision"]);
        if (m_sample_parameters.count("n_radial"))
            n_radial = boost::any_cast<int>(m_sample_parameters["n_radial"]);
    } catch (...) {
        // Use defaults
    }
    
    std::cout << "Parameters:" << std::endl;
    std::cout << "  NS radius: " << R_star << " code units" << std::endl;
    std::cout << "  Central density: " << rho_c << " code units" << std::endl;
    std::cout << "  Initial separation: " << separation << " code units" << std::endl;
    std::cout << "  Collision velocity: " << v_collision << " c" << std::endl;
    std::cout << "  Particles per radius: " << n_radial << std::endl;
    
    // Polytropic constant K from central conditions
    // For n=1 polytrope: P_c = K * rho_c^2, and mass-radius relation gives K
    // We'll set K such that the star is in hydrostatic equilibrium
    const real K = 0.1;  // Polytropic constant (adjust for desired compactness)
    
    std::vector<SPHParticle> particles;
    
    // Particle spacing
    const real dr = R_star / n_radial;
    
    // Initial smoothing length (will be iteratively refined)
    const real h_initial = 1.5 * dr;
    
    // Star centers
    vec_t center1 = {-separation / 2.0, 0.0};
    vec_t center2 = { separation / 2.0, 0.0};
    
    // Velocities (head-on collision along x-axis)
    vec_t vel1 = { v_collision, 0.0};
    vec_t vel2 = {-v_collision, 0.0};
    
    int particle_id = 0;
    
    // Function to create particles for one star
    auto create_star_particles = [&](const vec_t& center, const vec_t& velocity) {
        // Use concentric rings for 2D
        for (int ir = 0; ir <= n_radial; ++ir) {
            real r = ir * dr;
            
            if (ir == 0) {
                // Central particle
                real rho = lane_emden_n1_density(0, R_star, rho_c);
                if (rho < 1e-6 * rho_c) continue;
                
                real pres = polytropic_pressure(rho, K, gamma);
                real ene = pres / ((gamma - 1.0) * rho);
                
                // Area element for central particle
                real area = M_PI * (dr/2) * (dr/2);
                real mass = rho * area;
                
                // Baryon number (for SR-GSPH): nu = n * V = rho * area (in rest frame)
                real nu = rho * area;
                
                SPHParticle p;
                p.pos = center;
                p.vel = velocity;
                p.dens = rho;
                p.pres = pres;
                p.ene = ene;
                p.mass = mass;
                p.nu = nu;  // Baryon number for SR-GSPH
                p.sml = h_initial;  // Initial smoothing length
                p.id = particle_id++;
                particles.push_back(p);
            } else {
                // Ring of particles at radius r
                real rho = lane_emden_n1_density(r, R_star, rho_c);
                if (rho < 1e-6 * rho_c) continue;
                
                real pres = polytropic_pressure(rho, K, gamma);
                real ene = pres / ((gamma - 1.0) * rho);
                
                // Number of particles in this ring (circumference / spacing)
                int n_theta = std::max(6, static_cast<int>(2 * M_PI * r / dr));
                real dtheta = 2 * M_PI / n_theta;
                
                // Area element for ring particle
                real area = r * dtheta * dr;
                real mass = rho * area;
                
                // Baryon number for SR-GSPH
                real nu = rho * area;
                
                for (int it = 0; it < n_theta; ++it) {
                    real theta = it * dtheta;
                    
                    SPHParticle p;
                    p.pos[0] = center[0] + r * std::cos(theta);
                    p.pos[1] = center[1] + r * std::sin(theta);
                    p.vel = velocity;
                    p.dens = rho;
                    p.pres = pres;
                    p.ene = ene;
                    p.mass = mass;
                    p.nu = nu;  // Baryon number for SR-GSPH
                    p.sml = h_initial;  // Initial smoothing length
                    p.id = particle_id++;
                    particles.push_back(p);
                }
            }
        }
    };
    
    // Create both stars
    std::cout << "Creating star 1 at x = " << center1[0] << std::endl;
    create_star_particles(center1, vel1);
    int n_star1 = particles.size();
    
    std::cout << "Creating star 2 at x = " << center2[0] << std::endl;
    create_star_particles(center2, vel2);
    int n_star2 = particles.size() - n_star1;
    
    std::cout << "Star 1: " << n_star1 << " particles" << std::endl;
    std::cout << "Star 2: " << n_star2 << " particles" << std::endl;
    std::cout << "Total: " << particles.size() << " particles" << std::endl;
    
    // Calculate total mass for verification
    real total_mass = 0;
    for (const auto& p : particles) {
        total_mass += p.mass;
    }
    std::cout << "Total mass: " << total_mass << " code units" << std::endl;
    
    // Set particles in simulation
    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());
    
    std::cout << "=== NS Merger 2D initialization complete ===" << std::endl;
#endif
}

} // namespace sph
