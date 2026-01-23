// 3D Shock Tube - Column-based approach (Kitajima style)
// Using uniform HCP lattice with NON-EQUAL MASS particles
// Domain: from rangeMin to rangeMax (set in config)
// Discontinuity at x=0 (domain should be centered at origin)
//
// Key difference from equal-mass approach:
// - Particles have UNIFORM SPACING throughout the domain
// - Density difference achieved by varying PARTICLE MASS
// - Left side: higher mass particles (rho=1.0)
// - Right side: lower mass particles (rho=0.125)
//
// Left state: rho=1.0, P=1.0, v=0
// Right state: rho=0.125, P=0.1, v=0

#include "../../include/solver.hpp"
#include "../../include/simulation.hpp"
#include "../../include/particle.hpp"
#include "../../include/exception.hpp"
#include "../../include/parameters.hpp"
#include <cmath>
#include <iostream>
#include <vector>

namespace sph {

// HCP lattice particle position formula (same as shock_tube_3d.cpp)
static vec_t get_hcp_pos_column(real dr, int i, int j, int k, const vec_t& box_min)
{
    vec_t r;
    r[0] = 2 * i + ((j + k) % 2);
    r[1] = std::sqrt(3.0) * (j + (1.0/3.0) * (k % 2));
    r[2] = 2.0 * std::sqrt(6.0) * k / 3.0;
    r[0] *= dr;
    r[1] *= dr;
    r[2] *= dr;
    r[0] += box_min[0];
    r[1] += box_min[1];
    r[2] += box_min[2];
    return r;
}

void Solver::make_shock_tube_3d_column()
{
#if DIM != 3
    THROW_ERROR("DIM != 3");
#else
    // Physical parameters
    const real gamma = m_param->physics.gamma;

    const real rho_L = 1.0;      // Left density
    const real rho_R = 0.125;    // Right density
    const real P_L = 1.0;        // Left pressure
    const real P_R = 0.1;        // Right pressure

    // Internal energy from pressure: u = P / ((gamma-1) * rho)
    const real u_L = P_L / ((gamma - 1.0) * rho_L);
    const real u_R = P_R / ((gamma - 1.0) * rho_R);

    // Get domain from config (should be centered at origin)
    const real x_min = m_param->periodic.range_min[0];
    const real x_max = m_param->periodic.range_max[0];
    const real y_min = m_param->periodic.range_min[1];
    const real y_max = m_param->periodic.range_max[1];
    const real z_min = m_param->periodic.range_min[2];
    const real z_max = m_param->periodic.range_max[2];

    const real Lx = x_max - x_min;
    const real Ly = y_max - y_min;
    const real Lz = z_max - z_min;

    // Resolution parameter (particles in x-direction total)
    int Nx = 100;
    try {
        Nx = boost::any_cast<int>(m_sample_parameters["Nx"]);
    } catch(...) {}

    // UNIFORM spacing throughout entire domain
    // For HCP, effective x-spacing is 2*dr
    const real dr = Lx / (2.0 * Nx);

    std::cout << "\n=== SHOCK TUBE 3D SETUP (COLUMN/NON-EQUAL MASS) ===" << std::endl;
    std::cout << "Domain: x in [" << x_min << ", " << x_max << "]" << std::endl;
    std::cout << "        y in [" << y_min << ", " << y_max << "]" << std::endl;
    std::cout << "        z in [" << z_min << ", " << z_max << "]" << std::endl;
    std::cout << "Discontinuity at x = 0" << std::endl;
    std::cout << "UNIFORM dr = " << dr << " (same for entire domain)" << std::endl;
    std::cout << "Density ratio achieved via MASS variation" << std::endl;

    std::vector<SPHParticle> particles;

    // Compute lattice counts for entire domain
    vec_t box_min_vec = {x_min, y_min, z_min};
    vec_t box_max_vec = {x_max, y_max, z_max};
    vec_t box_dim = box_max_vec - box_min_vec;

    int ix = static_cast<int>(std::ceil(box_dim[0] / (2.0 * dr)));
    int iy = static_cast<int>(std::ceil(box_dim[1] / (std::sqrt(3.0) * dr)));
    int iz = static_cast<int>(std::ceil(box_dim[2] / (2.0 * std::sqrt(6.0) * dr / 3.0)));

    std::cout << "HCP lattice: (" << ix << ", " << iy << ", " << iz << ")" << std::endl;

    // Compute particle volume (same for all particles due to uniform spacing)
    // For HCP: volume per particle = dr^3 * sqrt(2)
    const real vol_per_particle = dr * dr * dr * std::sqrt(2.0);

    // Compute masses based on local density
    const real mass_L = rho_L * vol_per_particle;
    const real mass_R = rho_R * vol_per_particle;

    std::cout << "Volume per particle: " << vol_per_particle << std::endl;
    std::cout << "Mass left (rho=" << rho_L << "): " << mass_L << std::endl;
    std::cout << "Mass right (rho=" << rho_R << "): " << mass_R << std::endl;
    std::cout << "Mass ratio: " << mass_L / mass_R << " (should be " << rho_L/rho_R << ")" << std::endl;

    int count_left = 0, count_right = 0;

    for(int k = 0; k < iz; ++k) {
        for(int j = 0; j < iy; ++j) {
            for(int i = 0; i < ix; ++i) {
                vec_t r = get_hcp_pos_column(dr, i, j, k, box_min_vec);

                // Only add if inside the domain
                if(r[0] >= box_min_vec[0] && r[0] < box_max_vec[0] &&
                   r[1] >= box_min_vec[1] && r[1] < box_max_vec[1] &&
                   r[2] >= box_min_vec[2] && r[2] < box_max_vec[2]) {

                    SPHParticle p;
                    p.pos = r;
                    p.vel = vec_t(0.0);
                    p.id = static_cast<int>(particles.size());

                    // Assign properties based on position (left or right of discontinuity)
                    if (r[0] < 0.0) {
                        // LEFT SIDE: high density, high pressure
                        p.dens = rho_L;
                        p.pres = P_L;
                        p.ene = u_L;
                        p.mass = mass_L;
                        count_left++;
                    } else {
                        // RIGHT SIDE: low density, low pressure
                        p.dens = rho_R;
                        p.pres = P_R;
                        p.ene = u_R;
                        p.mass = mass_R;
                        count_right++;
                    }

                    particles.push_back(p);
                }
            }
        }
    }

    // Compute total mass for verification
    real total_mass_actual = 0.0;
    for (const auto& p : particles) {
        total_mass_actual += p.mass;
    }

    // Expected total mass
    const real vol_left = (Lx / 2.0) * Ly * Lz;
    const real vol_right = (Lx / 2.0) * Ly * Lz;
    const real expected_mass = rho_L * vol_left + rho_R * vol_right;

    std::cout << "Left side particles: " << count_left << std::endl;
    std::cout << "Right side particles: " << count_right << std::endl;
    std::cout << "Total particles: " << particles.size() << std::endl;
    std::cout << "Total mass (actual): " << total_mass_actual << std::endl;
    std::cout << "Total mass (expected): " << expected_mass << std::endl;
    std::cout << "============================================\n" << std::endl;

    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());
#endif
}

} // namespace sph
