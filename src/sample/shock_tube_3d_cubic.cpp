// 3D Shock Tube (Sod shock tube extended to 3D)
// Using CUBIC lattice for particle placement
// Domain: from rangeMin to rangeMax (set in config)
// Discontinuity at x=0 (domain should be centered at origin)
//
// Left state: rho=1.0, P=1.0, v=0
// Right state: rho=0.125, P=0.1, v=0
// Spacing ratio: (rho_L/rho_R)^(1/3) = 2.0

#include "../../include/solver.hpp"
#include "../../include/simulation.hpp"
#include "../../include/particle.hpp"
#include "../../include/exception.hpp"
#include "../../include/parameters.hpp"
#include <cmath>
#include <iostream>
#include <vector>

namespace sph {

// Cubic lattice particle position formula
static vec_t get_cubic_pos(real dr, int i, int j, int k, const vec_t& box_min)
{
    vec_t r;
    r[0] = (i + 0.5) * dr + box_min[0];
    r[1] = (j + 0.5) * dr + box_min[1];
    r[2] = (k + 0.5) * dr + box_min[2];
    return r;
}

void Solver::make_shock_tube_3d_cubic()
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

    // Spacing factor: fact = (rho_L/rho_R)^(1/3) = 2.0
    const real fact = std::cbrt(rho_L / rho_R);

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

    // Resolution parameter (particles in x-direction on left side)
    int Nx_left = 100;
    try {
        Nx_left = boost::any_cast<int>(m_sample_parameters["Nx"]);
    } catch(...) {}

    // Left side spacing: fit Nx_left particles in half the domain
    const real dr_left = (Lx / 2.0) / Nx_left;  // Cubic spacing
    const real dr_right = dr_left * fact;  // Right side has larger spacing

    std::cout << "\n=== SHOCK TUBE 3D SETUP (CUBIC LATTICE) ===" << std::endl;
    std::cout << "Domain: x in [" << x_min << ", " << x_max << "]" << std::endl;
    std::cout << "        y in [" << y_min << ", " << y_max << "]" << std::endl;
    std::cout << "        z in [" << z_min << ", " << z_max << "]" << std::endl;
    std::cout << "Discontinuity at x = 0" << std::endl;
    std::cout << "dr_left = " << dr_left << ", dr_right = " << dr_right << std::endl;
    std::cout << "Spacing factor: " << fact << " (rho_L/rho_R = " << rho_L/rho_R << ")" << std::endl;

    std::vector<SPHParticle> particles;

    // LEFT SIDE: Cubic lattice with spacing dr_left
    {
        vec_t box_min_left = {x_min, y_min, z_min};
        vec_t box_max_left = {0.0, y_max, z_max};

        int nx = static_cast<int>(std::ceil((box_max_left[0] - box_min_left[0]) / dr_left));
        int ny = static_cast<int>(std::ceil((box_max_left[1] - box_min_left[1]) / dr_left));
        int nz = static_cast<int>(std::ceil((box_max_left[2] - box_min_left[2]) / dr_left));

        std::cout << "Left side (dr=" << dr_left << "): lattice (" << nx << ", " << ny << ", " << nz << ")" << std::endl;

        for(int k = 0; k < nz; ++k) {
            for(int j = 0; j < ny; ++j) {
                for(int i = 0; i < nx; ++i) {
                    vec_t r = get_cubic_pos(dr_left, i, j, k, box_min_left);

                    if(r[0] >= box_min_left[0] && r[0] < box_max_left[0] &&
                       r[1] >= box_min_left[1] && r[1] < box_max_left[1] &&
                       r[2] >= box_min_left[2] && r[2] < box_max_left[2]) {
                        SPHParticle p;
                        p.pos = r;
                        p.vel = vec_t(0.0);
                        p.dens = rho_L;
                        p.pres = P_L;
                        p.ene = u_L;
                        p.id = static_cast<int>(particles.size());
                        particles.push_back(p);
                    }
                }
            }
        }
    }

    int num_left = particles.size();
    std::cout << "Left side particles: " << num_left << std::endl;

    // RIGHT SIDE: Cubic lattice with spacing dr_right
    {
        vec_t box_min_right = {0.0, y_min, z_min};
        vec_t box_max_right = {x_max, y_max, z_max};

        int nx = static_cast<int>(std::ceil((box_max_right[0] - box_min_right[0]) / dr_right));
        int ny = static_cast<int>(std::ceil((box_max_right[1] - box_min_right[1]) / dr_right));
        int nz = static_cast<int>(std::ceil((box_max_right[2] - box_min_right[2]) / dr_right));

        std::cout << "Right side (dr=" << dr_right << "): lattice (" << nx << ", " << ny << ", " << nz << ")" << std::endl;

        for(int k = 0; k < nz; ++k) {
            for(int j = 0; j < ny; ++j) {
                for(int i = 0; i < nx; ++i) {
                    vec_t r = get_cubic_pos(dr_right, i, j, k, box_min_right);

                    if(r[0] >= box_min_right[0] && r[0] < box_max_right[0] &&
                       r[1] >= box_min_right[1] && r[1] < box_max_right[1] &&
                       r[2] >= box_min_right[2] && r[2] < box_max_right[2]) {
                        SPHParticle p;
                        p.pos = r;
                        p.vel = vec_t(0.0);
                        p.dens = rho_R;
                        p.pres = P_R;
                        p.ene = u_R;
                        p.id = static_cast<int>(particles.size());
                        particles.push_back(p);
                    }
                }
            }
        }
    }

    int num_right = particles.size() - num_left;
    std::cout << "Right side particles: " << num_right << std::endl;

    // Set particle mass (equal mass for all)
    const real vol_left = (Lx / 2.0) * Ly * Lz;
    const real vol_right = (Lx / 2.0) * Ly * Lz;
    const real total_mass = rho_L * vol_left + rho_R * vol_right;
    const real pmass = total_mass / particles.size();

    std::cout << "Total mass: " << total_mass << std::endl;
    std::cout << "Particle mass: " << pmass << std::endl;
    std::cout << "Total particles: " << particles.size() << std::endl;
    std::cout << "============================================\n" << std::endl;

    for(auto& p : particles) {
        p.mass = pmass;
    }

    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());
#endif
}

} // namespace sph
