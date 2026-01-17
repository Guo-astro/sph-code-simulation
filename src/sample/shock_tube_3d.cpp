// 3D Shock Tube (Sod shock tube extended to 3D)
// Classic Riemann problem in 3D: left and right states separated at x_discontinuity
// Domain: [rangeMin[0], rangeMax[0]] x [0,Ly] x [0,Lz]
// Left state: rho=1.0, P=1.0, v=0
// Right state: rho=0.125, P=0.1, v=0
//
// CUBIC LATTICE approach (isotropic on both sides):
// - Left side: cubic lattice with spacing d_L = dx_L = dy_L = dz_L
// - Right side: cubic lattice with spacing d_R = dx_R = dy_R = dz_R
// - For equal mass: d_R = d_L * (rho_L / rho_R)^(1/3)
//   With rho_L/rho_R = 8, we get d_R = 2 * d_L
//
// This ensures SPH kernels see isotropic neighbor distributions on both sides.

#include "../../include/solver.hpp"
#include "../../include/simulation.hpp"
#include "../../include/particle.hpp"
#include "../../include/exception.hpp"
#include "../../include/parameters.hpp"
#include <cmath>
#include <iostream>

namespace sph {

void Solver::make_shock_tube_3d()
{
#if DIM != 3
    THROW_ERROR("DIM != 3");
#else
    // Get domain size from rangeMin/rangeMax parameters
    real x_min = m_param->periodic.range_min[0];
    real x_max = m_param->periodic.range_max[0];
    real Lx = x_max - x_min;

    real Ly = 0.2;
    real Lz = 0.2;
    try {
        Ly = boost::any_cast<real>(m_sample_parameters["Ly"]);
        Lz = boost::any_cast<real>(m_sample_parameters["Lz"]);
    } catch(...) {}

    // Physical parameters
    const real gamma = m_param->physics.gamma;

    // Discontinuity at midpoint of domain
    const real x_discontinuity = (x_min + x_max) / 2.0;

    // Left state (x < x_discontinuity)
    const real rho_left = 1.0;
    const real pres_left = 1.0;

    // Right state (x >= x_discontinuity)
    const real rho_right = 0.125;
    const real pres_right = 0.1;

    // ================================================================
    // CUBIC LATTICE: ISOTROPIC ON BOTH SIDES
    // ================================================================
    // For proper 3D SPH, both sides need cubic (isotropic) lattice.
    // Equal mass condition with cubic lattice:
    //   m = rho_L * d_L^3 = rho_R * d_R^3
    //   d_R = d_L * (rho_L / rho_R)^(1/3)
    // With rho_L/rho_R = 8, d_R = 2 * d_L
    // ================================================================

    const real density_ratio = rho_left / rho_right;  // = 8
    const real spacing_ratio = std::cbrt(density_ratio);  // = 2

    // Get resolution from Nx parameter (particles in x-direction for LEFT side)
    int Nx_param;
    try {
        Nx_param = boost::any_cast<int>(m_sample_parameters["Nx"]);
    } catch(...) {
        Nx_param = 100;  // Default: 100 particles in left half
    }

    // Left side length and right side length
    const real Lx_left = x_discontinuity - x_min;
    const real Lx_right = x_max - x_discontinuity;

    // LEFT SIDE: cubic lattice with spacing d_L
    const int Nx_left = std::max(1, static_cast<int>(Nx_param * (Lx_left / 0.5) + 0.5));
    const real d_left = Lx_left / Nx_left;  // Cubic spacing on left
    const int Ny_left = std::max(1, static_cast<int>(Ly / d_left + 0.5));
    const int Nz_left = std::max(1, static_cast<int>(Lz / d_left + 0.5));
    const real dy_left = Ly / Ny_left;
    const real dz_left = Lz / Nz_left;

    // RIGHT SIDE: cubic lattice with spacing d_R = d_L * spacing_ratio
    const real d_right = d_left * spacing_ratio;  // 2x larger spacing
    const int Nx_right = std::max(1, static_cast<int>(Lx_right / d_right + 0.5));
    const int Ny_right = std::max(1, static_cast<int>(Ly / d_right + 0.5));
    const int Nz_right = std::max(1, static_cast<int>(Lz / d_right + 0.5));
    const real dx_right = Lx_right / Nx_right;
    const real dy_right = Ly / Ny_right;
    const real dz_right = Lz / Nz_right;

    // Particle mass (same for all particles)
    // Based on left side: m = rho_L * d_L^3
    const real particle_mass = rho_left * d_left * dy_left * dz_left;

    // Verify mass consistency on right side
    const real mass_right_check = rho_right * dx_right * dy_right * dz_right;

    const int num_left = Nx_left * Ny_left * Nz_left;
    const int num_right = Nx_right * Ny_right * Nz_right;
    const int num_total = num_left + num_right;

    std::vector<SPHParticle> p;
    p.reserve(num_total);

    std::cout << "\n=== SHOCK TUBE 3D SETUP (CUBIC LATTICE - ISOTROPIC) ===" << std::endl;
    std::cout << "Domain: [" << x_min << ", " << x_max << "] x [0, " << Ly << "] x [0, " << Lz << "]" << std::endl;
    std::cout << "Discontinuity at x = " << x_discontinuity << std::endl;
    std::cout << "Density ratio: " << density_ratio << ", Spacing ratio: " << spacing_ratio << std::endl;
    std::cout << "Left:  d=" << d_left << " (Nx=" << Nx_left << " Ny=" << Ny_left << " Nz=" << Nz_left << ")" << std::endl;
    std::cout << "       actual: dx=" << d_left << " dy=" << dy_left << " dz=" << dz_left << std::endl;
    std::cout << "Right: d=" << d_right << " (Nx=" << Nx_right << " Ny=" << Ny_right << " Nz=" << Nz_right << ")" << std::endl;
    std::cout << "       actual: dx=" << dx_right << " dy=" << dy_right << " dz=" << dz_right << std::endl;
    std::cout << "Particle mass: " << particle_mass << std::endl;
    std::cout << "Mass check (right): " << mass_right_check << " (ratio: " << mass_right_check/particle_mass << ")" << std::endl;
    std::cout << "Total particles: " << num_total << " (left: " << num_left << ", right: " << num_right << ")" << std::endl;
    std::cout << "========================================================\n" << std::endl;

    int id = 0;

    // Left region: cubic lattice with spacing d_left
    for(int iz = 0; iz < Nz_left; ++iz) {
        real z = dz_left * (0.5 + iz);
        for(int iy = 0; iy < Ny_left; ++iy) {
            real y = dy_left * (0.5 + iy);
            for(int ix = 0; ix < Nx_left; ++ix) {
                real x = x_min + d_left * (0.5 + ix);

                SPHParticle p_i;
                p_i.pos[0] = x;
                p_i.pos[1] = y;
                p_i.pos[2] = z;
                p_i.vel[0] = 0.0;
                p_i.vel[1] = 0.0;
                p_i.vel[2] = 0.0;
                p_i.dens = rho_left;
                p_i.pres = pres_left;
                p_i.mass = particle_mass;
                p_i.ene = pres_left / ((gamma - 1.0) * rho_left);
                p_i.id = id++;
                p.push_back(p_i);
            }
        }
    }

    // Right region: cubic lattice with spacing d_right (2x larger)
    for(int iz = 0; iz < Nz_right; ++iz) {
        real z = dz_right * (0.5 + iz);
        for(int iy = 0; iy < Ny_right; ++iy) {
            real y = dy_right * (0.5 + iy);
            for(int ix = 0; ix < Nx_right; ++ix) {
                real x = x_discontinuity + dx_right * (0.5 + ix);

                SPHParticle p_i;
                p_i.pos[0] = x;
                p_i.pos[1] = y;
                p_i.pos[2] = z;
                p_i.vel[0] = 0.0;
                p_i.vel[1] = 0.0;
                p_i.vel[2] = 0.0;
                p_i.dens = rho_right;
                p_i.pres = pres_right;
                p_i.mass = particle_mass;
                p_i.ene = pres_right / ((gamma - 1.0) * rho_right);
                p_i.id = id++;
                p.push_back(p_i);
            }
        }
    }

    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
#endif
}

} // namespace sph
