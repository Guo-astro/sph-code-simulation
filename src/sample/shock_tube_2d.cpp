// 2D Shock Tube (Sod shock tube extended to 2D)
// Classic Riemann problem in 2D: left and right states separated at x=0.5
// Domain: [0,1] x [0,Ly]
// Left state: rho=1.0, P=1.0, v=0
// Right state: rho=0.125, P=0.1, v=0
//
// EQUAL MASS approach (following shamrock/PySPH convention):
// All particles have the same mass m.
// Particle spacing varies with density:
//   - In 2D: spacing ratio = (rho_L / rho_R)^(1/2) = sqrt(8) ≈ 2.83
//   - Left side (high density): spacing = dr
//   - Right side (low density): spacing = dr * fact, where fact = (rho_L/rho_R)^(1/2)

#include "../../include/solver.hpp"
#include "../../include/simulation.hpp"
#include "../../include/particle.hpp"
#include "../../include/exception.hpp"
#include "../../include/parameters.hpp"
#include <cmath>
#include <iostream>

namespace sph {

void Solver::make_shock_tube_2d()
{
#if DIM != 2
    THROW_ERROR("DIM != 2");
#else
    // Get domain size from parameters
    real Lx = 1.0;
    real Ly = 0.2;
    try {
        auto rangeMin = boost::any_cast<std::vector<real>>(m_sample_parameters["rangeMin"]);
        auto rangeMax = boost::any_cast<std::vector<real>>(m_sample_parameters["rangeMax"]);
        Lx = rangeMax[0] - rangeMin[0];
        Ly = rangeMax[1] - rangeMin[1];
    } catch(...) {}

    // Physical parameters
    const real gamma = m_param->physics.gamma;
    const real x_discontinuity = Lx / 2.0;  // Discontinuity at midpoint

    // Left state (x < x_discontinuity)
    const real rho_left = 1.0;
    const real pres_left = 1.0;

    // Right state (x >= x_discontinuity)
    const real rho_right = 0.125;
    const real pres_right = 0.1;

    // ================================================================
    // EQUAL MASS PARTICLE SETUP (shamrock/PySPH convention)
    // ================================================================
    // For 2D: spacing ratio = (rho_L / rho_R)^(1/2)
    // This ensures equal mass particles when:
    //   m = rho * dx * dy (2D area element)
    //   m_L = rho_L * dx_L * dy = m_R = rho_R * dx_R * dy
    //   => dx_R / dx_L = rho_L / rho_R  (if dy is same)
    //
    // But for isotropic spacing (dx = dy = dr):
    //   m = rho * dr^2
    //   m_L = rho_L * dr_L^2 = m_R = rho_R * dr_R^2
    //   => dr_R / dr_L = sqrt(rho_L / rho_R)
    // ================================================================

    const real density_ratio = rho_left / rho_right;  // = 8
    const real spacing_factor = std::pow(density_ratio, 0.5);  // = sqrt(8) ≈ 2.83 for 2D

    // Get resolution from Nx parameter (particles in x-direction for LEFT side)
    int Nx_left;
    try {
        Nx_left = boost::any_cast<int>(m_sample_parameters["Nx"]);
    } catch(...) {
        Nx_left = 200;  // Default
    }

    // Base spacing for left side (high density region)
    const real dr_left = (Lx / 2.0) / Nx_left;

    // Spacing for right side (low density region)
    const real dr_right = dr_left * spacing_factor;

    // Number of particles in each direction
    const int Ny_left = static_cast<int>(Ly / dr_left + 0.5);
    const int Nx_right = static_cast<int>((Lx / 2.0) / dr_right + 0.5);
    const int Ny_right = static_cast<int>(Ly / dr_right + 0.5);

    // Actual spacing after rounding
    const real dx_left = (Lx / 2.0) / Nx_left;
    const real dy_left = Ly / Ny_left;
    const real dx_right = (Lx / 2.0) / Nx_right;
    const real dy_right = Ly / Ny_right;

    // Particle mass (same for all particles)
    // Using left side as reference: m = rho_L * dx_L * dy_L
    const real particle_mass = rho_left * dx_left * dy_left;

    // Verify mass consistency
    const real mass_right_check = rho_right * dx_right * dy_right;

    const int num_left = Nx_left * Ny_left;
    const int num_right = Nx_right * Ny_right;
    const int num_total = num_left + num_right;

    std::vector<SPHParticle> p;
    p.reserve(num_total);

    std::cout << "\n=== SHOCK TUBE 2D SETUP (EQUAL MASS - SHAMROCK STYLE) ===" << std::endl;
    std::cout << "Domain: [0, " << Lx << "] x [0, " << Ly << "]" << std::endl;
    std::cout << "Density ratio (rho_L/rho_R): " << density_ratio << std::endl;
    std::cout << "Spacing factor (2D): " << spacing_factor << " (= sqrt(" << density_ratio << "))" << std::endl;
    std::cout << "Left side:  Nx=" << Nx_left << ", Ny=" << Ny_left
              << ", dx=" << dx_left << ", dy=" << dy_left << std::endl;
    std::cout << "Right side: Nx=" << Nx_right << ", Ny=" << Ny_right
              << ", dx=" << dx_right << ", dy=" << dy_right << std::endl;
    std::cout << "Particle mass: " << particle_mass << std::endl;
    std::cout << "Mass check (right): " << mass_right_check << " (ratio: " << mass_right_check/particle_mass << ")" << std::endl;
    std::cout << "Total particles: " << num_total << " (left: " << num_left << ", right: " << num_right << ")" << std::endl;
    std::cout << "=========================================================\n" << std::endl;

    int id = 0;

    // Left region: x in [0, Lx/2)
    for(int iy = 0; iy < Ny_left; ++iy) {
        real y = dy_left * (0.5 + iy);
        for(int ix = 0; ix < Nx_left; ++ix) {
            real x = dx_left * (0.5 + ix);

            SPHParticle p_i;
            p_i.pos[0] = x;
            p_i.pos[1] = y;
            p_i.vel[0] = 0.0;
            p_i.vel[1] = 0.0;
            p_i.dens = rho_left;
            p_i.pres = pres_left;
            p_i.mass = particle_mass;
            p_i.ene = pres_left / ((gamma - 1.0) * rho_left);
            p_i.id = id++;
            p.push_back(p_i);
        }
    }

    // Right region: x in [Lx/2, Lx]
    for(int iy = 0; iy < Ny_right; ++iy) {
        real y = dy_right * (0.5 + iy);
        for(int ix = 0; ix < Nx_right; ++ix) {
            real x = x_discontinuity + dx_right * (0.5 + ix);

            SPHParticle p_i;
            p_i.pos[0] = x;
            p_i.pos[1] = y;
            p_i.vel[0] = 0.0;
            p_i.vel[1] = 0.0;
            p_i.dens = rho_right;
            p_i.pres = pres_right;
            p_i.mass = particle_mass;
            p_i.ene = pres_right / ((gamma - 1.0) * rho_right);
            p_i.id = id++;
            p.push_back(p_i);
        }
    }

    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
#endif
}

} // namespace sph
