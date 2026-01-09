#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"

#include <cmath>
#include <vector>

namespace sph
{

/**
 * Orszag-Tang Vortex: 2D MHD turbulence benchmark
 *
 * This is a standard 2D MHD test problem that develops complex shock
 * interactions and current sheets from smooth initial conditions.
 * There is NO analytic solution - results are compared with high-resolution
 * grid MHD codes.
 *
 * Initial conditions (periodic box [0, 1] x [0, 1]):
 *   rho = gamma^2 / (4*pi) = 25 / (36*pi) ~ 0.221
 *   P   = gamma / (4*pi)   = 5 / (12*pi) ~ 0.133
 *   vx  = -sin(2*pi*y)
 *   vy  = sin(2*pi*x)
 *   Bx  = -B0 * sin(2*pi*y)
 *   By  = B0 * sin(4*pi*x)
 *   Bz  = 0
 *
 * where B0 = 1/sqrt(4*pi) ~ 0.282 and gamma = 5/3
 *
 * Key features:
 * - Periodic boundaries in both x and y
 * - Develops shocks at t ~ 0.5
 * - Forms current sheets and magnetic reconnection regions
 * - Standard end time: t = 0.5 (shock formation) or t = 1.0 (developed turbulence)
 *
 * Reference: Orszag & Tang (1979), Dahlburg & Picone (1989)
 */
void Solver::make_mhd_orszag_tang()
{
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real gamma = m_param->physics.gamma;

    const real sqrt_4pi = std::sqrt(4.0 * M_PI);
    const real B0 = 1.0 / sqrt_4pi;

    const real rho0 = gamma * gamma / (4.0 * M_PI);
    const real P0 = gamma / (4.0 * M_PI);

    const real x_min = 0.0, x_max = 1.0;
    const real y_min = 0.0, y_max = 1.0;
    const real Lx = x_max - x_min;
    const real Ly = y_max - y_min;

    const int Nx = static_cast<int>(std::sqrt(static_cast<real>(N) * Lx / Ly));
    const int Ny = static_cast<int>(std::sqrt(static_cast<real>(N) * Ly / Lx));
    const int num_total = Nx * Ny;

    const real dx = Lx / Nx;
    const real dy = Ly / Ny;
    const real mass = rho0 * dx * dy;

    std::vector<SPHParticle> particles(num_total);

    int idx = 0;
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            const real x = x_min + (ix + 0.5) * dx;
            const real y = y_min + (iy + 0.5) * dy;

            const real Bx = -B0 * std::sin(2.0 * M_PI * y);
            const real By = B0 * std::sin(4.0 * M_PI * x);

            auto& p = particles[idx];
#if DIM >= 2
            const real vx = -std::sin(2.0 * M_PI * y);
            const real vy = std::sin(2.0 * M_PI * x);
#if DIM == 2
            p.pos = vec_t(x, y);
            p.vel = vec_t(vx, vy);
#else
            p.pos = vec_t(x, y, 0.0);
            p.vel = vec_t(vx, vy, 0.0);
#endif
#else
            (void)x; (void)y;
            THROW_ERROR("Orszag-Tang vortex requires DIM >= 2");
#endif
            p.dens = rho0;
            p.pres = P0;
            p.mass = mass;
            p.ene = P0 / ((gamma - 1.0) * rho0);
            p.B = vec3_t(Bx, By, 0.0);
            p.id = idx;

            ++idx;
        }
    }

    particles.resize(idx);

    std::cout << "[MHD Orszag-Tang Vortex]" << std::endl;
    std::cout << "  Resolution: " << Nx << " x " << Ny << " = " << idx << " particles" << std::endl;
    std::cout << "  Domain: [" << x_min << ", " << x_max << "] x [" << y_min << ", " << y_max << "]" << std::endl;
    std::cout << "  dx = " << dx << ", dy = " << dy << std::endl;
    std::cout << "  rho0 = " << rho0 << ", P0 = " << P0 << std::endl;
    std::cout << "  B0 = " << B0 << std::endl;
    std::cout << "  gamma = " << gamma << std::endl;

    m_sim->set_particles(particles);
    m_sim->set_particle_num(particles.size());
}

}
