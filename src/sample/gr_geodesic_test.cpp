/**
 * GR Geodesic Test
 *
 * Tests the gravitational source terms by simulating test particle motion
 * in Schwarzschild spacetime. This is equivalent to solving the geodesic equation.
 *
 * Based on Liptai & Price (2019) Section 6.2:
 * - Radial infall from rest
 * - Circular orbits
 *
 * The exact solution for radial infall starting at rest from r_0:
 *   v^r(r) = (1 - 2M/r) / sqrt(1 - 2M/r_0) * sqrt(2M * (1/r - 1/r_0))
 *
 * For circular orbits at radius r:
 *   Omega = sqrt(M / r^3)
 */

#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"
#include "logger.hpp"
#include <cmath>
#include <string>

namespace sph
{

void Solver::make_gr_geodesic_test()
{
    // Black hole mass
    real M = 1.0;
    if (m_sample_parameters.count("blackHoleMass")) {
        M = boost::any_cast<real>(m_sample_parameters["blackHoleMass"]);
    }

    // Test type: "radial_infall" or "circular_orbit"
    std::string test_type = "radial_infall";
    if (m_sample_parameters.count("testType")) {
        test_type = boost::any_cast<std::string>(m_sample_parameters["testType"]);
    }

    // Number of test particles
    int N = 10;
    if (m_sample_parameters.count("N")) {
        N = boost::any_cast<int>(m_sample_parameters["N"]);
    }

    const real gamma = m_param->physics.gamma;

    WRITE_LOG << "GR Geodesic Test:";
    WRITE_LOG << "  Black hole mass M = " << M;
    WRITE_LOG << "  Test type: " << test_type;
    WRITE_LOG << "  Number of particles: " << N;

    auto& particles = m_sim->get_particles();
    particles.clear();
    particles.reserve(N);

    if (test_type == "radial_infall") {
        // Radial infall from rest at various radii
        // Domain: r_0 in [4M, 40M] evenly spaced
        real r_min = 4.0 * M;
        real r_max = 40.0 * M;
        if (m_sample_parameters.count("rMin")) {
            r_min = boost::any_cast<real>(m_sample_parameters["rMin"]);
        }
        if (m_sample_parameters.count("rMax")) {
            r_max = boost::any_cast<real>(m_sample_parameters["rMax"]);
        }

        WRITE_LOG << "  Radial infall test: r_0 in [" << r_min << ", " << r_max << "]";

        for (int i = 0; i < N; ++i) {
            SPHParticle p;
            p.id = i;

            // Initial radius (evenly spaced)
            const real r0 = r_min + (r_max - r_min) * i / (N - 1);
            p.pos[0] = r0;
#if DIM >= 2
            p.pos[1] = 0.0;
#endif
#if DIM == 3
            p.pos[2] = 0.0;
#endif

            // Start at rest (v = 0)
            p.vel[0] = 0.0;
#if DIM >= 2
            p.vel[1] = 0.0;
#endif
#if DIM == 3
            p.vel[2] = 0.0;
#endif
            p.vel_t = 0.0;

            // Negligible pressure (geodesic = pressure-less)
            p.pres = 1.0e-10;
            p.dens = 1.0;
            p.ene = p.pres / ((gamma - 1.0) * p.dens);
            p.sound = std::sqrt(gamma * p.pres / p.dens);
            p.gamma_lor = 1.0;
            p.N = p.dens;
            p.nu = 1.0;
            p.mass = 1.0 / N;  // Equal mass particles
            // Smoothing length should be large enough for kernel overlap
            // With particles spaced ~4M apart, need sml > spacing for neighbors
            const real dr = (r_max - r_min) / (N - 1);  // Particle spacing
            p.sml = 2.0 * dr;  // Ensure neighbor overlap
            p.gradh = 1.0;     // Initialize gradh to 1.0
            p.is_ghost = false;

            particles.push_back(p);
        }

    } else if (test_type == "circular_orbit") {
        // Circular orbit test - create an annulus of particles
        real r_orbit = 10.0 * M;
        if (m_sample_parameters.count("rOrbit")) {
            r_orbit = boost::any_cast<real>(m_sample_parameters["rOrbit"]);
        }

        // Width of the annulus (for SPH neighbor finding)
        real dr = 0.5 * M;  // Annulus half-width
        if (m_sample_parameters.count("annulusWidth")) {
            dr = boost::any_cast<real>(m_sample_parameters["annulusWidth"]) / 2.0;
        }

        WRITE_LOG << "  Circular orbit test: r = " << r_orbit;
        WRITE_LOG << "  Annulus width: " << (2.0 * dr);

        // Number of radial layers
        int n_radial = 5;
        if (m_sample_parameters.count("nRadial")) {
            n_radial = boost::any_cast<int>(m_sample_parameters["nRadial"]);
        }

        // Distribute N particles across n_radial layers
        int particles_per_layer = N / n_radial;
        real layer_spacing = (n_radial > 1) ? (2.0 * dr) / (n_radial - 1) : 0.0;

        int particle_id = 0;
        for (int layer = 0; layer < n_radial; ++layer) {
            // Radius for this layer
            real r = r_orbit - dr + layer * layer_spacing;
            if (n_radial == 1) r = r_orbit;

            // Circular orbit in Schwarzschild coordinates
            // Using Newtonian orbital velocity: v = sqrt(M/r)
            // The GR source term with omega (gradh) correction should provide proper balance
            const real Omega = std::sqrt(M / (r * r * r));
            const real v_phi = Omega * r;  // = sqrt(M/r)

            // Number of particles in this layer (adjust for circumference)
            int n_phi = particles_per_layer;
            if (layer == 0) {
                WRITE_LOG << "  Orbital velocity v_phi = " << v_phi << " at r = " << r;
            }

            for (int i = 0; i < n_phi; ++i) {
                SPHParticle p;
                p.id = particle_id++;

                // Initial position (equally spaced in phi, offset alternate layers)
                real phi_offset = (layer % 2 == 0) ? 0.0 : M_PI / n_phi;
                const real phi = 2.0 * M_PI * i / n_phi + phi_offset;
                p.pos[0] = r * std::cos(phi);
#if DIM >= 2
                p.pos[1] = r * std::sin(phi);
#endif
#if DIM == 3
                p.pos[2] = 0.0;
#endif

                // Tangential velocity (perpendicular to radius)
                p.vel[0] = -v_phi * std::sin(phi);
#if DIM >= 2
                p.vel[1] = v_phi * std::cos(phi);
#endif
#if DIM == 3
                p.vel[2] = 0.0;
#endif
                p.vel_t = 0.0;

                // Set small but non-negligible density for SPH
                p.dens = 1.0;
                p.pres = 1.0e-6;  // Small pressure
                p.ene = p.pres / ((gamma - 1.0) * p.dens);
                p.sound = std::sqrt(gamma * p.pres / p.dens);

                // Lorentz factor for circular orbit
                const real v2 = v_phi * v_phi;
                p.gamma_lor = 1.0 / std::sqrt(1.0 - v2);

                p.N = p.dens * p.gamma_lor;
                p.nu = 1.0;
                p.mass = 1.0 / (n_radial * n_phi);  // Equal mass particles
                p.sml = 1.5 * layer_spacing;  // Smoothing length based on spacing
                if (p.sml < 0.1 * M) p.sml = 0.1 * M;
                p.is_ghost = false;

                particles.push_back(p);
            }
        }

        WRITE_LOG << "  Created " << particles.size() << " particles in annulus";

    } else {
        THROW_ERROR("Unknown geodesic test type: " + test_type);
    }

    m_sim->set_particle_num(static_cast<int>(particles.size()));

    // Store parameters for output
    m_sample_parameters["M"] = M;
    m_sample_parameters["test_type"] = test_type;

    WRITE_LOG << "Geodesic test initialized with " << particles.size() << " particles";
}

} // namespace sph
