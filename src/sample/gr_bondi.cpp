/**
 * Bondi Accretion in Schwarzschild Spacetime
 *
 * Steady-state spherical accretion onto a Schwarzschild black hole.
 * Tests GR-GSPH in curved spacetime with gravitational source terms.
 *
 * Based on Liptai & Price (2019) MNRAS 485, 819 - Figures 14 and 18.
 *
 * Test types:
 *   "geodesic" - Pressure-less radial infall (Fig 14)
 *                v^r = -sqrt(2M/r) * (1 - 2M/r)
 *                n = n_0 / r^2 * sqrt(r / 2M)
 *
 *   "sonic"    - Transonic flow through sonic point (Fig 18)
 *                Critical point at r_c = 8M for gamma=5/3
 *
 * Domain: r ∈ [r_in, r_out] where r_in > 2M (outside horizon)
 * Particles flow inward from r_out toward the black hole.
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

void Solver::make_gr_bondi()
{
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real gamma = m_param->physics.gamma;
    const real c_speed = m_param->srgsph.c_speed;

    // Black hole mass
    real M = 1.0;
    if (m_sample_parameters.count("blackHoleMass")) {
        M = boost::any_cast<real>(m_sample_parameters["blackHoleMass"]);
    }

    // Test type: "geodesic" or "sonic"
    std::string test_type = "geodesic";
    if (m_sample_parameters.count("testType")) {
        test_type = boost::any_cast<std::string>(m_sample_parameters["testType"]);
    }

    // Domain (in units of M)
    real r_in = 2.5 * M;   // Inner edge (above horizon at 2M)
    real r_out = 18.0 * M; // Outer edge

    if (m_sample_parameters.count("rIn")) {
        r_in = boost::any_cast<real>(m_sample_parameters["rIn"]);
    }
    if (m_sample_parameters.count("rOut")) {
        r_out = boost::any_cast<real>(m_sample_parameters["rOut"]);
    }

    // Validate domain
    const real r_horizon = 2.0 * M;
    if (r_in <= r_horizon) {
        THROW_ERROR("Bondi accretion: r_in must be > 2M (horizon)");
    }

    WRITE_LOG << "GR-GSPH Bondi Accretion (" << test_type << "):";
    WRITE_LOG << "  Black hole mass M = " << M;
    WRITE_LOG << "  Horizon r_s = " << r_horizon;
    WRITE_LOG << "  Domain: r ∈ [" << r_in << ", " << r_out << "]";

    // Reference density at r_out
    const real n_0 = 1.0;

    // ============================================================================
    // Compute Bondi solution at each particle position
    // ============================================================================
    auto compute_geodesic = [M, gamma](real r, real n_0_ref, real r_ref) {
        // Geodesic (pressure-less) radial infall from Hawley, Smarr & Wilson (1984)
        // v^r (coordinate velocity) = -sqrt(2M/r) * (1 - 2M/r)
        // This is the velocity measured by a coordinate observer at infinity
        const real sqrt_2M_r = std::sqrt(2.0 * M / r);
        const real alpha_sq = 1.0 - 2.0 * M / r;  // lapse squared
        const real v_r = -sqrt_2M_r * alpha_sq;   // Coordinate velocity (inward)

        // Density profile from Liptai & Price (2019) / Hawley et al. (1984):
        // rho = rho_0 / r^2 * sqrt(r / 2M)
        const real n = n_0_ref / (r * r) * std::sqrt(r / (2.0 * M));

        // Thermal energy for geodesic flow (Hawley et al. 1984):
        // u = u_0 * (r^2 * sqrt(2M/r))^(1-gamma)
        // This follows from adiabatic compression during infall
        const real u_0 = 1.0e-9;  // Reference internal energy
        const real u = u_0 * std::pow(r * r * sqrt_2M_r, 1.0 - gamma);

        // Pressure from ideal gas: P = (gamma-1) * rho * u
        const real P = (gamma - 1.0) * n * u;

        return std::make_tuple(n, v_r, P, u);
    };

    // Sonic point flow constants (Liptai & Price 2019, equations 966-976)
    // Critical radius r_c = 8M for gamma = 5/3
    const real r_c = 8.0 * M;
    const real poly_n = 1.0 / (gamma - 1.0);  // Polytropic index n = 3/2 for gamma=5/3
    const real K_0 = 1.0;  // Entropy normalization

    // Critical point values
    const real U_c = std::sqrt(M / (2.0 * r_c));
    const real v_c = U_c / std::sqrt(1.0 - 3.0 * U_c * U_c);
    const real T_c = v_c * v_c * poly_n / ((1.0 + poly_n) * (1.0 - poly_n * v_c * v_c));

    // Constants C_1 and C_2
    const real C_1 = U_c * r_c * r_c * std::pow(T_c, poly_n);
    const real C_2 = std::pow(1.0 + (1.0 + poly_n) * T_c, 2.0) *
                     (1.0 - 2.0 * M / r_c + C_1 * C_1 / (std::pow(r_c, 4) * std::pow(T_c, 2.0 * poly_n)));

    WRITE_LOG << "  Sonic point: r_c = " << r_c << ", T_c = " << T_c;
    WRITE_LOG << "  Constants: C_1 = " << C_1 << ", C_2 = " << C_2;

    // Bisection solver for implicit T equation
    // Uses same approach as Python: search in [1e-10, 2*T_c] for inflow branch
    auto solve_T = [M, poly_n, C_1, C_2, T_c](real r) -> real {
        // Implicit equation: C_2 = (1 + (n+1)*T)^2 * (1 - 2M/r + C_1^2/(r^4 * T^(2n)))
        auto equation = [M, poly_n, C_1, C_2, r](real T) -> real {
            if (T <= 0.0) return 1.0e10;
            const real term1 = std::pow(1.0 + (poly_n + 1.0) * T, 2.0);
            const real term2 = 1.0 - 2.0 * M / r + C_1 * C_1 / (std::pow(r, 4) * std::pow(T, 2.0 * poly_n));
            return term1 * term2 - C_2;
        };

        // For inflow solution, search in narrow range [1e-10, 2*T_c]
        // Same as Python: brentq(equation, 1e-10, T_c * 2)
        real T_lo = 1.0e-10;
        real T_hi = T_c * 2.0;

        real f_lo = equation(T_lo);
        real f_hi = equation(T_hi);

        // If no sign change (no root in range), return T_c as fallback
        // This matches Python's except clause behavior
        if (f_lo * f_hi > 0) {
            return T_c;
        }

        // Bisection iteration
        for (int iter = 0; iter < 100; ++iter) {
            real T_mid = 0.5 * (T_lo + T_hi);
            real f_mid = equation(T_mid);

            if (std::abs(f_mid) < 1.0e-12 || (T_hi - T_lo) < 1.0e-14) {
                return T_mid;
            }

            if (f_lo * f_mid < 0) {
                T_hi = T_mid;
                f_hi = f_mid;
            } else {
                T_lo = T_mid;
                f_lo = f_mid;
            }
        }

        return 0.5 * (T_lo + T_hi);
    };

    auto compute_sonic = [M, gamma, poly_n, K_0, C_1, &solve_T](real r) {
        // Solve for temperature T(r)
        const real T = solve_T(r);

        // Density: rho = K_0 * T^n
        const real rho = K_0 * std::pow(T, poly_n);

        // Four-velocity radial component: U^r = C_1 / (r^2 * T^n)
        const real U_r = C_1 / (r * r * std::pow(T, poly_n));

        // Convert U^r to coordinate velocity v^r
        // U^0 = sqrt(1 - 2M/r + (U^r)^2) / (1 - 2M/r)
        const real alpha_sq = 1.0 - 2.0 * M / r;
        const real U_0 = std::sqrt(alpha_sq + U_r * U_r) / alpha_sq;
        const real v_r = -U_r / U_0;  // Negative for infall

        // Specific internal energy: u = n * T
        const real u = poly_n * T;

        // Pressure: P = rho * T = (gamma - 1) * rho * u
        const real P = rho * T;

        return std::make_tuple(rho, v_r, P, u);
    };

    // ============================================================================
    // Particle distribution
    // ============================================================================
    // Use equal-mass particles, spacing varies with density
    const real dr = (r_out - r_in) / N;

    auto& particles = m_sim->get_particles();
    particles.clear();
    particles.reserve(N + 20);  // Extra for ghost particles

    int particle_id = 0;

    for (int i = 0; i < N; ++i) {
        SPHParticle p;
        p.id = particle_id++;

        // Position
        const real r = r_in + (i + 0.5) * dr;
        p.pos[0] = r;
#if DIM >= 2
        p.pos[1] = 0.0;
#endif
#if DIM == 3
        p.pos[2] = 0.0;
#endif

        // Get Bondi solution at this radius
        real n, v_r, P, u_thermal;
        if (test_type == "geodesic") {
            std::tie(n, v_r, P, u_thermal) = compute_geodesic(r, n_0, r_out);
        } else {
            std::tie(n, v_r, P, u_thermal) = compute_sonic(r);
        }

        // Primitive variables
        p.dens = n;
        p.pres = P;
        p.vel[0] = v_r;  // Radial (inward) velocity
#if DIM >= 2
        p.vel[1] = 0.0;
#endif
#if DIM == 3
        p.vel[2] = 0.0;
#endif
        p.vel_t = 0.0;

        // Baryon number per particle (constant for equal-mass)
        // nu = n * dV where dV = dr for 1D
        p.nu = n * dr;
        p.mass = p.nu;

        // Specific internal energy (use analytic value for geodesic)
        p.ene = u_thermal;

        // Sound speed
        const real h = 1.0 + gamma / (gamma - 1.0) * P / n;
        p.sound = std::sqrt(gamma * P / (n * h));

        // Lorentz factor
        const real v2 = v_r * v_r;
        p.gamma_lor = 1.0 / std::sqrt(1.0 - v2);

        // Lab-frame density
        p.N = n * p.gamma_lor;

        // Smoothing length
        p.sml = 1.2 * dr;

        p.is_ghost = false;

        particles.push_back(p);
    }

    m_sim->set_particle_num(static_cast<int>(particles.size()));

    // ============================================================================
    // Ghost particles at boundaries
    // ============================================================================
    int ghost_layers = 6;
    if (m_sample_parameters.count("ghostLayers")) {
        ghost_layers = boost::any_cast<int>(m_sample_parameters["ghostLayers"]);
    }

    const int num_real = m_sim->get_particle_num();

    // Inner boundary ghosts (outflow - particles leaving domain)
    for (int i = 0; i < ghost_layers; ++i) {
        SPHParticle ghost;
        ghost.id = particle_id++;

        const real r = r_in - (i + 0.5) * dr;
        ghost.pos[0] = r;
#if DIM >= 2
        ghost.pos[1] = 0.0;
#endif
#if DIM == 3
        ghost.pos[2] = 0.0;
#endif

        // Use innermost particle values (outflow BC)
        real n, v_r, P, u_thermal;
        if (test_type == "geodesic") {
            std::tie(n, v_r, P, u_thermal) = compute_geodesic(r_in, n_0, r_out);
        } else {
            std::tie(n, v_r, P, u_thermal) = compute_sonic(r_in);
        }

        ghost.dens = n;
        ghost.pres = P;
        ghost.vel[0] = v_r;
#if DIM >= 2
        ghost.vel[1] = 0.0;
#endif
#if DIM == 3
        ghost.vel[2] = 0.0;
#endif
        ghost.vel_t = 0.0;
        ghost.nu = n * dr;
        ghost.mass = ghost.nu;
        ghost.ene = u_thermal;
        ghost.sound = particles[0].sound;
        ghost.gamma_lor = particles[0].gamma_lor;
        ghost.N = n * ghost.gamma_lor;
        ghost.sml = dr;
        ghost.is_ghost = true;

        particles.push_back(ghost);
    }

    // Outer boundary ghosts (inflow - fixed Bondi solution)
    for (int i = 0; i < ghost_layers; ++i) {
        SPHParticle ghost;
        ghost.id = particle_id++;

        const real r = r_out + (i + 0.5) * dr;
        ghost.pos[0] = r;
#if DIM >= 2
        ghost.pos[1] = 0.0;
#endif
#if DIM == 3
        ghost.pos[2] = 0.0;
#endif

        // Use exact Bondi solution at outer boundary (inflow BC)
        real n, v_r, P, u_thermal;
        if (test_type == "geodesic") {
            std::tie(n, v_r, P, u_thermal) = compute_geodesic(r, n_0, r_out);
        } else {
            std::tie(n, v_r, P, u_thermal) = compute_sonic(r);
        }

        ghost.dens = n;
        ghost.pres = P;
        ghost.vel[0] = v_r;
#if DIM >= 2
        ghost.vel[1] = 0.0;
#endif
#if DIM == 3
        ghost.vel[2] = 0.0;
#endif
        ghost.vel_t = 0.0;
        ghost.nu = n * dr;
        ghost.mass = ghost.nu;
        ghost.ene = u_thermal;
        ghost.sound = particles[num_real - 1].sound;

        const real v2 = v_r * v_r;
        ghost.gamma_lor = 1.0 / std::sqrt(1.0 - v2);
        ghost.N = n * ghost.gamma_lor;
        ghost.sml = dr;
        ghost.is_ghost = true;

        particles.push_back(ghost);
    }

    WRITE_LOG << "  Added " << (2 * ghost_layers) << " ghost particles at boundaries";
    WRITE_LOG << "Bondi accretion (" << test_type << ") initialized with " << particles.size() << " total particles";

    // Store parameters
    m_sample_parameters["r_in"] = r_in;
    m_sample_parameters["r_out"] = r_out;
    m_sample_parameters["bondi_type"] = test_type;
}

} // namespace sph
