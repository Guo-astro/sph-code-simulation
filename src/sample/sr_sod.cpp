// Special Relativistic Shock Tube Tests
// Implementation of Kitajima et al. (2025) arXiv:2510.18251v1 Section 3
//
// Test Types (set via "testType" parameter):
//   "sod"           - Sod shock tube (Section 3.1.1)
//   "blast_wave"    - Standard relativistic blast wave (Section 3.1.2)
//   "strong_blast"  - Strong relativistic blast wave (Section 3.1.3)
//
// Particle distribution: equal baryon number per particle (equal_nu)
// This follows Fig. 1 style from the paper where ν is constant across all particles.

#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"
#include "logger.hpp"
#include "srgsph/sr_primitive_recovery.hpp"
#include <cmath>
#include <string>

namespace sph
{

void Solver::make_sr_sod()
{
#if DIM != 1
    THROW_ERROR("SR Sod test requires DIM == 1");
#else

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real gamma = m_param->physics.gamma;
    const real c_speed = m_param->srgsph.c_speed;
    const real c2 = c_speed * c_speed;

    // ============================================================================
    // Test type selection (default: sod)
    // ============================================================================
    std::string test_type = "sod";
    if (m_sample_parameters.count("testType")) {
        test_type = boost::any_cast<std::string>(m_sample_parameters["testType"]);
    }

    // ============================================================================
    // Initial conditions based on test type
    // From Kitajima et al. (2025) Section 3.1
    // All tests use equal_nu mode (equal baryon number per particle)
    // ============================================================================
    real P_left, n_left, v_left;
    real P_right, n_right, v_right;
    int N_left, N_right;

    if (test_type == "sod") {
        // Section 3.1.1 Sod Problem (Eqs. 74-75)
        // Left:  (P, n, v) = (1.0, 1.0, 0)
        // Right: (P, n, v) = (0.1, 0.125, 0)
        P_left = 1.0;
        n_left = 1.0;
        v_left = 0.0;
        P_right = 0.1;
        n_right = 0.125;
        v_right = 0.0;

        // Figure 1: 3200 left + 400 right (equal baryon number per particle)
        N_left = (N > 0) ? static_cast<int>(N * 8.0 / 9.0) : 3200;
        N_right = (N > 0) ? static_cast<int>(N * 1.0 / 9.0) : 400;
        WRITE_LOG << "Test type: SOD (Section 3.1.1)";

    } else if (test_type == "blast_wave") {
        // Section 3.1.2 Standard Relativistic Blast Wave
        // Left:  (P, n, v) = (40/3, 10.0, 0)
        // Right: (P, n, v) = (10^-6, 1.0, 0)
        P_left = 40.0 / 3.0;
        n_left = 10.0;
        v_left = 0.0;
        P_right = 1.0e-6;
        n_right = 1.0;
        v_right = 0.0;

        // 5000 left + 500 right (equal baryon number per particle)
        N_left = (N > 0) ? static_cast<int>(N * 10.0 / 11.0) : 5000;
        N_right = (N > 0) ? static_cast<int>(N * 1.0 / 11.0) : 500;
        WRITE_LOG << "Test type: STANDARD BLAST WAVE (Section 3.1.2)";

    } else if (test_type == "strong_blast") {
        // Section 3.1.3 Strong Relativistic Blast Wave
        // Left:  (P, n, v) = (1000, 1.0, 0)
        // Right: (P, n, v) = (0.01, 1.0, 0)
        P_left = 1000.0;
        n_left = 1.0;
        v_left = 0.0;
        P_right = 0.01;
        n_right = 1.0;
        v_right = 0.0;

        // Equal density on both sides, so equal count gives equal ν
        N_left = N;
        N_right = N;
        WRITE_LOG << "Test type: STRONG BLAST WAVE (Section 3.1.3)";

    } else if (test_type == "ultra_relat") {
        // Section 3.2 Ultra-Relativistic Shock
        // Left:  (P, n, v) = (1.0, 1.0, v_left) - moving at relativistic speed
        // Right: (P, n, v) = (1.0, 1.0, 0) - at rest
        // Paper tests: v = 0.9, 0.99, 0.999, 0.9999, 0.999999999
        P_left = 1.0;
        n_left = 1.0;
        // Get v_left from parameters (default 0.9)
        if (m_sample_parameters.count("v_left")) {
            v_left = boost::any_cast<real>(m_sample_parameters["v_left"]);
        } else {
            v_left = 0.9;  // Default: 0.9c
        }
        P_right = 1.0;
        n_right = 1.0;
        v_right = 0.0;

        // 1000 particles on each side (Fig. 7)
        N_left = N;
        N_right = N;
        WRITE_LOG << "Test type: ULTRA-RELATIVISTIC SHOCK (Section 3.2), v_left=" << v_left << "c";

    } else {
        THROW_ERROR("Unknown test type: " + test_type);
    }

    // ============================================================================
    // Ghost particle configuration
    // Ghost particles extend beyond boundaries to provide proper neighbor support
    // ============================================================================
    bool use_ghost_particles = true;  // Enable ghost particles for better conservation
    if (m_sample_parameters.count("useGhostParticles")) {
        use_ghost_particles = boost::any_cast<bool>(m_sample_parameters["useGhostParticles"]);
    }
    
    // Number of ghost particle layers
    // Need enough layers so boundary particles have same neighbor count as interior
    // For ultra-relativistic test: particles travel v*t = 0.9*0.3 = 0.27
    // Need ghost particles to extend at least that far for proper inflow
    // At dx=0.0005, need ~540 layers. Using 600 to be safe.
    int ghost_layers = 6;  // Default: 6 layers for standard tests
    if (m_sample_parameters.count("ghostLayers")) {
        ghost_layers = boost::any_cast<int>(m_sample_parameters["ghostLayers"]);
    }
    
    // For ultra-relativistic test, we need MANY more ghost particles
    // to simulate continuous inflow over the simulation time
    if (std::abs(v_left) > 0.1) {
        // Calculate layers needed: v_left * end_time / dx
        // Assuming end_time ~ 0.3 and N ~ 1000, dx ~ 0.5/1000 = 0.0005
        // For v=0.9, need 0.9*0.3/0.0005 = 540 layers
        const real end_time = 0.3;  // Approximate simulation time
        const real dx_approx = 0.5 / N;  // Approximate particle spacing
        const int needed_layers = static_cast<int>(std::ceil(v_left * end_time / dx_approx)) + 10;
        ghost_layers = std::max(ghost_layers, needed_layers);
        WRITE_LOG << "Ultra-relativistic: using " << ghost_layers << " ghost layers for inflow";
    }

    const int N_ghost_left = use_ghost_particles ? ghost_layers : 0;
    const int N_ghost_right = use_ghost_particles ? ghost_layers : 0;
    const int num = N_left + N_right + N_ghost_left + N_ghost_right;

    // Domain boundaries from config (rangeMin/rangeMax), discontinuity at x=0
    const real x_left_start = m_param->periodic.range_min[0];
    const real x_left_end = 0.0;
    const real x_right_start = 0.0;
    const real x_right_end = m_param->periodic.range_max[0];

    const real dx_left = (x_left_end - x_left_start) / N_left;
    const real dx_right = (x_right_end - x_right_start) / N_right;

    // Lorentz factors for initial states
    const real gamma_left = 1.0 / std::sqrt(1.0 - v_left * v_left / c2);
    const real gamma_right = 1.0 / std::sqrt(1.0 - v_right * v_right / c2);

    // Baryon number per particle (ν)
    // In the lab frame where we set up particles:
    //   - Particles are placed at spacing dx (lab-frame spacing)
    //   - For moving fluid, lab-frame density is N = γ * n
    //   - Baryon number is ν = N * V_lab = γ * n * dx
    // This ensures that when the kernel sum gives N_lab, the rest-frame
    // density recovered is n = N_lab / γ = correct initial value
    const real nu_left = gamma_left * n_left * dx_left;
    const real nu_right = gamma_right * n_right * dx_right;

    std::vector<SPHParticle> p(num);

    // ============================================================================
    // Initialize left state particles
    // ============================================================================
    for (int i = 0; i < N_left; ++i) {
        auto& p_i = p[i];
        p_i.id = i;
        p_i.pos[0] = x_left_start + (i + 0.5) * dx_left;

        // Baryon number per particle
        p_i.nu = nu_left;
        p_i.mass = nu_left;

        vec_t vel;
        vel[0] = v_left;

        // Compute Lorentz factor
        const real v2 = vel[0] * vel[0];
        p_i.gamma_lor = 1.0 / std::sqrt(1.0 - v2 / c2);

        // Thermodynamic quantities
        const real u_init = P_left / ((gamma - 1.0) * n_left);
        const real H_init = 1.0 + u_init / c2 + P_left / (n_left * c2);
        p_i.sound = std::sqrt((gamma - 1.0) * (H_init - 1.0) / H_init) * c_speed;
        p_i.enthalpy = H_init;

        // Conserved variables
        const real N_conserved = p_i.gamma_lor * n_left;
        const vec_t S_conserved = vel * (p_i.gamma_lor * H_init);
        const real X = gamma / (gamma - 1.0);
        const real e_conserved = (H_init * (X * p_i.gamma_lor * p_i.gamma_lor - 1.0) + 1.0) / (X * p_i.gamma_lor);

        p_i.N = N_conserved;
        p_i.S = S_conserved;
        p_i.e = e_conserved;
        p_i.vel = vel;
        p_i.ene = u_init;
        p_i.pres = P_left;
        p_i.dens = N_conserved;

        // Initialize derivatives
        p_i.dS = vec_t(0.0);
        p_i.de = 0.0;
        p_i.dS_old = vec_t(0.0);
        p_i.de_old = 0.0;

        // Smoothing length initial guess
        p_i.sml = dx_left;
    }

    // ============================================================================
    // Initialize right state particles
    // ============================================================================
    for (int i = 0; i < N_right; ++i) {
        auto& p_i = p[N_left + i];
        p_i.id = N_left + i;
        p_i.pos[0] = x_right_start + (i + 0.5) * dx_right;

        // Baryon number per particle
        p_i.nu = nu_right;
        p_i.mass = nu_right;

        vec_t vel;
        vel[0] = v_right;

        // Compute Lorentz factor
        const real v2 = vel[0] * vel[0];
        p_i.gamma_lor = 1.0 / std::sqrt(1.0 - v2 / c2);

        // Thermodynamic quantities
        const real u_init = P_right / ((gamma - 1.0) * n_right);
        const real H_init = 1.0 + u_init / c2 + P_right / (n_right * c2);
        p_i.sound = std::sqrt((gamma - 1.0) * (H_init - 1.0) / H_init) * c_speed;
        p_i.enthalpy = H_init;

        // Conserved variables
        const real N_conserved = p_i.gamma_lor * n_right;
        const vec_t S_conserved = vel * (p_i.gamma_lor * H_init);
        const real X = gamma / (gamma - 1.0);
        const real e_conserved = (H_init * (X * p_i.gamma_lor * p_i.gamma_lor - 1.0) + 1.0) / (X * p_i.gamma_lor);

        p_i.N = N_conserved;
        p_i.S = S_conserved;
        p_i.e = e_conserved;
        p_i.vel = vel;
        p_i.ene = u_init;
        p_i.pres = P_right;
        p_i.dens = N_conserved;

        // Initialize derivatives
        p_i.dS = vec_t(0.0);
        p_i.de = 0.0;
        p_i.dS_old = vec_t(0.0);
        p_i.de_old = 0.0;

        // Smoothing length initial guess
        p_i.sml = dx_right;
    }

    // ============================================================================
    // Initialize ghost particles (left boundary - extend left state)
    // Ghost particles mirror the boundary conditions for proper neighbor support
    // ============================================================================
    if (use_ghost_particles) {
        int ghost_id = N_left + N_right;
        
        // Left ghost particles (extend left state beyond x = -0.5)
        // For INFLOW boundary (ultra-relativistic test): use same velocity as interior
        // For REFLECTING boundary (standard Sod): use reflected velocity
        // Since v_left=0 for Sod test, -v_left = 0 anyway
        // For ultra-relativistic test with v_left=0.9, we want inflow (same velocity)
        const bool use_inflow = (std::abs(v_left) > 0.1);  // Inflow if v_left is significant
        
        for (int i = 0; i < N_ghost_left; ++i) {
            auto& p_i = p[ghost_id];
            p_i.id = ghost_id;
            // Place ghost particles beyond left boundary, mirroring spacing
            p_i.pos[0] = x_left_start - (i + 0.5) * dx_left;
            p_i.is_ghost = true;  // Mark as ghost particle

            // Copy left state properties
            p_i.nu = nu_left;
            p_i.mass = nu_left;

            vec_t vel;
            // INFLOW: same velocity (for ultra-relativistic with moving left state)
            // REFLECTING: negate velocity (for Sod with v_left=0, no difference)
            vel[0] = use_inflow ? v_left : -v_left;

            const real v2 = vel[0] * vel[0];
            p_i.gamma_lor = 1.0 / std::sqrt(1.0 - v2 / c2);

            const real u_init = P_left / ((gamma - 1.0) * n_left);
            const real H_init = 1.0 + u_init / c2 + P_left / (n_left * c2);
            p_i.sound = std::sqrt((gamma - 1.0) * (H_init - 1.0) / H_init) * c_speed;
            p_i.enthalpy = H_init;

            const real N_conserved = p_i.gamma_lor * n_left;
            const vec_t S_conserved = vel * (p_i.gamma_lor * H_init);
            const real X = gamma / (gamma - 1.0);
            const real e_conserved = (H_init * (X * p_i.gamma_lor * p_i.gamma_lor - 1.0) + 1.0) / (X * p_i.gamma_lor);

            p_i.N = N_conserved;
            p_i.S = S_conserved;
            p_i.e = e_conserved;
            p_i.vel = vel;
            p_i.ene = u_init;
            p_i.pres = P_left;
            p_i.dens = N_conserved;

            p_i.dS = vec_t(0.0);
            p_i.de = 0.0;
            p_i.dS_old = vec_t(0.0);
            p_i.de_old = 0.0;
            p_i.sml = dx_left;

            ghost_id++;
        }

        // Right ghost particles (extend right state beyond x = 0.5)
        for (int i = 0; i < N_ghost_right; ++i) {
            auto& p_i = p[ghost_id];
            p_i.id = ghost_id;
            // Place ghost particles beyond right boundary
            p_i.pos[0] = x_right_end + (i + 0.5) * dx_right;
            p_i.is_ghost = true;  // Mark as ghost particle

            // Copy right state properties
            p_i.nu = nu_right;
            p_i.mass = nu_right;

            vec_t vel;
            vel[0] = -v_right;  // Reflect velocity for wall condition

            const real v2 = vel[0] * vel[0];
            p_i.gamma_lor = 1.0 / std::sqrt(1.0 - v2 / c2);

            const real u_init = P_right / ((gamma - 1.0) * n_right);
            const real H_init = 1.0 + u_init / c2 + P_right / (n_right * c2);
            p_i.sound = std::sqrt((gamma - 1.0) * (H_init - 1.0) / H_init) * c_speed;
            p_i.enthalpy = H_init;

            const real N_conserved = p_i.gamma_lor * n_right;
            const vec_t S_conserved = vel * (p_i.gamma_lor * H_init);
            const real X = gamma / (gamma - 1.0);
            const real e_conserved = (H_init * (X * p_i.gamma_lor * p_i.gamma_lor - 1.0) + 1.0) / (X * p_i.gamma_lor);

            p_i.N = N_conserved;
            p_i.S = S_conserved;
            p_i.e = e_conserved;
            p_i.vel = vel;
            p_i.ene = u_init;
            p_i.pres = P_right;
            p_i.dens = N_conserved;

            p_i.dS = vec_t(0.0);
            p_i.de = 0.0;
            p_i.dS_old = vec_t(0.0);
            p_i.de_old = 0.0;
            p_i.sml = dx_right;

            ghost_id++;
        }
    }

    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());

    WRITE_LOG << "SR shock tube initialized (Kitajima et al. 2025):";
    WRITE_LOG << "  Left:  " << N_left << " particles, P=" << P_left << ", n=" << n_left << ", v=" << v_left;
    WRITE_LOG << "  Right: " << N_right << " particles, P=" << P_right << ", n=" << n_right << ", v=" << v_right;
    if (use_ghost_particles) {
        WRITE_LOG << "  Ghost: " << (N_ghost_left + N_ghost_right) << " particles (" << N_ghost_left << " left, " << N_ghost_right << " right)";
    }
    WRITE_LOG << "  Baryon numbers: nu_L=" << nu_left << ", nu_R=" << nu_right;
    WRITE_LOG << "  dx_left=" << dx_left << ", dx_right=" << dx_right;
#endif
}

}
