/**
 * @file bns_cocoon_2d.cpp
 * @brief 2D axisymmetric cocoon shock breakout simulation for BNS mergers
 * 
 * This sample creates a 2D simulation of a jet-driven cocoon propagating
 * through BNS merger ejecta with angular dependence. Physics based on:
 * - Gutiérrez et al. (2024) arXiv:2408.15973v3 for cocoon breakout
 * - Radice et al. (2018) for NR ejecta profiles
 * - Shibata & Hotokezaka (2019) for GRB 170817A constraints
 * 
 * Two ejecta profile options:
 * 1. Sharp cutoff: ρ → 0 at v_max
 * 2. Extended tail: ρ ∝ (Γβ)^(-α) beyond v_break (preferred for GRB 170817A)
 * 
 * Two initialization modes:
 * 1. Direct cocoon: Hot bubble at center (fast, simplified)
 * 2. Jet injection: Continuous jet launching (realistic, slower)
 * 
 * Physics: Special Relativistic Godunov SPH (SRGSPH)
 */

#include "solver.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "openmp.hpp"

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <random>

namespace sph
{

namespace {

// Physical constants in code units (c = 1)
constexpr real C_LIGHT = 1.0;
constexpr real PI = 3.14159265358979323846;

/**
 * @brief Compute Lorentz factor from velocity
 */
inline real lorentz_factor(real v) {
    return 1.0 / std::sqrt(1.0 - v * v);
}

/**
 * @brief Compute Γβ = v / sqrt(1 - v²)
 */
inline real gamma_beta(real v) {
    return v / std::sqrt(1.0 - v * v);
}

/**
 * @brief Ejecta profile types following Gutiérrez+2024
 */
enum class EjectaProfileType {
    SHARP_CUTOFF,    // ρ → 0 at v_max
    EXTENDED_TAIL    // ρ ∝ (Γβ)^(-α) beyond v_break
};

/**
 * @brief 2D axisymmetric ejecta profile with angular dependence
 * 
 * Based on NR simulations from Radice+2018 and fits from Gutiérrez+2024.
 * 
 * Radial profile:
 *   For v < v_break: ρ ∝ r^(-n) where n = β + 2 (homologous expansion)
 *   For v > v_break (extended tail): ρ ∝ (Γβ)^(-α)
 * 
 * Angular distribution:
 *   Polar region (θ < θ_polar): Lower density (jet corridor)
 *   Equatorial (θ ~ 90°): Maximum density
 *   Distribution: ρ(θ) ∝ sin^(angular_power)(θ)
 */
struct EjectaProfile2D {
    // Ejecta parameters from NR simulations
    real M_ej;              // Total ejecta mass (code units, ~M_sun)
    real v_min;             // Minimum ejecta velocity (c)
    real v_max;             // Maximum ejecta velocity (c) - from NR: 0.6-0.8c
    real v_break;           // Break velocity for extended tail (c)
    real beta;              // Radial power-law index in inner region
    real alpha_tail;        // Power-law index for extended tail: ρ ∝ (Γβ)^(-α)
    real t0;                // Reference time (simulation start time)
    real gamma_eos;         // Adiabatic index (4/3 for radiation dominated)
    real rho_floor;         // Minimum density (ambient medium)
    
    // Angular distribution parameters
    real theta_polar;       // Polar opening angle (radians) - lower density region
    real polar_density_factor;  // Density reduction in polar region
    real angular_power;     // Power for angular distribution sin^n(θ)
    
    // Profile type
    EjectaProfileType profile_type;
    
    // Derived quantities
    real r_min, r_max, r_break;
    real rho_norm;
    real rho_at_break;      // Density at break radius for tail matching
    real gb_at_break;       // Γβ at break velocity
    
    void compute_derived() {
        r_min = v_min * t0;
        r_max = v_max * t0;
        r_break = v_break * t0;
        gb_at_break = gamma_beta(v_break);
        
        // Normalization: integrate ρ(r,θ) over volume = M_ej
        real n = beta + 2.0;
        real radial_integral;
        
        // Integrate r² ρ(r) dr from r_min to r_max
        if (profile_type == EjectaProfileType::SHARP_CUTOFF) {
            if (std::abs(n - 3.0) < 1e-10) {
                radial_integral = std::log(r_max / r_min);
            } else {
                radial_integral = (std::pow(r_max, 3.0 - n) - std::pow(r_min, 3.0 - n)) / (3.0 - n);
            }
        } else {
            // Extended tail: integrate in two parts
            if (std::abs(n - 3.0) < 1e-10) {
                radial_integral = std::log(r_break / r_min);
            } else {
                radial_integral = (std::pow(r_break, 3.0 - n) - std::pow(r_min, 3.0 - n)) / (3.0 - n);
            }
            // Add tail contribution (approximate)
            radial_integral += 0.2 * std::pow(r_break, 3.0);  // Rough estimate
        }
        
        // Angular integral: ∫ sin^(n+1)(θ) dθ from 0 to π
        // For sin²(θ): integral = 4/3
        // For sin⁴(θ): integral = 16/15
        real angular_integral;
        if (angular_power <= 0.1) {
            angular_integral = 2.0;  // Isotropic
        } else if (std::abs(angular_power - 2.0) < 0.1) {
            angular_integral = 4.0 / 3.0;
        } else if (std::abs(angular_power - 4.0) < 0.1) {
            angular_integral = 16.0 / 15.0;
        } else {
            // Approximate for other values
            angular_integral = 2.0 / (1.0 + angular_power / 2.0);
        }
        
        // Account for polar density reduction
        real polar_fraction = theta_polar / PI;
        real effective_angular = angular_integral * (1.0 - polar_fraction * (1.0 - polar_density_factor));
        
        rho_norm = M_ej / (4.0 * PI * radial_integral * effective_angular);
        
        // Compute density at break radius for tail matching
        rho_at_break = rho_norm * std::pow(r_break / r_min, -n);
    }
    
    /**
     * @brief Angular density factor
     * 
     * From Radice+2018 and Gutiérrez+2024:
     * - Polar region (θ < θ_polar): reduced density for jet corridor
     * - Equatorial: maximum density (dynamical + disk wind ejecta)
     */
    real angular_factor(real theta) const {
        // Normalize theta to [0, π]
        theta = std::fmod(std::abs(theta), PI);
        if (theta > PI) theta = 2.0 * PI - theta;
        
        // Base angular distribution: peaks at equator
        real sin_theta = std::sin(theta);
        real ang_dist = std::pow(std::abs(sin_theta), angular_power);
        
        // Polar density reduction (for jet corridor)
        if (theta < theta_polar || theta > PI - theta_polar) {
            ang_dist *= polar_density_factor;
        }
        
        return std::max(ang_dist, 0.01);  // Minimum 1% to avoid zero density
    }
    
    /**
     * @brief Density at given radius and angle
     */
    real density(real r, real theta) const {
        if (r < r_min) return rho_floor;
        
        real rho;
        if (r <= r_break || profile_type == EjectaProfileType::SHARP_CUTOFF) {
            if (r > r_max) return rho_floor;
            // Power-law in inner region
            real n = beta + 2.0;
            rho = rho_norm * std::pow(r / r_min, -n);
        } else {
            // Extended tail: ρ ∝ (Γβ)^(-α)
            real v = r / t0;
            if (v >= C_LIGHT) return rho_floor;  // Cannot exceed c
            real gb = gamma_beta(v);
            rho = rho_at_break * std::pow(gb / gb_at_break, -alpha_tail);
        }
        
        return std::max(rho_floor, rho * angular_factor(theta));
    }
    
    /**
     * @brief Velocity at given radius (homologous expansion)
     */
    real velocity(real r) const {
        real v = r / t0;
        return std::min(v, 0.99 * C_LIGHT);  // Cap below c
    }
    
    /**
     * @brief Pressure from density (cold EOS)
     */
    real pressure(real rho) const {
        // Cold ejecta: P ∝ ρ^γ with small K for subsonic expansion
        real cs_target = 0.01 * C_LIGHT;  // Sound speed << expansion velocity
        real K = cs_target * cs_target / (gamma_eos * std::pow(std::max(rho, rho_floor), gamma_eos - 1.0));
        return K * std::pow(std::max(rho, rho_floor), gamma_eos);
    }
    
    /**
     * @brief Internal energy from pressure and density
     */
    real internal_energy(real pres, real rho) const {
        return pres / ((gamma_eos - 1.0) * std::max(rho, rho_floor));
    }
};

/**
 * @brief Cocoon/Jet parameters following Gutiérrez+2024
 */
struct CocoonParams {
    // Injection mode
    bool jet_injection;         // true: continuous jet, false: initial hot bubble
    
    // Cocoon geometry
    real r_cocoon;              // Cocoon outer radius
    real theta_jet;             // Jet half-opening angle (radians) - paper: 8-16°
    
    // Cocoon thermodynamics
    real E_cocoon;              // Total cocoon energy
    real Gamma_cocoon;          // Lorentz factor at injection
    real rho_cocoon;            // Cocoon density (relative to ejecta)
    
    // Jet injection parameters (if jet_injection = true)
    real L_jet;                 // Jet luminosity
    real t_jet_start;           // Jet launch delay (t_del in paper)
    real t_jet_duration;        // Jet active duration
    real Gamma_jet;             // Jet Lorentz factor (~200 for SGRB)
};

}  // anonymous namespace

void Solver::make_bns_cocoon_2d()
{
#if DIM != 2
    THROW_ERROR("BNS Cocoon 2D requires DIM == 2");
#else
    std::cout << "=== Initializing 2D BNS Cocoon Shock Breakout ===" << std::endl;
    std::cout << "Physics based on Gutiérrez et al. (2024) arXiv:2408.15973v3" << std::endl;
    
    // =========================================================================
    // Default parameters from paper Table 1 (SLy_M145-125 model)
    // =========================================================================
    
    // Ejecta parameters (NR-derived defaults)
    real M_ej = 0.01;           // Total ejecta mass (code units ~ M_sun)
    real v_min = 0.05;          // Minimum velocity (c)
    real v_max = 0.7;           // Maximum velocity (c) - paper: 0.58-0.78
    real v_break = 0.5;         // Break velocity for extended tail (c)
    real beta = 0.5;            // Radial power-law index
    real alpha_tail = 6.0;      // Extended tail power-law index
    real t0 = 1.0;              // Reference time (code units)
    real gamma_eos = 4.0/3.0;   // Radiation-dominated EOS
    real rho_floor = 1e-12;     // Ambient medium density
    
    // Profile type: 0 = sharp cutoff, 1 = extended tail
    int profile_type_int = 1;   // Extended tail preferred for GRB 170817A
    
    // Angular distribution (from NR fits)
    real theta_polar = 15.0 * PI / 180.0;   // 15° polar opening
    real polar_density_factor = 0.1;        // 10% density in polar region
    real angular_power = 2.0;               // sin²(θ) distribution
    
    // Cocoon/jet parameters (from paper)
    bool jet_injection = false;             // Start with hot bubble mode
    real E_cocoon = 1e-4;                   // Cocoon energy (code units)
    real r_cocoon_frac = 0.3;               // Cocoon radius as fraction of r_min
    real theta_jet = 12.0 * PI / 180.0;     // 12° half-opening (paper: 8-16°)
    real Gamma_cocoon = 2.0;                // Initial Lorentz factor
    
    // Resolution parameters
    int n_radial = 200;         // Radial resolution
    int n_angular = 160;        // Angular resolution
    
    // =========================================================================
    // Read parameters from config
    // =========================================================================
    
    auto read_real = [this](const std::string& key, real& val) {
        if (m_sample_parameters.count(key))
            val = boost::any_cast<real>(m_sample_parameters[key]);
    };
    auto read_int = [this](const std::string& key, int& val) {
        if (m_sample_parameters.count(key))
            val = boost::any_cast<int>(m_sample_parameters[key]);
    };
    auto read_bool = [this](const std::string& key, bool& val) {
        if (m_sample_parameters.count(key))
            val = boost::any_cast<bool>(m_sample_parameters[key]);
    };
    
    // Ejecta
    read_real("M_ej", M_ej);
    read_real("v_min", v_min);
    read_real("v_max", v_max);
    read_real("v_break", v_break);
    read_real("beta", beta);
    read_real("alpha_tail", alpha_tail);
    read_real("t0", t0);
    read_real("gamma", gamma_eos);
    read_real("rho_floor", rho_floor);
    read_int("profile_type", profile_type_int);
    
    // Angular
    read_real("theta_polar", theta_polar);
    read_real("polar_density_factor", polar_density_factor);
    read_real("angular_power", angular_power);
    
    // Cocoon
    read_bool("jet_injection", jet_injection);
    read_real("E_cocoon", E_cocoon);
    read_real("r_cocoon_frac", r_cocoon_frac);
    read_real("theta_jet", theta_jet);
    read_real("Gamma_cocoon", Gamma_cocoon);
    
    // Resolution
    read_int("n_radial", n_radial);
    read_int("n_angular", n_angular);
    
    // =========================================================================
    // Setup ejecta profile
    // =========================================================================
    
    EjectaProfile2D ejecta;
    ejecta.M_ej = M_ej;
    ejecta.v_min = v_min;
    ejecta.v_max = v_max;
    ejecta.v_break = v_break;
    ejecta.beta = beta;
    ejecta.alpha_tail = alpha_tail;
    ejecta.t0 = t0;
    ejecta.gamma_eos = gamma_eos;
    ejecta.rho_floor = rho_floor;
    ejecta.theta_polar = theta_polar;
    ejecta.polar_density_factor = polar_density_factor;
    ejecta.angular_power = angular_power;
    ejecta.profile_type = (profile_type_int == 0) ? 
        EjectaProfileType::SHARP_CUTOFF : EjectaProfileType::EXTENDED_TAIL;
    ejecta.compute_derived();
    
    real r_min = ejecta.r_min;
    real r_max = ejecta.r_max;
    real r_cocoon = r_cocoon_frac * r_min;
    
    // =========================================================================
    // Print configuration
    // =========================================================================
    
    std::cout << "\n--- Ejecta Configuration ---" << std::endl;
    std::cout << "  M_ej = " << M_ej << " (code units ~ M_sun)" << std::endl;
    std::cout << "  v_min = " << v_min << " c, v_max = " << v_max << " c" << std::endl;
    std::cout << "  r_min = " << r_min << ", r_max = " << r_max << std::endl;
    std::cout << "  Profile: " << (profile_type_int == 0 ? "Sharp cutoff" : "Extended tail") << std::endl;
    if (profile_type_int != 0) {
        std::cout << "  v_break = " << v_break << " c, alpha_tail = " << alpha_tail << std::endl;
    }
    std::cout << "  Angular: sin^" << angular_power << "(θ)" << std::endl;
    std::cout << "  Polar opening: " << theta_polar * 180.0 / PI << "°, factor = " << polar_density_factor << std::endl;
    
    std::cout << "\n--- Cocoon Configuration ---" << std::endl;
    std::cout << "  Mode: " << (jet_injection ? "Jet injection" : "Hot bubble") << std::endl;
    std::cout << "  E_cocoon = " << E_cocoon << std::endl;
    std::cout << "  r_cocoon = " << r_cocoon << " (" << r_cocoon_frac * 100 << "% of r_min)" << std::endl;
    std::cout << "  θ_jet = " << theta_jet * 180.0 / PI << "°" << std::endl;
    std::cout << "  Γ_cocoon = " << Gamma_cocoon << std::endl;
    
    std::cout << "\n--- Resolution ---" << std::endl;
    std::cout << "  n_radial = " << n_radial << ", n_angular = " << n_angular << std::endl;
    std::cout << "  Estimated particles: ~" << n_radial * n_angular * 2 << std::endl;
    
    // =========================================================================
    // Create particles in 2D (r-z plane, axisymmetric)
    // =========================================================================
    
    auto& particles = m_sim->get_particles();
    particles.clear();
    
    // Domain extends beyond r_max for ambient medium
    real domain_max = 1.2 * r_max;
    
    // Spacing
    real dr = domain_max / n_radial;
    real dtheta = PI / n_angular;
    
    // Estimate particle mass (total mass / number of particles)
    real particle_mass = M_ej / (n_radial * n_angular);
    real initial_sml = 2.0 * dr;
    
    std::cout << "\n--- Creating particles ---" << std::endl;
    std::cout << "  Domain: [0, " << domain_max << "]" << std::endl;
    std::cout << "  dr = " << dr << ", dθ = " << dtheta * 180.0 / PI << "°" << std::endl;
    
    int id = 0;
    int n_ejecta = 0, n_cocoon = 0, n_ambient = 0;
    
    // Create particles in polar coordinates, map to Cartesian
    for (int ir = 0; ir < n_radial; ++ir) {
        real r = (ir + 0.5) * dr;
        
        for (int itheta = 0; itheta < n_angular; ++itheta) {
            real theta = (itheta + 0.5) * dtheta;  // 0 to π
            
            // Skip very small angles to avoid singularity at pole
            if (theta < 0.01 || theta > PI - 0.01) continue;
            
            // Cartesian coordinates (x = r sin θ, y = r cos θ)
            // y is the polar (jet) axis
            real x = r * std::sin(theta);
            real y = r * std::cos(theta);
            
            SPHParticle p;
            p.id = id++;
            p.pos[0] = x;
            p.pos[1] = y;
            p.vel[0] = 0.0;
            p.vel[1] = 0.0;
            p.acc[0] = 0.0;
            p.acc[1] = 0.0;
            p.mass = particle_mass;
            p.sml = initial_sml;
            
            // Determine region
            bool in_cocoon = (r < r_cocoon) && (theta < theta_jet || theta > PI - theta_jet);
            bool in_ejecta = (r >= r_min && r <= r_max) || 
                             (r < r_min && !in_cocoon);
            
            if (in_cocoon) {
                // ==== Cocoon region: hot, relativistic ====
                n_cocoon++;
                
                // Cocoon volume (bipolar cones)
                real V_cocoon = (2.0/3.0) * PI * std::pow(r_cocoon, 3) * (1.0 - std::cos(theta_jet)) * 2.0;
                
                // Hot, low-density cocoon
                p.dens = 0.01 * ejecta.density(r_min, PI/2.0);  // 1% of inner ejecta density
                p.pres = E_cocoon / V_cocoon;
                p.ene = ejecta.internal_energy(p.pres, p.dens);
                
                // Relativistic outward velocity
                real v_cocoon = std::sqrt(1.0 - 1.0/(Gamma_cocoon * Gamma_cocoon));
                p.vel[0] = v_cocoon * std::sin(theta);
                p.vel[1] = v_cocoon * std::cos(theta);
                
            } else if (in_ejecta && r >= r_min) {
                // ==== Ejecta region: homologous expansion ====
                n_ejecta++;
                
                p.dens = ejecta.density(r, theta);
                p.pres = ejecta.pressure(p.dens);
                p.ene = ejecta.internal_energy(p.pres, p.dens);
                
                // Homologous velocity: v = r/t
                real v_rad = ejecta.velocity(r);
                p.vel[0] = v_rad * std::sin(theta);
                p.vel[1] = v_rad * std::cos(theta);
                
            } else if (r < r_min && !in_cocoon) {
                // ==== Inner transition region ====
                n_ejecta++;
                
                // Smooth transition from center to ejecta
                real frac = r / r_min;
                real rho_ref = ejecta.density(r_min, theta);
                p.dens = frac * rho_ref + (1.0 - frac) * 0.1 * rho_ref;
                p.pres = ejecta.pressure(p.dens);
                p.ene = ejecta.internal_energy(p.pres, p.dens);
                
                real v_rad = frac * ejecta.velocity(r_min);
                p.vel[0] = v_rad * std::sin(theta);
                p.vel[1] = v_rad * std::cos(theta);
                
            } else {
                // ==== Ambient medium ====
                n_ambient++;
                
                p.dens = rho_floor;
                p.pres = ejecta.pressure(p.dens);
                p.ene = ejecta.internal_energy(p.pres, p.dens);
                p.vel[0] = 0.0;
                p.vel[1] = 0.0;
            }
            
            // SR-GSPH requires baryon number
            p.nu = p.mass;
            
            // Initialize other SPH fields
            p.gradh = 1.0;
            p.alpha = m_param->av.alpha;
            p.balsara = 1.0;
            p.phi = 0.0;
            p.sound = std::sqrt(gamma_eos * p.pres / std::max(p.dens, rho_floor));
            
            particles.push_back(p);
        }
    }
    
    // Mirror to negative y hemisphere (except equatorial particles)
    int n_upper = particles.size();
    for (int i = 0; i < n_upper; ++i) {
        if (std::abs(particles[i].pos[1]) > 1e-10 * domain_max) {
            SPHParticle p = particles[i];
            p.id = id++;
            p.pos[1] = -p.pos[1];
            p.vel[1] = -p.vel[1];
            particles.push_back(p);
        }
    }
    
    // =========================================================================
    // Summary
    // =========================================================================
    
    std::cout << "\n--- Particle Summary ---" << std::endl;
    std::cout << "  Total particles: " << particles.size() << std::endl;
    std::cout << "  Cocoon: " << n_cocoon * 2 << std::endl;
    std::cout << "  Ejecta: " << n_ejecta * 2 << std::endl;
    std::cout << "  Ambient: " << n_ambient * 2 << std::endl;
    
    // Set particle count in simulation
    m_sim->set_particle_num(particles.size());
    
    // Resolution quality estimate
    real h_over_r_min = initial_sml / r_min;
    real particles_per_jet_angle = theta_jet / dtheta;
    std::cout << "\n--- Resolution Quality ---" << std::endl;
    std::cout << "  h/r_min = " << h_over_r_min << " (should be < 0.1 for good resolution)" << std::endl;
    std::cout << "  Particles across jet angle: " << particles_per_jet_angle 
              << " (should be > 3)" << std::endl;
    
    if (h_over_r_min > 0.2 || particles_per_jet_angle < 3) {
        std::cout << "  WARNING: Resolution may be too low for accurate breakout physics!" << std::endl;
        std::cout << "           Consider increasing n_radial and n_angular." << std::endl;
    }
    
    std::cout << "\n=== BNS Cocoon 2D initialization complete ===" << std::endl;
    
#endif  // DIM == 2
}

}  // namespace sph
