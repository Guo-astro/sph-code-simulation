/**
 * @file bns_cocoon_1d.cpp
 * @brief 1D spherical cocoon shock breakout simulation for BNS mergers
 * 
 * This sample creates a 1D spherically symmetric simulation of a cocoon
 * propagating through BNS merger ejecta. Physics based on:
 * - Gutiérrez et al. (2024) arXiv:2408.15973v3 for cocoon breakout
 * - Radice et al. (2018) for NR ejecta profiles
 * - Nakar & Piran (2017) for analytic cocoon models
 * 
 * Two ejecta profile options:
 * 1. Sharp cutoff: ρ → 0 at v_max
 * 2. Extended tail: ρ ∝ (Γβ)^(-α) beyond v_break (preferred for GRB 170817A)
 * 
 * Use case: Rapid testing and parameter exploration before 2D/3D runs.
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
 * @brief 1D spherical ejecta profile for BNS merger
 * 
 * Based on NR simulations from Radice+2018 and fits from Gutiérrez+2024.
 * 
 * Radial profile:
 *   For v < v_break: ρ ∝ r^(-n) where n = β + 2 (homologous expansion)
 *   For v > v_break (extended tail): ρ ∝ (Γβ)^(-α)
 */
struct EjectaProfile1D {
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
        
        // Normalization: 4π ∫ r² ρ(r) dr = M_ej
        real n = beta + 2.0;
        real radial_integral;
        
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
            radial_integral += 0.2 * std::pow(r_break, 3.0);
        }
        
        rho_norm = M_ej / (4.0 * PI * radial_integral);
        
        // Compute density at break radius for tail matching
        rho_at_break = rho_norm * std::pow(r_break / r_min, -n);
    }
    
    /**
     * @brief Density at given radius
     */
    real density(real r) const {
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
            if (v >= C_LIGHT) return rho_floor;
            real gb = gamma_beta(v);
            rho = rho_at_break * std::pow(gb / gb_at_break, -alpha_tail);
        }
        
        return std::max(rho_floor, rho);
    }
    
    /**
     * @brief Velocity at given radius (homologous expansion)
     */
    real velocity(real r) const {
        real v = r / t0;
        return std::min(v, 0.99 * C_LIGHT);
    }
    
    /**
     * @brief Pressure from density (cold EOS)
     */
    real pressure(real rho) const {
        real cs_target = 0.01 * C_LIGHT;
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

}  // anonymous namespace

void Solver::make_bns_cocoon_1d()
{
#if DIM != 1
    THROW_ERROR("BNS Cocoon 1D requires DIM == 1");
#else
    std::cout << "=== Initializing 1D BNS Cocoon Shock Breakout ===" << std::endl;
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
    
    // Cocoon parameters
    real E_cocoon = 1e-4;       // Cocoon energy (code units)
    real r_cocoon_frac = 0.3;   // Cocoon radius as fraction of r_min
    real Gamma_cocoon = 2.0;    // Initial Lorentz factor
    
    // Resolution
    int n_particles = 500;      // Number of particles
    
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
    read_real("E_cocoon", E_cocoon);
    read_real("r_cocoon_frac", r_cocoon_frac);
    read_real("Gamma_cocoon", Gamma_cocoon);
    read_int("n_particles", n_particles);
    
    // =========================================================================
    // Setup ejecta profile
    // =========================================================================
    
    EjectaProfile1D ejecta;
    ejecta.M_ej = M_ej;
    ejecta.v_min = v_min;
    ejecta.v_max = v_max;
    ejecta.v_break = v_break;
    ejecta.beta = beta;
    ejecta.alpha_tail = alpha_tail;
    ejecta.t0 = t0;
    ejecta.gamma_eos = gamma_eos;
    ejecta.rho_floor = rho_floor;
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
    
    std::cout << "\n--- Cocoon Configuration ---" << std::endl;
    std::cout << "  E_cocoon = " << E_cocoon << std::endl;
    std::cout << "  r_cocoon = " << r_cocoon << std::endl;
    std::cout << "  Γ_cocoon = " << Gamma_cocoon << std::endl;
    
    std::cout << "\n--- Resolution ---" << std::endl;
    std::cout << "  n_particles = " << n_particles << std::endl;
    
    // =========================================================================
    // Create particles
    // =========================================================================
    
    auto& particles = m_sim->get_particles();
    particles.clear();
    
    // Domain extends beyond r_max
    real domain_max = 1.2 * r_max;
    real dx = domain_max / n_particles;
    real particle_mass = M_ej / n_particles;
    real initial_sml = 2.0 * dx;
    
    std::cout << "\n--- Creating particles ---" << std::endl;
    std::cout << "  Domain: [0, " << domain_max << "]" << std::endl;
    std::cout << "  dx = " << dx << std::endl;
    
    int n_ejecta = 0, n_cocoon = 0, n_ambient = 0;
    
    for (int i = 0; i < n_particles; ++i) {
        real r = (i + 0.5) * dx;
        
        SPHParticle p;
        p.id = i;
        p.pos[0] = r;
        p.vel[0] = 0.0;
        p.acc[0] = 0.0;
        p.mass = particle_mass;
        p.sml = initial_sml;
        
        if (r < r_cocoon) {
            // ==== Cocoon region ====
            n_cocoon++;
            
            real V_cocoon = (4.0/3.0) * PI * std::pow(r_cocoon, 3);
            p.dens = 0.01 * ejecta.density(r_min);
            p.pres = E_cocoon / V_cocoon;
            p.ene = ejecta.internal_energy(p.pres, p.dens);
            
            real v_cocoon = std::sqrt(1.0 - 1.0/(Gamma_cocoon * Gamma_cocoon));
            p.vel[0] = v_cocoon;
            
        } else if (r < r_min) {
            // ==== Inner transition ====
            n_ejecta++;
            
            real frac = r / r_min;
            p.dens = frac * ejecta.density(r_min) + (1.0 - frac) * 0.1 * ejecta.density(r_min);
            p.pres = ejecta.pressure(p.dens);
            p.ene = ejecta.internal_energy(p.pres, p.dens);
            p.vel[0] = frac * ejecta.velocity(r_min);
            
        } else if (r <= r_max) {
            // ==== Ejecta region ====
            n_ejecta++;
            
            p.dens = ejecta.density(r);
            p.pres = ejecta.pressure(p.dens);
            p.ene = ejecta.internal_energy(p.pres, p.dens);
            p.vel[0] = ejecta.velocity(r);
            
        } else {
            // ==== Ambient medium ====
            n_ambient++;
            
            p.dens = rho_floor;
            p.pres = ejecta.pressure(p.dens);
            p.ene = ejecta.internal_energy(p.pres, p.dens);
            p.vel[0] = 0.0;
        }
        
        // SR-GSPH fields
        p.nu = p.mass;
        p.gradh = 1.0;
        p.alpha = m_param->av.alpha;
        p.balsara = 1.0;
        p.phi = 0.0;
        p.sound = std::sqrt(gamma_eos * p.pres / std::max(p.dens, rho_floor));
        
        particles.push_back(p);
    }
    
    std::cout << "\n--- Particle Summary ---" << std::endl;
    std::cout << "  Total: " << particles.size() << std::endl;
    std::cout << "  Cocoon: " << n_cocoon << std::endl;
    std::cout << "  Ejecta: " << n_ejecta << std::endl;
    std::cout << "  Ambient: " << n_ambient << std::endl;
    
    m_sim->set_particle_num(particles.size());
    
    // Resolution quality
    real h_over_r_min = initial_sml / r_min;
    std::cout << "\n--- Resolution Quality ---" << std::endl;
    std::cout << "  h/r_min = " << h_over_r_min << " (should be < 0.1)" << std::endl;
    
    std::cout << "\n=== BNS Cocoon 1D initialization complete ===" << std::endl;
    
#endif  // DIM == 1
}

}  // namespace sph
