/**
 * @file point_mass_bh.cpp
 * @brief Implementation of fixed point-mass black hole with sink accretion
 *
 * Implements:
 *   - Plummer-softened gravity: a = -G M_BH r / (r² + ε²)^(3/2)
 *   - Sink accretion with inflow + bound criteria
 *   - Accretion diagnostics tracking
 */

#include "external_forces/point_mass_bh.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "logger.hpp"
#include <cmath>
#include <algorithm>
#include <vector>
#include <functional>

namespace sph
{
namespace external_forces
{

PointMassBlackHole::PointMassBlackHole()
    : m_position(0.0)
    , m_mass(1e5)
    , m_softening_length(0.001)
    , m_softening_squared(1e-6)
    , m_sink_radius(0.005)
    , m_sink_radius_squared(2.5e-5)
    , m_G_constant(1.0)
    , m_is_moving(false)
    , m_velocity(0.0)
    , m_enable_sink(true)
    , m_initialized(false)
    , m_accreted_mass(0.0)
    , m_accreted_count(0)
{
}

void PointMassBlackHole::initialize(const PointMassBHParams& params)
{
    m_position = params.position;
    m_mass = params.mass;
    m_softening_length = params.softening_length;
    m_softening_squared = m_softening_length * m_softening_length;
    m_sink_radius = params.sink_radius;
    m_sink_radius_squared = m_sink_radius * m_sink_radius;
    m_G_constant = params.G_constant;
    m_is_moving = params.is_moving;
    m_velocity = params.velocity;
    m_enable_sink = params.enable_sink;
    m_initialized = true;

    // Reset accretion counters on initialization
    m_accreted_mass = 0.0;
    m_accreted_count = 0;

    WRITE_LOG << "PointMassBlackHole initialized:";
    WRITE_LOG << "  Mass = " << m_mass << " (code units)";
    WRITE_LOG << "  Position = (" << m_position[0] << ", " << m_position[1]
#if DIM == 3
              << ", " << m_position[2]
#endif
              << ")";
    WRITE_LOG << "  Softening epsilon = " << m_softening_length;
    WRITE_LOG << "  Sink radius = " << m_sink_radius;
    WRITE_LOG << "  Sink enabled = " << (m_enable_sink ? "yes" : "no");
    WRITE_LOG << "  Is moving = " << (m_is_moving ? "yes" : "no");
}

void PointMassBlackHole::initialize(std::shared_ptr<SPHParameters> param)
{
    // Read from SPHParameters external_bh section
    PointMassBHParams bh_params;
    bh_params.G_constant = param->gravity.constant;
    bh_params.mass = param->external_bh.mass;
    bh_params.softening_length = param->external_bh.softening;
    bh_params.position = param->external_bh.position;
    // Use default sink radius (can be extended in SPHParameters if needed)
    bh_params.sink_radius = 0.005;
    bh_params.enable_sink = true;
    bh_params.is_moving = false;
    bh_params.velocity = vec_t(0.0);

    initialize(bh_params);
}

void PointMassBlackHole::calculation(std::shared_ptr<Simulation> sim)
{
    if (!m_initialized) {
        return;
    }
    
    auto& particles = sim->get_particles();
    const int num = sim->get_particle_num();
    
#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        auto& p_i = particles[i];
        
        // Vector from BH to particle
        const vec_t r_vec = p_i.pos - m_position;
        const real r_squared = abs2(r_vec);
        
        // Softened force: F = -G M_BH r_vec / (r² + ε²)^(3/2)
        const real force_factor = softened_force_factor(r_squared);
        const vec_t acc_bh = -r_vec * force_factor;
        
        // Add to particle acceleration
        p_i.acc += acc_bh;
        
        // NOTE: Do NOT modify p_i.phi here!
        // p.phi should only contain self-gravity potential (from BH tree).
        // External BH potential is computed separately in energy calculations
        // with the correct factor (no 0.5 for external potential).
        // See: Solver::compute_energy() and OutputManager::compute_energies()
    }
}

void PointMassBlackHole::update_position(real dt)
{
    if (m_is_moving) {
        m_position += m_velocity * dt;
    }
}

int PointMassBlackHole::apply_sink_accretion(std::shared_ptr<Simulation> sim)
{
    if (!m_initialized || !m_enable_sink) {
        return 0;
    }

    auto& particles = sim->get_particles();
    int num = sim->get_particle_num();

    // Collect indices of particles to remove
    std::vector<int> to_remove;
    to_remove.reserve(100);  // Preallocate for efficiency

    for (int i = 0; i < num; ++i) {
        const auto& p_i = particles[i];

        // Skip ghost particles
        if (p_i.is_ghost) continue;

        // Vector from BH to particle
        const vec_t r_vec = p_i.pos - m_position;
        const real r_squared = abs2(r_vec);

        // === NaN/Inf guard: skip corrupted particles ===
        // When cooling causes runaway collapse, positions can become NaN.
        // NaN comparisons always return false, so without this guard,
        // NaN-position particles would pass all checks and be incorrectly accreted.
        if (std::isnan(r_squared) || std::isinf(r_squared)) {
            continue;  // Skip corrupted particles
        }

        // === Distance check: r < r_sink ===
        if (r_squared >= m_sink_radius_squared) {
            continue;  // Outside sink radius
        }

        const real r = std::sqrt(r_squared);

        // === Inflow check: v_rad = v · r̂ < 0 ===
        // Radial velocity: positive = outward, negative = inward
        const real v_rad = inner_product(p_i.vel, r_vec) / std::max(r, 1e-20);
        if (v_rad >= 0.0) {
            continue;  // Moving outward, do not accrete
        }

        // === Bound check: E = 0.5|v|² + Φ_BH < 0 ===
        const real v_squared = abs2(p_i.vel);
        const real phi_bh = potential(p_i.pos);
        const real E_specific = 0.5 * v_squared + phi_bh;

        if (E_specific >= 0.0) {
            continue;  // Unbound, may escape
        }

        // Particle passes all checks: mark for removal
        to_remove.push_back(i);
    }

    // Remove particles (in reverse order to maintain index validity)
    const int num_removed = static_cast<int>(to_remove.size());

    if (num_removed > 0) {
        // Track accreted mass before removal
        for (int idx : to_remove) {
            m_accreted_mass += particles[idx].mass;
        }
        m_accreted_count += num_removed;

        // Sort indices in descending order for safe removal
        std::sort(to_remove.begin(), to_remove.end(), std::greater<int>());

        // Remove particles by swapping with last and popping
        for (int idx : to_remove) {
            if (idx < static_cast<int>(particles.size()) - 1) {
                std::swap(particles[idx], particles.back());
            }
            particles.pop_back();
        }

        // Update particle count in simulation
        sim->set_particle_num(static_cast<int>(particles.size()));

        WRITE_LOG << "[SINK] Removed " << num_removed << " particles"
                  << " (total accreted: " << m_accreted_count
                  << ", mass: " << m_accreted_mass << ")";
    }

    return num_removed;
}

vec_t PointMassBlackHole::acceleration(const vec_t& r) const
{
    const vec_t r_vec = r - m_position;
    const real r_squared = abs2(r_vec);
    const real force_factor = softened_force_factor(r_squared);
    return -r_vec * force_factor;
}

real PointMassBlackHole::potential(const vec_t& r) const
{
    const vec_t r_vec = r - m_position;
    const real r_squared = abs2(r_vec);
    const real r_soft = std::sqrt(r_squared + m_softening_squared);
    return -m_G_constant * m_mass / r_soft;
}

real PointMassBlackHole::get_enclosed_mass(std::shared_ptr<Simulation> sim, real r) const
{
    const auto& particles = sim->get_particles();
    const int num = sim->get_particle_num();

    real enclosed_mass = 0.0;
    const real r_squared = r * r;

    for (int i = 0; i < num; ++i) {
        const vec_t r_vec = particles[i].pos - m_position;
        const real distance_squared = abs2(r_vec);

        if (distance_squared < r_squared) {
            enclosed_mass += particles[i].mass;
        }
    }

    return enclosed_mass;
}

int PointMassBlackHole::count_particles_within(std::shared_ptr<Simulation> sim, real r_acc) const
{
    const auto& particles = sim->get_particles();
    const int num = sim->get_particle_num();

    int count = 0;
    const real r_acc_squared = r_acc * r_acc;

    for (int i = 0; i < num; ++i) {
        const vec_t r_vec = particles[i].pos - m_position;
        const real distance_squared = abs2(r_vec);

        if (distance_squared < r_acc_squared) {
            count++;
        }
    }

    return count;
}

} // namespace external_forces
} // namespace sph
