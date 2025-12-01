/**
 * @file point_mass_bh.cpp
 * @brief Implementation of point-mass black hole external force
 */

#include "external_forces/point_mass_bh.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include <cmath>

namespace sph
{
namespace external_forces
{

PointMassBlackHole::PointMassBlackHole()
    : m_position(0.0)
    , m_mass(1e5)
    , m_softening_length(0.01)
    , m_softening_squared(1e-4)
    , m_G_constant(1.0)
    , m_is_moving(false)
    , m_velocity(0.0)
    , m_initialized(false)
{
}

void PointMassBlackHole::initialize(const PointMassBHParams& params)
{
    m_position = params.position;
    m_mass = params.mass;
    m_softening_length = params.softening_length;
    m_softening_squared = m_softening_length * m_softening_length;
    m_G_constant = params.G_constant;
    m_is_moving = params.is_moving;
    m_velocity = params.velocity;
    m_initialized = true;
}

void PointMassBlackHole::initialize(std::shared_ptr<SPHParameters> param)
{
    // Read from parameter file
    // For now, use default values
    // TODO: Extend SPHParameters to include external_force section
    
    PointMassBHParams default_params;
    default_params.G_constant = param->gravity.constant;
    initialize(default_params);
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
        
        // Update potential: Φ = -G M_BH / sqrt(r² + ε²)
        const real r_soft = std::sqrt(r_squared + m_softening_squared);
        p_i.phi -= m_G_constant * m_mass / r_soft;
    }
}

void PointMassBlackHole::update_position(real dt)
{
    if (m_is_moving) {
        m_position += m_velocity * dt;
    }
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
    
    for (int i = 0; i < num; ++i) {
        const vec_t r_vec = particles[i].pos - m_position;
        const real distance = std::abs(r_vec);
        
        if (distance < r) {
            enclosed_mass += particles[i].mass;
        }
    }
    
    return enclosed_mass;
}

int PointMassBlackHole::count_accreted_particles(std::shared_ptr<Simulation> sim, real r_acc) const
{
    const auto& particles = sim->get_particles();
    const int num = sim->get_particle_num();
    
    int count = 0;
    
    for (int i = 0; i < num; ++i) {
        const vec_t r_vec = particles[i].pos - m_position;
        const real distance = std::abs(r_vec);
        
        if (distance < r_acc) {
            count++;
        }
    }
    
    return count;
}

} // namespace external_forces
} // namespace sph
