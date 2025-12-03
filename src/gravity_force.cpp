#include "defines.hpp"
#include "gravity_force.hpp"
#include "particle.hpp"
#include "periodic.hpp"
#include "simulation.hpp"
#include "bhtree.hpp"

#include <iostream>
#include <atomic>

#ifdef EXHAUSTIVE_SEARCH
#include "exhaustive_search.hpp"
#endif

namespace sph
{

 // Hernquist & Katz (1989)
inline real f(const real r, const real h)
{
    const real e = h * 0.5;
    const real u = r / e;
    
    if(u < 1.0) {
        return (-0.5 * u * u * (1.0 / 3.0 - 3.0 / 20 * u * u + u * u * u / 20) + 1.4) / e;
    } else if(u < 2.0) {
        return -1.0 / (15 * r) + (-u * u * (4.0 / 3.0 - u + 0.3 * u * u - u * u * u / 30) + 1.6) / e;
    } else {
        return 1 / r;
    }
}

inline real g(const real r, const real h)
{
    const real e = h * 0.5;
    const real u = r / e;
    
    if(u < 1.0) {
        return (4.0 / 3.0 - 1.2 * u * u + 0.5 * u * u * u) / (e * e * e);
    } else if(u < 2.0) {
        return (-1.0 / 15 + 8.0 / 3 * u * u * u - 3 * u * u * u * u + 1.2 * u * u * u * u * u - u * u * u * u * u * u / 6.0) / (r * r * r);
    } else {
        return 1 / (r * r * r);
    }
}

void GravityForce::initialize(std::shared_ptr<SPHParameters> param)
{
    m_is_valid = param->gravity.is_valid;
    if(m_is_valid) {
        m_constant = param->gravity.constant;
    }
}

void GravityForce::calculation(std::shared_ptr<Simulation> sim)
{
    if(!m_is_valid) {
        return;
    }

    auto & particles = sim->get_particles();
    const int num = sim->get_particle_num();
    
    // Debug output
    static std::atomic<int> grav_debug_count{0};
    bool should_debug = grav_debug_count.load() < 3;
    
#ifdef EXHAUSTIVE_SEARCH
    auto * periodic = sim->get_periodic().get();
#else
    auto * tree = sim->get_tree().get();
#endif

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        
#ifdef EXHAUSTIVE_SEARCH
        real phi = 0.0;
        vec_t force(0.0);
        const vec_t & r_i = p_i.pos;

        for(int j = 0; j < num; ++j) {
            const auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
            const real r = std::abs(r_ij);
            phi -= m_constant * p_j.mass * (f(r, p_i.sml) + f(r, p_j.sml)) * 0.5;
            force -= r_ij * (m_constant * p_j.mass * (g(r, p_i.sml) + g(r, p_j.sml)) * 0.5);
        }

        p_i.grav_acc = force;  // Store gravity acceleration separately
        // Note: grav_acc is NOT added to acc here - that's done after fluid force
        p_i.phi = phi;
#else
        tree->tree_force(p_i);
#endif
    }
    
    // Debug: Check particle 0's grav_acc after computation
    if (should_debug) {
        int expected = grav_debug_count.load();
        if (grav_debug_count.compare_exchange_strong(expected, expected + 1)) {
            std::cerr << "[GRAVITY-CALC] #" << expected 
                      << " t=" << sim->get_time()
                      << " grav_acc[0]=(" << particles[0].grav_acc[0] 
                      << "," << particles[0].grav_acc[1] 
                      << "," << particles[0].grav_acc[2] << ")"
                      << " |grav_acc|=" << std::abs(particles[0].grav_acc)
                      << std::endl;
        }
    }
}

}
