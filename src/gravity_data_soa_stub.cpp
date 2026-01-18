/**
 * @file gravity_data_soa_stub.cpp
 * @brief Implementation of Structure-of-Arrays gravity data
 *
 * Provides SoA layout for cache-efficient gravity computation.
 */

#include "gravity_data_soa.hpp"
#include "analytic_gravity_test.hpp"  // For GravitySofteningType
#include <cmath>

namespace sph {

GravityDataSoA::GravityDataSoA(size_t n) {
    if (n > 0) {
        resize(n);
    }
}

GravityDataSoA GravityDataSoA::from_aos(const std::vector<SPHParticle>& particles) {
    size_t n = particles.size();
    GravityDataSoA soa(n);

    for (size_t i = 0; i < n; ++i) {
        soa.pos_x[i] = particles[i].pos[0];
#if DIM >= 2
        soa.pos_y[i] = particles[i].pos[1];
#endif
#if DIM == 3
        soa.pos_z[i] = particles[i].pos[2];
#endif
        soa.mass[i] = particles[i].mass;
        soa.sml[i] = particles[i].sml;
    }

    return soa;
}

std::vector<SPHParticle> GravityDataSoA::to_aos() const {
    size_t n = size();
    std::vector<SPHParticle> particles(n);

    for (size_t i = 0; i < n; ++i) {
        particles[i].pos[0] = pos_x[i];
#if DIM >= 2
        particles[i].pos[1] = pos_y[i];
#endif
#if DIM == 3
        particles[i].pos[2] = pos_z[i];
#endif
        particles[i].mass = mass[i];
        particles[i].sml = sml[i];
        particles[i].grav_acc[0] = grav_acc_x[i];
#if DIM >= 2
        particles[i].grav_acc[1] = grav_acc_y[i];
#endif
#if DIM == 3
        particles[i].grav_acc[2] = grav_acc_z[i];
#endif
        particles[i].phi = phi[i];
    }

    return particles;
}

void GravityDataSoA::copy_results_to_aos(std::vector<SPHParticle>& particles) const {
    size_t n = std::min(size(), particles.size());
    for (size_t i = 0; i < n; ++i) {
        particles[i].grav_acc[0] = grav_acc_x[i];
#if DIM >= 2
        particles[i].grav_acc[1] = grav_acc_y[i];
#endif
#if DIM == 3
        particles[i].grav_acc[2] = grav_acc_z[i];
#endif
        particles[i].phi = phi[i];
    }
}

void GravityDataSoA::resize(size_t n) {
    pos_x.resize(n);
#if DIM >= 2
    pos_y.resize(n);
#endif
#if DIM == 3
    pos_z.resize(n);
#endif
    mass.resize(n);
    sml.resize(n);
    grav_acc_x.resize(n);
#if DIM >= 2
    grav_acc_y.resize(n);
#endif
#if DIM == 3
    grav_acc_z.resize(n);
#endif
    phi.resize(n);
}

void GravityDataSoA::clear() {
    pos_x.clear();
#if DIM >= 2
    pos_y.clear();
#endif
#if DIM == 3
    pos_z.clear();
#endif
    mass.clear();
    sml.clear();
    grav_acc_x.clear();
#if DIM >= 2
    grav_acc_y.clear();
#endif
#if DIM == 3
    grav_acc_z.clear();
#endif
    phi.clear();
}

void GravityDataSoA::allocate_aligned(size_t n) {
    resize(n);
}

vec_t compute_gravity_soa_single(int target_idx,
                                  const GravityDataSoA& data,
                                  real softening)
{
    vec_t force(0.0);
    const real G = 1.0;
    size_t n = data.size();

    real target_x = data.pos_x[target_idx];
#if DIM >= 2
    real target_y = data.pos_y[target_idx];
#endif
#if DIM == 3
    real target_z = data.pos_z[target_idx];
#endif

    for (size_t j = 0; j < n; ++j) {
        if (static_cast<int>(j) == target_idx) continue;

        real dx = target_x - data.pos_x[j];
#if DIM >= 2
        real dy = target_y - data.pos_y[j];
#endif
#if DIM == 3
        real dz = target_z - data.pos_z[j];
#endif

#if DIM == 1
        real r2 = dx * dx;
#elif DIM == 2
        real r2 = dx * dx + dy * dy;
#elif DIM == 3
        real r2 = dx * dx + dy * dy + dz * dz;
#endif
        (void)r2; // r2 used directly in r_soft calculation

        // Softened inverse square law
        real r_soft = std::sqrt(r2 + softening * softening);
        real factor = G * data.mass[j] / (r_soft * r_soft * r_soft);

        force[0] -= dx * factor;
#if DIM >= 2
        force[1] -= dy * factor;
#endif
#if DIM == 3
        force[2] -= dz * factor;
#endif
    }

    return force;
}

vec_t compute_gravity_soa_at_position(const vec_t& pos,
                                       const GravityDataSoA& data,
                                       real softening)
{
    vec_t force(0.0);
    const real G = 1.0;
    size_t n = data.size();

    for (size_t j = 0; j < n; ++j) {
        real dx = pos[0] - data.pos_x[j];
#if DIM >= 2
        real dy = pos[1] - data.pos_y[j];
#endif
#if DIM == 3
        real dz = pos[2] - data.pos_z[j];
#endif

#if DIM == 1
        real r2 = dx * dx;
#elif DIM == 2
        real r2 = dx * dx + dy * dy;
#elif DIM == 3
        real r2 = dx * dx + dy * dy + dz * dz;
#endif

        // Softened inverse square law
        real r_soft = std::sqrt(r2 + softening * softening);
        real factor = G * data.mass[j] / (r_soft * r_soft * r_soft);

        force[0] -= dx * factor;
#if DIM >= 2
        force[1] -= dy * factor;
#endif
#if DIM == 3
        force[2] -= dz * factor;
#endif
    }

    return force;
}

void compute_gravity_soa(GravityDataSoA& data, real G,
                         GravitySofteningType type, real softening)
{
    (void)type;  // Using simple softening for now
    size_t n = data.size();

    // Zero output arrays
    for (size_t i = 0; i < n; ++i) {
        data.grav_acc_x[i] = 0.0;
#if DIM >= 2
        data.grav_acc_y[i] = 0.0;
#endif
#if DIM == 3
        data.grav_acc_z[i] = 0.0;
#endif
        data.phi[i] = 0.0;
    }

    // O(N^2) direct summation
    for (size_t i = 0; i < n; ++i) {
        real pos_i_x = data.pos_x[i];
#if DIM >= 2
        real pos_i_y = data.pos_y[i];
#endif
#if DIM == 3
        real pos_i_z = data.pos_z[i];
#endif

        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;

            real dx = pos_i_x - data.pos_x[j];
#if DIM >= 2
            real dy = pos_i_y - data.pos_y[j];
#endif
#if DIM == 3
            real dz = pos_i_z - data.pos_z[j];
#endif

#if DIM == 1
            real r2 = dx * dx;
#elif DIM == 2
            real r2 = dx * dx + dy * dy;
#elif DIM == 3
            real r2 = dx * dx + dy * dy + dz * dz;
#endif

            real r_soft = std::sqrt(r2 + softening * softening);
            real factor = G * data.mass[j] / (r_soft * r_soft * r_soft);

            data.grav_acc_x[i] -= dx * factor;
#if DIM >= 2
            data.grav_acc_y[i] -= dy * factor;
#endif
#if DIM == 3
            data.grav_acc_z[i] -= dz * factor;
#endif

            // Potential
            data.phi[i] -= G * data.mass[j] / r_soft;
        }
    }
}

} // namespace sph
