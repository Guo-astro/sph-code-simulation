#include "srgsph/sr_fluid_force.hpp"
#include "particle.hpp"
#include "exception.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

template <typename T>
bool nearly_equal(T a, T b, T tol = static_cast<T>(1.0e-12))
{
    return std::abs(a - b) <= tol;
}

void given_particle_with_derivatives(sph::SPHParticle & particle, real nu)
{
    particle = sph::SPHParticle{};
    particle.nu = nu;
    particle.de = 6.0;
    particle.dene = -1.0; // sentinel
    particle.acc = vec_t(-1.0);

    for(int dim = 0; dim < DIM; ++dim) {
        particle.dS[dim] = static_cast<real>((dim + 1) * 4.0);
        particle.vel[dim] = 0.0;
    }
}

} // namespace

int main()
{
    using sph::SPHParticle;
    using sph::srgsph::FluidForce;

    // Given: SR particle derivatives expressed in <nu dS> form
    SPHParticle particle{};
    given_particle_with_derivatives(particle, 2.0);

    vec_t original_dS = particle.dS;
    const real original_de = particle.de;

    // When: we normalize by baryon number (per paper Eq. 425)
    FluidForce::normalize_sr_derivatives(particle);

    const real expected_scale = 0.5; // 1/nu
    bool pass = true;
    for(int dim = 0; dim < DIM; ++dim) {
        if(!nearly_equal(particle.dS[dim], original_dS[dim] * expected_scale)) {
            std::cerr << "FAIL: dS component " << dim << " mismatch" << std::endl;
            pass = false;
        }
        if(!nearly_equal(particle.acc[dim], original_dS[dim] * expected_scale)) {
            std::cerr << "FAIL: acc alias component " << dim << " mismatch" << std::endl;
            pass = false;
        }
    }
    if(!nearly_equal(particle.de, original_de * expected_scale)) {
        std::cerr << "FAIL: de scaling mismatch" << std::endl;
        pass = false;
    }
    if(!nearly_equal(particle.dene, original_de * expected_scale)) {
        std::cerr << "FAIL: dene alias mismatch" << std::endl;
        pass = false;
    }

    if(!pass) {
        return EXIT_FAILURE;
    }

    // Then: non-positive baryon numbers trigger a descriptive exception
    bool threw = false;
    try {
        SPHParticle invalid{};
        given_particle_with_derivatives(invalid, 0.0);
        FluidForce::normalize_sr_derivatives(invalid);
    } catch(const sph::SPHException &) {
        threw = true;
    }

    if(!threw) {
        std::cerr << "FAIL: expected SPHException for nu <= 0" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "SRGSPH derivative scaling test passed." << std::endl;
    return EXIT_SUCCESS;
}
