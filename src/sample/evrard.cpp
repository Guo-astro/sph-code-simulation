#include <random>

#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"

namespace sph
{

// Evrard collapse (Evrard 1988)
void Solver::make_evrard()
{
    std::cout << "make_evrard() called, DIM = " << DIM << std::endl;
    std::cout.flush();
#if DIM != 3
    std::cout << "ERROR: DIM != 3, throwing error" << std::endl;
    std::cout.flush();
    THROW_ERROR("DIM != 3");
#else

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    std::cout << "N from parameters = " << N << std::endl;
    std::cout.flush();
    auto & p = m_sim->get_particles();
    const real dx = 2.0 / N;

    std::cout << "Evrard: N = " << N << ", dx = " << dx << std::endl;
    std::cout.flush();


    int count_total = 0;
    int count_added = 0;

    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            for(int k = 0; k < N; ++k) {
                count_total++;
                vec_t r = {
                    (i + 0.5) * dx - 1.0,
                    (j + 0.5) * dx - 1.0,
                    (k + 0.5) * dx - 1.0
                };
                const real r_0 = std::abs(r);
                if(r_0 > 1.0) {
                    continue;
                }

                if(r_0 > 0.0) {
                    const real r_abs = std::pow(r_0, 1.5);
                    r *= r_abs / r_0;
                }

                SPHParticle p_i;
                p_i.pos = r;
                p.emplace_back(p_i);
                count_added++;
            }
        }
    }

    std::cout << "Evrard: Grid points = " << count_total << ", Added = " << count_added << std::endl;

    const real mass = 1.0 / p.size();
    const real gamma = m_param->physics.gamma;
    const real G = m_param->gravity.constant;
    const real u = 0.05 * G;

    std::cout << "Evrard: Created " << p.size() << " particles, mass = " << mass << std::endl;

    int i = 0;
    for(auto & p_i : p) {
        p_i.vel = 0.0;
        p_i.vel = 0.0;
        p_i.mass = mass;
        p_i.dens = 1.0 / (2.0 * M_PI * std::abs(p_i.pos));
        p_i.ene = u;
        p_i.pres = (gamma - 1.0) * p_i.dens * u;
        p_i.id = i;
        i++;
    }

    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
#endif
}

}
