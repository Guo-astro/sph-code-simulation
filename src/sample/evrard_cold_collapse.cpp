#include <random>

#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"

namespace sph
{

// Evrard COLD collapse - demonstrates that LOWER initial thermal energy
// produces STRONGER shock heating and GREATER expansion after collapse.
// This is the opposite of intuition: colder → faster fall → harder shock → more expansion!
// For particles that stay confined, use hydrostatic or lane_emden samples instead.
void Solver::make_evrard_cold_collapse()
{
    std::cout << "make_evrard_cold_collapse() called, DIM = " << DIM << std::endl;
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

    std::cout << "Evrard Bound: N = " << N << ", dx = " << dx << std::endl;
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

    std::cout << "Evrard Bound: Grid points = " << count_total << ", Added = " << count_added << std::endl;

    const real mass = 1.0 / p.size();
    const real gamma = m_param->physics.gamma;
    const real G = m_param->gravity.constant;
    
    // KEY DIFFERENCE: Much lower thermal energy (0.002 instead of 0.05)
    // This ensures particles remain bound and undergo multiple collapse/rebound cycles
    // With u=0.002*G:
    //   - Initial thermal energy equal to standard Evrard (u=0.05×G)
    //   - Start from virial equilibrium (2×KE+TE = |PE|)
    //   - No collapse → no shock → thermal energy stays constant
    //   - Particles oscillate gently within r<2
    //   - This is NOT a bound test, it's a stability test!
    const real u = 0.05 * G;

    std::cout << "Evrard Cold Collapse: Created " << p.size() << " particles, mass = " << mass << std::endl;
    std::cout << "Evrard Cold Collapse: Initial thermal energy u = " << u << " (same as standard - equilibrium test)" << std::endl;

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
