#pragma once

#include "vector_type.hpp"

namespace sph
{

class SPHParticle {
public:
    vec_t pos;    // position
    vec_t vel;    // velocity
    vec_t vel_p;  // velocity at t + dt/2
    vec_t acc;    // acceleration
    real mass;    // mass
    real dens;    // mass density
    real pres;    // pressure
    real ene;     // internal energy
    real ene_p;   // internal energy at t + dt/2
    real dene;    // du/dt
    real sml;     // smoothing length
    real sound;   // sound speed

    real balsara; // balsara switch
    real alpha;   // AV coefficient

    real gradh;   // grad-h term

    real phi = 0.0; // potential

    // Special Relativistic variables (only used when SPHType == SRGSPH)
    vec_t S;         // canonical momentum (S = γHv)
    real e;          // canonical energy (e = γH - P/(Nc²))
    real N;          // baryon number density (lab frame)
    real gamma_lor;  // Lorentz factor γ = 1/√(1-v²/c²)
    real enthalpy;   // specific enthalpy H
    real nu;         // baryon number in particle
    vec_t dS;        // dS/dt (time derivative of canonical momentum)
    real de;         // de/dt (time derivative of canonical energy)

    int id;
    int neighbor;
    SPHParticle *next = nullptr;
};

}