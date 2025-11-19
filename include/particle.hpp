#pragma once

#include "vector_type.hpp"

namespace sph
{

class SPHParticle {
public:
    // === Standard fields (used by all SPH methods) ===
    vec_t pos;    // position
    vec_t vel;    // velocity (SRGSPH: primitive velocity v, recovered from S)
    vec_t vel_p;  // velocity at t + dt/2
    vec_t acc;    // acceleration (SRGSPH: ALIAS for dS, do NOT use directly!)
    real mass;    // mass (SRGSPH: baryon number ν)
    real dens;    // mass density (SRGSPH: lab-frame density N = γn)
    real pres;    // pressure
    real ene;     // internal energy (SRGSPH: primitive u = P/[(γ-1)n])
    real ene_p;   // internal energy at t + dt/2
    real dene;    // du/dt (SRGSPH: ALIAS for de, do NOT use directly!)
    real sml;     // smoothing length
    real sound;   // sound speed

    real balsara; // balsara switch
    real alpha;   // AV coefficient

    real gradh;   // grad-h term

    real phi = 0.0; // potential

    // === Special Relativistic variables (SRGSPH only) ===
    // CONSERVED variables (time-evolved):
    vec_t S;         // canonical momentum S = γHv (per baryon)
    real e;          // canonical energy e = γH - P/(Nc²) (per baryon)
    real N;          // lab-frame baryon density N = γn
    
    // TIME DERIVATIVES of conserved variables:
    vec_t dS;        // dS/dt (canonical momentum derivative, COPIED to acc)
    real de;         // de/dt (canonical energy derivative, COPIED to dene)
    
    // DERIVED quantities:
    real gamma_lor;  // Lorentz factor γ = 1/√(1-v²/c²)
    real enthalpy;   // specific enthalpy H = 1 + u/c² + P/(nc²)
    real nu;         // baryon number per particle (constant)

    int id;
    int neighbor;
    SPHParticle *next = nullptr;
};

}