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

    // OLD DERIVATIVES for predictor-corrector:
    vec_t dS_old;    // dS/dt from previous timestep
    real de_old;     // de/dt from previous timestep

    // TANGENT VELOCITY for 1D SR simulations with transverse motion
    // Used in Section 3.1.5 tangent velocity tests (Kitajima et al. 2025)
    // This is tracked as a separate scalar because in 1D the vel vector has only one component
    real vel_t = 0.0;      // tangent velocity v^t (perpendicular to x-axis)
    real dS_t = 0.0;       // dS_t/dt (tangent momentum derivative)
    real dS_t_old = 0.0;   // dS_t/dt from previous timestep
    real S_t = 0.0;        // tangent momentum S_t = γ H v^t

    // DERIVED quantities:
    real gamma_lor;  // Lorentz factor γ = 1/√(1-v²/c²)
    real enthalpy;   // specific enthalpy H = 1 + u/c² + P/(nc²)
    real nu;         // baryon number per particle (constant)

    // === Particle type flags ===
    bool is_ghost = false;  // Ghost/boundary particle (fixed, not integrated)

    int id;
    int neighbor;
    SPHParticle *next = nullptr;
};

}