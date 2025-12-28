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

    real phi = 0.0;      // gravitational potential
    vec_t grav_acc;      // gravitational acceleration (for gravity-aware Riemann solver)

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

    // === MHD (Magnetohydrodynamics) variables ===
    // Based on Iwasaki & Inutsuka (2011) GSPMHD formulation
    // Uses vec3_t (always 3D) so MHD works correctly even with DIM=1
    vec3_t B;              // Magnetic field vector (always 3D)
    vec3_t dB;             // dB/dt time derivative
    vec3_t dB_old;         // dB/dt from previous timestep (for predictor-corrector)
    real div_B = 0.0;      // Divergence of B (for monitoring)

    // MHD velocity components (for 1D simulations that need 3D velocity)
    // In DIM=1, vel only has x component; vy, vz stored separately
    real vy_mhd = 0.0;     // y-component of velocity for MHD
    real vz_mhd = 0.0;     // z-component of velocity for MHD
    real acc_y_mhd = 0.0;  // y-component of acceleration for MHD
    real acc_z_mhd = 0.0;  // z-component of acceleration for MHD

    // === Particle type flags ===
    bool is_ghost = false;  // Ghost/boundary particle (fixed, not integrated)

    int id;
    int neighbor;
    SPHParticle *next = nullptr;
};

}