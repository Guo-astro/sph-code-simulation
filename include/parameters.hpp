#pragma once

#include "defines.hpp"

namespace sph
{

enum struct SPHType {
    SSPH,
    DISPH,
    GSPH,
    GDISPH,
    SRGSPH,  // Special Relativistic GSPH
};

enum struct KernelType {
    CUBIC_SPLINE,
    WENDLAND,
    GAUSSIAN,
    UNKNOWN,
};

enum struct RiemannSolverType {
    HLL,        // Harten-Lax-van Leer approximate solver (default for GSPH, fast, non-iterative)
    ITERATIVE,  // Iterative solver from van Leer (1997) (exact for ideal gas EOS)
    EXACT,      // Exact Riemann solver (Kitajima formulation for SR-GSPH, Brent's method)
    KITAJIMA,   // Kitajima-style iterative solver (Newton-Raphson with shock/rarefaction handling)
};

enum struct BoundaryType {
    REFLECTING,  // Wall boundary: velocity is reflected (v -> -v)
    OUTFLOW,     // Open/outflow boundary: velocity is copied (waves exit without reflection)
    INFLOW,      // Inflow boundary: ghosts keep their initial state (simulates infinite domain)
};

struct SPHParameters {

    struct Time {
        real start;
        real end;
        real output;
        real energy;
    } time;

    SPHType type;

    struct CFL {
        real sound;
        real force;
    } cfl;

    struct ArtificialViscosity {
        real alpha;
        bool use_balsara_switch;
        bool use_time_dependent_av;
        real alpha_max;
        real alpha_min;
        real epsilon; // tau = h / (epsilon * c)
    } av;

    struct ArtificialConductivity {
        real alpha;
        bool is_valid;
    } ac;

    struct Tree {
        int max_level;
        int leaf_particle_num;
    } tree;

    struct Physics {
        int neighbor_number;
        real gamma;
    } physics;

    KernelType kernel;

    bool iterative_sml;

    struct Periodic {
        bool is_valid;
        real range_max[DIM];
        real range_min[DIM];
        BoundaryType boundary_type;  // Type of ghost boundary: REFLECTING or OUTFLOW
    } periodic;

    struct Gravity {
        bool is_valid;
        real constant;
        real theta;
    } gravity;

    struct GSPH {
        bool is_2nd_order;
        RiemannSolverType riemann_solver;  // Riemann solver type: HLL (default) or ITERATIVE
    } gsph;

    struct SRGSPH {
        bool is_2nd_order;        // Enable MUSCL reconstruction
        real c_speed;             // Speed of light (default: 1.0)
        real c_shock;             // Shock detection parameter (default: 3.0)
        real c_cd;                // Contact discontinuity parameter (default: 0.2)
        real eta;                 // Smoothing length parameter (default: 1.0)
        real smoothing_length;    // Fixed smoothing length h (§2.2, constant for all particles)
    } srgsph;
};

}