#pragma once

#include "defines.hpp"
#include "vector_type.hpp"
#include <string>

namespace sph
{

enum struct SPHType {
    SSPH,
    DISPH,
    GSPH,
    GDISPH,
    SRGSPH,  // Special Relativistic GSPH
    GRGSPH,  // General Relativistic GSPH (GR-GSPH)
};

enum struct KernelType {
    CUBIC_SPLINE,
    WENDLAND,
    GAUSSIAN,
    UNKNOWN,
};

enum struct RiemannSolverType {
    HLL,        // Harten-Lax-van Leer approximate solver (default for GSPH, fast, non-iterative)
    HLLC,       // HLLC solver with contact wave (Mignone & Bodo 2005 for relativistic flows)
    ITERATIVE,  // Iterative solver from van Leer (1997) (exact for ideal gas EOS)
    EXACT,      // Exact Riemann solver (Kitajima formulation for SR-GSPH, Brent's method)
    KITAJIMA,   // Kitajima-style iterative solver (Newton-Raphson with shock/rarefaction handling)
};

enum struct BoundaryType {
    REFLECTING,  // Wall boundary: velocity is reflected (v -> -v)
    OUTFLOW,     // Open/outflow boundary: velocity is copied (waves exit without reflection)
    INFLOW,      // Inflow boundary: ghosts keep their initial state (simulates infinite domain)
};

enum struct GravitySofteningType {
    HERNQUIST_KATZ,  // Original: Hernquist & Katz (1989) spline, ε = h/2
    WENDLAND_C4,     // True kernel-convolved: use same Wendland C4 as SPH
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
        real c_smooth;  // Smoothing length expansion factor for h-adaptation (default=1.0, typical=2.0)
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
        bool use_fixed_softening;  // If true, use fixed softening instead of h-dependent
        real fixed_softening;      // Fixed softening length (only used if use_fixed_softening=true)
        GravitySofteningType softening_type;  // HERNQUIST_KATZ (default) or WENDLAND_C4
        bool use_kernel_gravity_1d;  // For 1D: use kernel-convolved gravity (default: true)
        bool use_kernel_gravity_2d;  // For 2D: use kernel-convolved gravity (default: true)
        bool use_kernel_gravity_planar_2d;  // For 2D planar slab: 1D gravity in y-direction
        bool use_kernel_gravity_cylinder_3d;  // For 3D cylinder: 2D radial gravity in xy-plane
    } gravity;

    struct GSPH {
        bool is_2nd_order;
        RiemannSolverType riemann_solver;  // Riemann solver type: HLL (default) or ITERATIVE
        bool use_gradh;  // Enable grad-h correction (default: true). Disabling causes core collapse in hydrostatic tests.
    } gsph;

    struct SRGSPH {
        bool is_2nd_order;        // Enable MUSCL reconstruction
        real c_speed;             // Speed of light (default: 1.0)
        real c_shock;             // Shock detection parameter (default: 3.0)
        real c_cd;                // Contact discontinuity parameter (default: 0.2)
        real eta;                 // Smoothing length parameter (default: 1.0)
        real smoothing_length;    // Fixed smoothing length h (§2.2, constant for all particles)
        RiemannSolverType riemann_solver;  // EXACT (default) or HLLC
    } srgsph;

    struct GRGSPH {
        std::string metric_type;  // "minkowski", "schwarzschild", or "kerr"
        real bh_mass;             // Black hole mass M (default: 1.0)
        real bh_spin;             // Black hole spin a (for Kerr, default: 0.0)
        bool geodesic_mode;       // Geodesic test mode: skip pairwise SPH forces
    } grgsph;

    struct Thermal {
        bool enable_cooling;      // Enable ISM cooling/heating
        real N_H_column;          // Column density [cm^-2] (1e19 or 1e20)
        real relaxation_time;     // Thermal relaxation timescale [code time units]
        real density_to_n_H;      // Conversion factor: code density -> n_H [cm^-3]
        real u_to_cgs;            // Conversion factor: code energy -> erg/g
        real t_to_cgs;            // Conversion factor: code time -> seconds
    } thermal;

    struct ExternalBH {
        bool enabled;             // Enable external point-mass BH
        real mass;                // BH mass in code units
        real softening;           // Softening length in code units
        vec_t position;           // BH position in code units
        real G_constant;          // Gravitational constant (usually 1.0)
    } external_bh;
};

}