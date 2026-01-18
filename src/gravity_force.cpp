#include "defines.hpp"
#include "gravity_force.hpp"
#include "particle.hpp"
#include "periodic.hpp"
#include "simulation.hpp"
#include "bhtree.hpp"
#include "hernquist_katz_lookup_table.hpp"
#include "softening_lookup_table.hpp"

#include <iostream>
#include <atomic>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>

#ifdef EXHAUSTIVE_SEARCH
#include "exhaustive_search.hpp"
#endif

namespace sph
{

// =============================================================================
// Hernquist & Katz (1989) softening kernels - using lookup table
// Softening length ε = h/2, support radius 2ε = h
// =============================================================================
inline real f(const real r, const real h)
{
    return HernquistKatzLookupTable::get_instance().f_full(r, h);
}

inline real g(const real r, const real h)
{
    return HernquistKatzLookupTable::get_instance().g_full(r, h);
}

// =============================================================================
// Wendland C4 kernel gravitational potential - using lookup table
// Reference: Price & Monaghan (2007), Dehnen & Aly (2012)
// =============================================================================

// Gravitational potential kernel for Wendland C4 (3D)
// Returns φ̃(r,h) such that φ_grav = -G m φ̃ / h
real GravityForce::wendland_c4_phi(const real r, const real h)
{
#if DIM == 3
    return WendlandC4LookupTable::get_instance().phi_full(r, h);
#else
    // For 1D/2D, fall back to point mass (not implemented)
    return 1.0 / (r + 1e-10);
#endif
}

// Force kernel for Wendland C4: g(r,h) = -dφ̃/dr / r - using lookup table
// For force: F = -G m_i m_j g(r,h) r_ij
real GravityForce::wendland_c4_g(const real r, const real h)
{
#if DIM == 3
    return WendlandC4LookupTable::get_instance().g_full(r, h);
#else
    return 1.0 / (r * r * r + 1e-30);
#endif
}

// =============================================================================
// 1D Kernel-Convolved Gravity Functions
// =============================================================================
// 
// For consistency with SPH, gravity should be computed using the same kernel
// that SPH uses for density estimation. The cumulative kernel function F(q)
// determines how much mass is "to the left" of position x when the mass is
// smoothed over the kernel support.
//
// g_i = -2πG Σ_j m_j [2F((x_i - x_j)/h_j) - 1]
//
// where F(q) = ∫_{-∞}^{q} W(q') dq' is the cumulative distribution
// =============================================================================

// Cumulative distribution function for 1D cubic spline kernel
// W(q) = (2/3)(1 - 1.5q² + 0.75|q|³) for |q| ≤ 1
// W(q) = (1/6)(2-|q|)³ for 1 < |q| ≤ 2
// W(q) = 0 for |q| > 2
//
// F(q) = ∫_{-∞}^{q} W(q') dq' normalized to [0, 1]
real GravityForce::cubic_spline_F_1d(const real q)
{
    if (q <= -2.0) {
        return 0.0;
    } else if (q <= -1.0) {
        // ∫ from -2 to q of (1/6)(2+q')³ dq'
        // = (1/24)(2+q')⁴ evaluated from -2 to q
        // = (1/24)(2+q)⁴
        const real t = 2.0 + q;  // t ranges from 0 to 1
        return (1.0/24.0) * t * t * t * t;
    } else if (q <= 1.0) {
        // ∫ from -2 to -1 of outer + ∫ from -1 to q of inner
        // outer integral = 1/24
        // inner: ∫ (2/3)(1 - 1.5q'² + 0.75q'³) dq' for q' from -1 to q (q' ≥ 0 → use +q'³)
        // For symmetric kernel, split into: q' < 0 and q' > 0
        // Actually for 1D, the kernel is symmetric, so we integrate:
        // (2/3)[q' - 0.5q'³ + 0.1875|q'|⁴] from -1 to q
        
        // At q' = -1: (2/3)[-1 - 0.5(-1) + 0.1875(1)] = (2/3)[-1 + 0.5 + 0.1875] = (2/3)(-0.3125) = -0.2083...
        // At q' = q: (2/3)[q - 0.5q³ + 0.1875q⁴ * sign(q)]
        
        // Simpler: just compute numerically once and fit polynomial
        // For cubic spline in 1D, normalization gives: ∫_{-2}^{2} W dq = 1
        // F(q) = ∫_{-2}^{q} W dq'
        
        // Direct integration:
        // F(q) for -1 ≤ q ≤ 1:
        // = 1/24 + (2/3) ∫_{-1}^{q} (1 - 1.5u² + 0.75|u|³) du
        //
        // For u < 0: |u|³ = -u³
        // For u > 0: |u|³ = +u³
        //
        // Split at u = 0:
        // ∫_{-1}^{0} (1 - 1.5u² - 0.75u³) du = [u - 0.5u³ - 0.1875u⁴]_{-1}^{0}
        //   = 0 - (-1 - 0.5(-1) - 0.1875(1)) = 0 - (-1 + 0.5 - 0.1875) = 0.6875
        //
        // ∫_{0}^{q} (1 - 1.5u² + 0.75u³) du = [u - 0.5u³ + 0.1875u⁴]_{0}^{q}
        //   = q - 0.5q³ + 0.1875q⁴
        //
        // For q < 0: ∫_{-1}^{q} (1 - 1.5u² - 0.75u³) du
        //   = [u - 0.5u³ - 0.1875u⁴]_{-1}^{q} = (q - 0.5q³ - 0.1875q⁴) - (-1 + 0.5 - 0.1875)
        //   = q - 0.5q³ - 0.1875q⁴ + 0.6875
        
        const real q2 = q * q;
        const real q3 = q2 * q;
        const real q4 = q2 * q2;
        
        if (q <= 0.0) {
            // F(q) = 1/24 + (2/3) * (q - 0.5q³ - 0.1875q⁴ + 0.6875)
            return 1.0/24.0 + (2.0/3.0) * (q - 0.5*q3 - 0.1875*q4 + 0.6875);
        } else {
            // F(q) = 1/24 + (2/3) * (0.6875 + q - 0.5q³ + 0.1875q⁴)
            return 1.0/24.0 + (2.0/3.0) * (0.6875 + q - 0.5*q3 + 0.1875*q4);
        }
    } else if (q <= 2.0) {
        // F(q) = F(1) + ∫_{1}^{q} (1/6)(2-q')³ dq'
        // F(1) = 1/24 + (2/3) * (0.6875 + 1 - 0.5 + 0.1875) = 1/24 + (2/3)(1.375) 
        //      = 0.04167 + 0.9167 = 0.9583
        // ∫_{1}^{q} (1/6)(2-q')³ dq' = -(1/24)(2-q')⁴ |_{1}^{q} = -(1/24)[(2-q)⁴ - 1]
        //                            = (1/24)[1 - (2-q)⁴]
        const real F1 = 1.0/24.0 + (2.0/3.0) * 1.375;  // F(1) ≈ 0.9583
        const real t = 2.0 - q;  // t ranges from 1 to 0
        return F1 + (1.0/24.0) * (1.0 - t*t*t*t);
    } else {
        return 1.0;
    }
}

// 1D kernel gravity contribution from particle j at position x_j to particle i at x_i
// Returns contribution to g_i (scalar in 1D)
real GravityForce::kernel_gravity_1d(const real x_ij, const real h)
{
    // x_ij = x_i - x_j
    // g_ij = -2πG m_j [2F(x_ij/h) - 1]
    // This function returns [2F(q) - 1] where q = x_ij / h
    const real q = x_ij / h;
    return 2.0 * cubic_spline_F_1d(q) - 1.0;
}

// =============================================================================
// 2D Kernel-Convolved Gravity Functions
// =============================================================================
// 
// For cylindrical 2D gravity with kernel convolution, we need:
//   F(q) = ∫_0^q W(q') 2πq' dq' / (total mass normalization)
// 
// where W(q) is the 2D cubic spline kernel:
//   W(q) = (10/7π) × w(q)  where w(q) is the radial part
//   w(q) = 1 - 1.5q² + 0.75q³  for 0 ≤ q ≤ 1
//   w(q) = 0.25(2-q)³         for 1 < q ≤ 2
//   w(q) = 0                   for q > 2
//
// The 2D kernel normalization: ∫_0^∞ W(q) 2πq dq = 1
// So F(q) = ∫_0^q W(q') 2πq' dq' represents the enclosed mass fraction.
//
// For gravity force: φ ~ ln(r) in 2D, so dφ/dr ~ 1/r
// Kernel-convolved: g(r) ~ F(q)/r where q = r/h
// =============================================================================

// Cumulative distribution function for 2D cubic spline kernel
// F(q) = ∫_0^q W(q') 2πq' dq' where W is normalized to integrate to 1
// Returns fraction of mass enclosed within radius q*h
real GravityForce::cubic_spline_F_2d(const real q)
{
    if (q <= 0.0) {
        return 0.0;
    } else if (q >= 2.0) {
        return 1.0;
    }
    
    // 2D cubic spline kernel (normalized):
    // W(q) = (10/(7π h²)) × w(q)
    // where w(q) = 1 - 1.5q² + 0.75q³    for 0 ≤ q ≤ 1
    //             = 0.25(2-q)³           for 1 < q ≤ 2
    //
    // Total normalization: ∫_0^2 w(q) 2πq dq = 7π/10 (so W integrates to 1/h²)
    //
    // F(q) = (10/7π) × 2π × ∫_0^q w(q') q' dq' / h²  × h² = (20/7) ∫_0^q w(q') q' dq'
    //
    // For 0 ≤ q ≤ 1:
    // ∫_0^q (1 - 1.5q'² + 0.75q'³) q' dq'
    // = ∫_0^q (q' - 1.5q'³ + 0.75q'⁴) dq'
    // = [q'²/2 - 1.5q'⁴/4 + 0.75q'⁵/5]_0^q
    // = q²/2 - 0.375q⁴ + 0.15q⁵
    //
    // For 1 < q ≤ 2:
    // F(q) = F(1) + (20/7) × ∫_1^q 0.25(2-q')³ q' dq'
    // Let u = 2-q', du = -dq', q' = 2-u
    // ∫_1^q 0.25(2-q')³ q' dq' = -0.25 ∫_{2-1}^{2-q} u³(2-u) du
    //                          = 0.25 ∫_{2-q}^{1} u³(2-u) du
    //                          = 0.25 ∫_{2-q}^{1} (2u³ - u⁴) du
    //                          = 0.25 [u⁴/2 - u⁵/5]_{2-q}^{1}
    //                          = 0.25 [(1/2 - 1/5) - ((2-q)⁴/2 - (2-q)⁵/5)]
    //                          = 0.25 [3/10 - (2-q)⁴/2 + (2-q)⁵/5]
    
    const real q2 = q * q;
    const real q4 = q2 * q2;
    const real q5 = q4 * q;
    
    if (q <= 1.0) {
        // F(q) = (20/7) × (q²/2 - 0.375q⁴ + 0.15q⁵)
        return (20.0/7.0) * (0.5*q2 - 0.375*q4 + 0.15*q5);
    } else {
        // F(1) = (20/7) × (0.5 - 0.375 + 0.15) = (20/7) × 0.275 = 5.5/7
        const real F1 = (20.0/7.0) * 0.275;
        const real t = 2.0 - q;  // t ranges from 1 to 0
        const real t4 = t * t * t * t;
        const real t5 = t4 * t;
        // Contribution from q=1 to q: (20/7) × 0.25 × [3/10 - t⁴/2 + t⁵/5]
        return F1 + (5.0/7.0) * (0.3 - 0.5*t4 + 0.2*t5);
    }
}

// 2D kernel gravity force factor
// For cylindrical 2D gravity: F = -2Gm × g(r,h) × r̂
// Point mass limit: g(r) = 1/r (from φ = -2Gm ln(r), so F = -dφ/dr × r̂ = 2Gm/r × r̂)
// 
// Kernel-convolved: g(r,h) = F(r/h) / r
// where F(q) is the enclosed mass fraction
//
// This ensures:
// - At r >> 2h: F(q) = 1, so g = 1/r (point mass)
// - At r → 0: F(q) ~ q² (smooth), so g ~ q²/r = q/h → 0 (no singularity)
real GravityForce::kernel_gravity_2d(const real r, const real h)
{
    if (r < 1e-30) {
        // At exact center, force is zero
        return 0.0;
    }
    
    const real q = r / h;
    
    if (q >= 2.0) {
        // Point mass limit
        return 1.0 / r;
    }
    
    // Kernel-convolved: F(q) / r
    return cubic_spline_F_2d(q) / r;
}

void GravityForce::initialize(std::shared_ptr<SPHParameters> param)
{
    m_is_valid = param->gravity.is_valid;
    if(m_is_valid) {
        m_constant = param->gravity.constant;
        m_use_fixed_softening = param->gravity.use_fixed_softening;
        m_fixed_softening = param->gravity.fixed_softening;
        m_softening_type = param->gravity.softening_type;
        m_use_kernel_gravity_1d = param->gravity.use_kernel_gravity_1d;
        m_use_kernel_gravity_2d = param->gravity.use_kernel_gravity_2d;
        m_use_kernel_gravity_planar_2d = param->gravity.use_kernel_gravity_planar_2d;
        m_use_kernel_gravity_cylinder_3d = param->gravity.use_kernel_gravity_cylinder_3d;
        
        // Debug output
        std::cerr << "[GRAVITY] Initialized: G=" << m_constant
                  << ", fixed_softening=" << m_use_fixed_softening;
        if (m_use_fixed_softening) {
            std::cerr << " (ε=" << m_fixed_softening << ")";
        }
        std::cerr << ", softening_type=";
        switch (m_softening_type) {
            case GravitySofteningType::HERNQUIST_KATZ:
                std::cerr << "HERNQUIST_KATZ";
                break;
            case GravitySofteningType::WENDLAND_C4:
                std::cerr << "WENDLAND_C4";
                break;
        }
#if DIM == 1
        std::cerr << ", kernel_gravity_1d=" << (m_use_kernel_gravity_1d ? "true" : "false");
#elif DIM == 2
        std::cerr << ", kernel_gravity_2d=" << (m_use_kernel_gravity_2d ? "true" : "false");
        std::cerr << ", kernel_gravity_planar_2d=" << (m_use_kernel_gravity_planar_2d ? "true" : "false");
#elif DIM == 3
        std::cerr << ", kernel_gravity_cylinder_3d=" << (m_use_kernel_gravity_cylinder_3d ? "true" : "false");
#endif
        std::cerr << std::endl;
    }
}

void GravityForce::calculation(std::shared_ptr<Simulation> sim)
{
    if(!m_is_valid) {
        return;
    }

    auto & particles = sim->get_particles();
    const int num = sim->get_particle_num();

#if DIM == 1
    // ==========================================================================
    // 1D Slab Gravity
    // ==========================================================================
    // Two options:
    // 1. Kernel-convolved gravity (default): consistent with SPH pressure
    //    g_i = -πG Σ_j m_j [2F((x_i - x_j)/h_j) - 1]
    //    where F is cumulative kernel function
    //
    //    Note: The factor is πG (not 2πG) because:
    //    - The kernel sum gives M_left - M_right = 2 × M(0→|x|)
    //    - The correct 1D slab formula is g = -2πG × M(0→|x|)
    //    - So we need: g = -2πG × (M_left - M_right)/2 = -πG × (M_left - M_right)
    //
    // 2. Analytical slab (Gauss's law): sharp step function
    //    g(x) = -2πG sign(x) Σ_enclosed
    // ==========================================================================
    
    const real two_pi_G = 2.0 * M_PI * m_constant;
    const real pi_G = M_PI * m_constant;  // For kernel gravity (factor of 2 correction)
    
    if (m_use_kernel_gravity_1d) {
        // Kernel-convolved gravity for consistency with SPH
        #pragma omp parallel for
        for (int i = 0; i < num; ++i) {
            auto & p_i = particles[i];
            const real x_i = p_i.pos[0];
            
            real g_sum = 0.0;
            real phi = 0.0;
            
            for (int j = 0; j < num; ++j) {
                const auto & p_j = particles[j];
                const real x_j = p_j.pos[0];
                const real x_ij = x_i - x_j;
                
                // Use average smoothing length for symmetry
                const real h_ij = 0.5 * (p_i.sml + p_j.sml);
                
                // Kernel gravity contribution
                // g_ij = 2F(q) - 1, where q = x_ij/h
                // F(q) is cumulative kernel, ranges from 0 (far left) to 1 (far right)
                // So 2F-1 ranges from -1 (x_j far right of x_i) to +1 (x_j far left of x_i)
                // The sum gives M_left - M_right = 2 × M(center→|x|)
                g_sum += p_j.mass * kernel_gravity_1d(x_ij, h_ij);
                
                // Potential: use kernel-softened |x| via cumulative integral
                // φ = -2πG Σ_j m_j ∫|x_i - x'| W(x' - x_j) dx'
                // Approximation: use |x_ij| with correction for overlap
                // More accurate: compute double integral. For now, use point mass.
                phi -= two_pi_G * p_j.mass * std::abs(x_ij);
            }
            
            // Use πG (not 2πG) because g_sum = M_left - M_right = 2 × M_enclosed
            p_i.grav_acc[0] = -pi_G * g_sum;
            p_i.phi = phi;
        }
    } else {
        // Original analytical slab gravity (Gauss's law)
        // Create index array sorted by |x|
        std::vector<std::pair<real, int>> sorted_particles(num);
        for (int i = 0; i < num; ++i) {
            sorted_particles[i] = {std::abs(particles[i].pos[0]), i};
        }
        std::sort(sorted_particles.begin(), sorted_particles.end());
        
        // Compute cumulative enclosed mass
        std::vector<real> enclosed_mass(num);
        real cumulative_mass = 0.0;
        
        for (int k = 0; k < num; ++k) {
            int idx = sorted_particles[k].second;
            enclosed_mass[idx] = cumulative_mass;
            cumulative_mass += particles[idx].mass;
        }
        
        #pragma omp parallel for
        for (int i = 0; i < num; ++i) {
            auto & p_i = particles[i];
            const real x_i = p_i.pos[0];
            
            // Gravitational acceleration: g = -2πG sign(x) Σ_enclosed
            real sign_x = (x_i > 0) ? 1.0 : ((x_i < 0) ? -1.0 : 0.0);
            p_i.grav_acc[0] = -two_pi_G * sign_x * enclosed_mass[i];
            
            // Gravitational potential
            real phi = 0.0;
            for (int j = 0; j < num; ++j) {
                const real x_j = particles[j].pos[0];
                phi -= two_pi_G * particles[j].mass * std::abs(x_i - x_j);
            }
            p_i.phi = phi;
        }
    }
    
#else  // DIM == 2 or DIM == 3
    
#if DIM == 2
    // ==========================================================================
    // 2D Gravity Options
    // ==========================================================================
    // 1. Planar 2D slab (useKernelGravityPlanar2D): 1D gravity in y-direction
    //    g_i = -πG Σ_j m_j [2F((y_i - y_j)/h_j) - 1]  (only y-component)
    //
    // 2. Cylindrical 2D disk (useKernelGravity2D): radial gravity in xy-plane
    //    F_i = -2G Σ_j m_j × kernel_gravity_2d(r_ij, h_ij) × r̂_ij
    //
    // 3. Softened point-mass (default fallback)
    // ==========================================================================
    
    if (m_use_kernel_gravity_planar_2d) {
        // ====================================================================
        // 2D PLANAR SLAB: 1D gravity in y-direction only
        // ====================================================================
        // This is for Lane-Emden planar slab in 2D: density varies in y,
        // uniform in x. Gravity acts only in y-direction (like 1D slab).
        const real pi_G = M_PI * m_constant;
        
#ifdef EXHAUSTIVE_SEARCH
        auto * periodic = sim->get_periodic().get();
#endif
        
        #pragma omp parallel for
        for (int i = 0; i < num; ++i) {
            auto & p_i = particles[i];
            const real y_i = p_i.pos[1];  // y-coordinate
            
            real g_sum = 0.0;
            real phi = 0.0;
            
            for (int j = 0; j < num; ++j) {
                const auto & p_j = particles[j];
#ifdef EXHAUSTIVE_SEARCH
                const vec_t r_ij = periodic->calc_r_ij(p_i.pos, p_j.pos);
                const real y_ij = r_ij[1];  // y-component of separation
#else
                const real y_ij = y_i - p_j.pos[1];
#endif
                
                // Average smoothing length
                const real h_ij = 0.5 * (p_i.sml + p_j.sml);
                
                // Kernel gravity contribution (same as 1D)
                g_sum += p_j.mass * kernel_gravity_1d(y_ij, h_ij);
                
                // Potential: φ = -2πG Σ_j m_j |y_i - y_j|
                phi -= 2.0 * M_PI * m_constant * p_j.mass * std::abs(y_ij);
            }
            
            // Gravity acts only in y-direction
            p_i.grav_acc[0] = 0.0;          // No x-acceleration
            p_i.grav_acc[1] = -pi_G * g_sum;  // y-acceleration
            p_i.phi = phi;
        }
    } else if (m_use_kernel_gravity_2d) {
        // Kernel-convolved gravity for 2D cylindrical geometry
        const real two_G = 2.0 * m_constant;  // Factor for 2D logarithmic potential
        
#ifdef EXHAUSTIVE_SEARCH
        auto * periodic = sim->get_periodic().get();
#endif

        #pragma omp parallel for
        for (int i = 0; i < num; ++i) {
            auto & p_i = particles[i];
            const vec_t & r_i = p_i.pos;
            
            real phi = 0.0;
            vec_t force(0.0);
            
            for (int j = 0; j < num; ++j) {
                const auto & p_j = particles[j];
#ifdef EXHAUSTIVE_SEARCH
                const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
#else
                const vec_t r_ij = r_i - p_j.pos;
#endif
                const real r = std::abs(r_ij);
                
                if (r < 1e-30) continue;  // Skip self
                
                // Average smoothing length
                const real h_ij = 0.5 * (p_i.sml + p_j.sml);
                
                // Kernel gravity: g = F(q)/r where F is enclosed mass fraction
                const real g = kernel_gravity_2d(r, h_ij);
                
                // Force: F = -2Gm × g × r̂ = -2Gm × g × r_ij/r
                force -= r_ij * (two_G * p_j.mass * g / r);
                
                // Potential: φ = -2Gm × ln(r) for point mass
                // For kernel: use approximation (not exact, but reasonable)
                phi -= two_G * p_j.mass * std::log(r + 1e-30);
            }
            
            p_i.grav_acc = force;
            p_i.phi = phi;
        }
    } else {
        // Standard softened gravity (fall back to 3D formulas)
#ifdef EXHAUSTIVE_SEARCH
        auto * periodic = sim->get_periodic().get();
#else
        auto * tree = sim->get_tree().get();
#endif

        #pragma omp parallel for
        for (int i = 0; i < num; ++i) {
            auto & p_i = particles[i];
            
#ifdef EXHAUSTIVE_SEARCH
            real phi = 0.0;
            vec_t force(0.0);
            const vec_t & r_i = p_i.pos;

            for (int j = 0; j < num; ++j) {
                const auto & p_j = particles[j];
                const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
                const real r = std::abs(r_ij);
                
                if (m_softening_type == GravitySofteningType::WENDLAND_C4) {
                    if (m_use_fixed_softening) {
                        phi -= m_constant * p_j.mass * wendland_c4_phi(r, m_fixed_softening);
                        force -= r_ij * (m_constant * p_j.mass * wendland_c4_g(r, m_fixed_softening));
                    } else {
                        const real h_ij = 0.5 * (p_i.sml + p_j.sml);
                        phi -= m_constant * p_j.mass * wendland_c4_phi(r, h_ij);
                        force -= r_ij * (m_constant * p_j.mass * wendland_c4_g(r, h_ij));
                    }
                } else {
                    if (m_use_fixed_softening) {
                        const real h_fixed = m_fixed_softening * 2.0;
                        phi -= m_constant * p_j.mass * f(r, h_fixed);
                        force -= r_ij * (m_constant * p_j.mass * g(r, h_fixed));
                    } else {
                        phi -= m_constant * p_j.mass * (f(r, p_i.sml) + f(r, p_j.sml)) * 0.5;
                        force -= r_ij * (m_constant * p_j.mass * (g(r, p_i.sml) + g(r, p_j.sml)) * 0.5);
                    }
                }
            }

            p_i.grav_acc = force;
            p_i.phi = phi;
#else
            tree->tree_force(p_i);
#endif
        }
    }

#else  // DIM == 3
    // ==========================================================================
    // 3D Gravity Options
    // ==========================================================================
    // 1. Cylindrical 3D (useKernelGravityCylinder3D): 2D radial gravity in xy-plane
    //    For infinite cylinder: density varies radially in xy, uniform in z
    //    F_i = -2G Σ_j m_j × g(r_⊥_ij, h_ij) × r̂_⊥_ij  (no z-component)
    //
    // 2. Spherical 3D (default): standard 1/r² point-mass gravity
    // ==========================================================================
    
    if (m_use_kernel_gravity_cylinder_3d) {
        // ====================================================================
        // 3D CYLINDRICAL: 2D radial gravity in xy-plane only
        // ====================================================================
        // This is for Lane-Emden cylinder: density varies radially in xy,
        // uniform in z. Gravity acts only in xy-plane (like 2D disk).
        const real two_G = 2.0 * m_constant;
        
#ifdef EXHAUSTIVE_SEARCH
        auto * periodic = sim->get_periodic().get();
#endif
        
        #pragma omp parallel for
        for (int i = 0; i < num; ++i) {
            auto & p_i = particles[i];
            const real x_i = p_i.pos[0];
            const real y_i = p_i.pos[1];
            
            real phi = 0.0;
            vec_t force(0.0);
            
            for (int j = 0; j < num; ++j) {
                const auto & p_j = particles[j];
#ifdef EXHAUSTIVE_SEARCH
                const vec_t r_ij = periodic->calc_r_ij(p_i.pos, p_j.pos);
                const real x_ij = r_ij[0];
                const real y_ij = r_ij[1];
#else
                const real x_ij = x_i - p_j.pos[0];
                const real y_ij = y_i - p_j.pos[1];
#endif
                
                // Radial distance in xy-plane (perpendicular to cylinder axis)
                const real r_perp = std::sqrt(x_ij * x_ij + y_ij * y_ij);
                
                if (r_perp < 1e-30) continue;  // Skip self
                
                // Average smoothing length
                const real h_ij = 0.5 * (p_i.sml + p_j.sml);
                
                // Kernel gravity: use 2D kernel function
                const real g = kernel_gravity_2d(r_perp, h_ij);
                
                // Force in xy-plane only (no z-component)
                const real force_factor = two_G * p_j.mass * g / r_perp;
                force[0] -= x_ij * force_factor;
                force[1] -= y_ij * force_factor;
                // force[2] = 0 (no z-gravity for cylinder)
                
                // Potential: φ = -2Gm × ln(r_⊥) for point mass
                phi -= two_G * p_j.mass * std::log(r_perp + 1e-30);
            }
            
            p_i.grav_acc = force;
            p_i.phi = phi;
        }
    } else {
        // ====================================================================
        // 3D SPHERICAL: standard point-mass gravity (existing implementation)
        // ====================================================================
    
#ifdef EXHAUSTIVE_SEARCH
    auto * periodic = sim->get_periodic().get();
#else
    auto * tree = sim->get_tree().get();
#endif

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        
#ifdef EXHAUSTIVE_SEARCH
        real phi = 0.0;
        vec_t force(0.0);
        const vec_t & r_i = p_i.pos;

        for(int j = 0; j < num; ++j) {
            const auto & p_j = particles[j];
            const vec_t r_ij = periodic->calc_r_ij(r_i, p_j.pos);
            const real r = std::abs(r_ij);
            
            if (m_softening_type == GravitySofteningType::WENDLAND_C4) {
                // True kernel-convolved gravity using Wendland C4
                if (m_use_fixed_softening) {
                    phi -= m_constant * p_j.mass * wendland_c4_phi(r, m_fixed_softening);
                    force -= r_ij * (m_constant * p_j.mass * wendland_c4_g(r, m_fixed_softening));
                } else {
                    // h-dependent: average over h_i and h_j
                    const real h_ij = 0.5 * (p_i.sml + p_j.sml);
                    phi -= m_constant * p_j.mass * wendland_c4_phi(r, h_ij);
                    force -= r_ij * (m_constant * p_j.mass * wendland_c4_g(r, h_ij));
                }
            } else {
                // Hernquist-Katz (original)
                if (m_use_fixed_softening) {
                    const real h_fixed = m_fixed_softening * 2.0;
                    phi -= m_constant * p_j.mass * f(r, h_fixed);
                    force -= r_ij * (m_constant * p_j.mass * g(r, h_fixed));
                } else {
                    phi -= m_constant * p_j.mass * (f(r, p_i.sml) + f(r, p_j.sml)) * 0.5;
                    force -= r_ij * (m_constant * p_j.mass * (g(r, p_i.sml) + g(r, p_j.sml)) * 0.5);
                }
            }
        }

        p_i.grav_acc = force;
        p_i.phi = phi;
#else
        tree->tree_force(p_i);
#endif
    }
    }  // end else (spherical 3D)
#endif  // DIM == 3
#endif  // DIM == 1
}

}  // namespace sph