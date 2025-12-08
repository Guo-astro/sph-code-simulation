#pragma once

#include "module.hpp"
#include "vector_type.hpp"
#include "parameters.hpp"

namespace sph
{
class SPHParticle;

class GravityForce : public Module {
    bool  m_is_valid;
    real m_constant;
    bool m_use_fixed_softening;
    real m_fixed_softening;
    GravitySofteningType m_softening_type;
    bool m_use_kernel_gravity_1d;  // For 1D: use kernel-convolved gravity
    bool m_use_kernel_gravity_2d;  // For 2D: use kernel-convolved gravity
    bool m_use_kernel_gravity_planar_2d;    // For 2D planar slab: 1D gravity in y-direction
    bool m_use_kernel_gravity_cylinder_3d;  // For 3D cylinder: 2D radial gravity in xy-plane

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
    
    // Gravitational potential kernel derived from Wendland C4
    // Solves ∇²φ̃ = -4πG W(r,h) for Wendland C4 kernel
    static real wendland_c4_phi(const real r, const real h);
    
    // Force kernel: -dφ̃/dr / r (for F = -m_j * g * r_ij)
    static real wendland_c4_g(const real r, const real h);
    
    // 1D kernel-convolved gravity functions
    // Cumulative kernel function F(q) = ∫_{-∞}^{q} W(q') dq' for cubic spline
    static real cubic_spline_F_1d(const real q);
    
    // 1D kernel gravity: g_i = -πG Σ_j m_j [2F((x_i - x_j)/h_j) - 1]
    static real kernel_gravity_1d(const real x_ij, const real h);
    
    // 2D kernel-convolved gravity functions (cylindrical geometry)
    // Cumulative kernel F(q) = ∫_0^q W(q') 2πq' dq' / ∫_0^∞ W(q') 2πq' dq'
    // Represents fraction of mass enclosed within radius q*h
    static real cubic_spline_F_2d(const real q);
    
    // 2D kernel gravity force magnitude: g(r,h) such that F = -2Gm * g(r,h) * r_hat
    // Returns force/mass factor (analogous to 1/r in point-mass limit)
    static real kernel_gravity_2d(const real r, const real h);
};
}
