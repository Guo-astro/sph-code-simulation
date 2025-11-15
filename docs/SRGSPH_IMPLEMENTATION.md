# Special Relativistic GSPH Implementation Plan

## Overview

Implementation of **Special Relativistic Godunov SPH (SRGSPH)** based on:
> Kitajima, K., Inutsuka, S., & Seno, I. (2025). "Special Relativistic Smoothed Particle Hydrodynamics Based on Riemann Solver"  
> arXiv:2510.18251v1

This extends the existing non-relativistic GSPH method to special relativistic hydrodynamics.

## Key Differences from Non-Relativistic GSPH

### 1. Basic Equations (Section 2.1)

**Non-relativistic GSPH:**
- Evolves density ρ, velocity v, energy u
- EOS: P = (γ-1)ρu

**Special Relativistic GSPH:**
- Evolves baryon number density N, canonical momentum S, canonical energy e
- Lagrangian derivatives:
  ```
  dN/dt = -N∇·v                           (Eq. 1)
  dS/dt = -(1/N)∇P                        (Eq. 2)
  de/dt = -(1/N)∇·(Pv)                    (Eq. 3)
  ```

Where:
- **Lorentz factor**: γ = 1/√(1-v²/c²)
- **Canonical momentum**: S = γHv  (Eq. 5)
- **Canonical energy**: e = γH - P/(Nc²)  (Eq. 6)
- **Enthalpy per baryon**: H = 1 + u/c² + P/(nc²)  (Eq. 8)
- **Rest frame density**: n (baryon number density in rest frame)
- **Lab frame density**: N = γn

### 2. Primitive Variable Recovery (Section 2.6)

**Critical difference**: In SR-GSPH, we evolve **conserved variables** (S, e) but need **primitive variables** (v, ρ, P) for the Riemann solver.

Recovery algorithm:
1. Solve quartic equation for γ (Eq. 67):
   ```
   (γ²-1)(Xe⟨γ⟩-1)² - S²(Xγ²-1)² = 0
   ```
   where X = γ_c/(γ_c-1)

2. Recover velocity from γ (Eq. 69):
   ```
   v = (Xγ²-1)/(γ(Xe⟨γ⟩-1)) S
   ```

3. Calculate pressure and other thermodynamic quantities from EOS

### 3. Volume-Based Particle Approach (Sections 2.3-2.4)

**Key innovation**: Instead of defining number density directly, define **particle volume**:

```
Vp(x) = [Σ_j W(x-x_j, h(x))]^(-1)     (Eq. 33)
```

Then number density:
```
N(x) = ν(x) / Vp(x)                     (Eq. 37)
```

**Advantages** (Section 2.4):
- Smoothing length varies smoothly across discontinuities
- Avoids overshooting/undershooting with non-uniform baryon numbers
- Better handling of low-density regions

### 4. Variable Smoothing Length (Section 2.3)

**Gather approach** with smoothing:
```
h(x) = η Vp*^(1/d)                      (Eq. 35)
Vp*(x) = [Σ_j W(x-x_j, C_smooth h(x))]^(-1)   (Eq. 36)
```

Parameters (from paper):
- η = 1.0
- C_smooth = 2.0 (smooths gradient ∇h, allows approximation ∇h ≈ 0)

### 5. Equations of Motion and Energy (Section 2.5)

**Fixed smoothing length** (Eqs. 31-32):
```
⟨Ṡ_i⟩ = -Σ_j ν P*_ij V²_ij,interp [∇_i W(x_i-x_j, 2h) - ∇_j W(x_i-x_j, 2h)]

⟨ė_i⟩ = -Σ_j ν P*_ij v*_ij · V²_ij,interp [∇_i W(x_i-x_j, 2h) - ∇_j W(x_i-x_j, 2h)]
```

**Variable smoothing length** (Eqs. 64-65):
```
⟨ν_i Ṡ_i⟩ = -Σ_j P*_ij V²_ij,interp [∇_i W(x_i-x_j, 2h_i) - ∇_j W(x_i-x_j, 2h_j)]

⟨ν_i ė_i⟩ = -Σ_j P*_ij v*_ij · V²_ij,interp [∇_i W(x_i-x_j, 2h_i) - ∇_j W(x_i-x_j, 2h_j)]
```

Key: P*_ij and v*_ij come from **Riemann solver** between particles i and j

### 6. MUSCL Reconstruction (Section 2.5.2)

**Second-order spatial accuracy** using van Leer limiter.

**Monotonicity constraint** (Eq. 66) - critical for stability:
```
Set gradients to zero if:
  C_shock e_ij·(v_i - v_j) > min(c_s,i, c_s,j)  OR
  |log10(P_i/P_j)| > C_c.d.
```

Where:
- C_shock = 3 (shock detection threshold)
- C_c.d. = 1 (contact discontinuity threshold)

This prevents overshooting at shocks and large tangential velocities.

## Implementation Architecture

### Directory Structure

Following project coding rules:
```
include/srgsph/
  sr_fluid_force.hpp
  sr_pre_interaction.hpp
  sr_primitive_recovery.hpp
  CMakeLists.txt

src/srgsph/
  sr_fluid_force.cpp
  sr_pre_interaction.cpp
  sr_primitive_recovery.cpp
  CMakeLists.txt
```

### Class Design

#### 1. SRGSPHFluidForce

Inherits from `sph::FluidForce`

**Members:**
```cpp
class FluidForce : public sph::FluidForce {
    bool m_is_2nd_order;
    real m_gamma_eos;  // EOS gamma
    real m_c_speed;    // speed of light
    real m_c_shock;    // shock detection parameter
    real m_c_cd;       // contact discontinuity parameter
    
    // Riemann solver: (left[4], right[4]) -> (pstar, vstar)
    std::function<void(const real[], const real[], real&, real&)> m_solver;
    
    void exact_riemann_solver();  // Pons et al. (2000)
    
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};
```

**Key functions:**
- `calculation()`: Main loop computing dS/dt and de/dt
- `exact_riemann_solver()`: SR Riemann solver
- MUSCL reconstruction with SR monotonicity constraints

#### 2. SRGSPHPreInteraction

Inherits from `sph::PreInteraction`

**Members:**
```cpp
class PreInteraction : public sph::PreInteraction {
    real m_eta;        // smoothing length parameter (1.0)
    real m_c_smooth;   // smoothing parameter (2.0)
    real m_c_speed;    // speed of light
    
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
    
private:
    real compute_volume(/* ... */);
    real compute_smoothing_length(/* ... */);
};
```

**Key functions:**
- `calculation()`: Compute Vp, h, N, gradients
- `compute_volume()`: Implements Eq. 33
- `compute_smoothing_length()`: Implements Eqs. 35-36 with iteration

#### 3. SRPrimitiveRecovery

Utility class for converting conserved ↔ primitive variables

**Functions:**
```cpp
namespace srgsph {
    
// Solve quartic for gamma (Eq. 67)
real solve_lorentz_factor(
    const vec_t& S,  // canonical momentum
    real e,          // canonical energy
    real N,          // baryon number density
    real gamma_eos   // EOS gamma
);

// Recover velocity from gamma (Eq. 69)
vec_t recover_velocity(
    const vec_t& S,
    real gamma,
    real e,
    real gamma_eos
);

// Compute thermodynamic quantities
struct PrimitiveVariables {
    vec_t vel;
    real density;
    real pressure;
    real sound_speed;
    real enthalpy;
};

PrimitiveVariables conserved_to_primitive(
    const vec_t& S,
    real e,
    real N,
    real gamma_eos,
    real c_light
);

}
```

### Particle Data Extensions

Two options:

**Option A: Extend SPHParticle**
```cpp
// Add to particle.hpp
class SPHParticle {
    // Existing members...
    
    // SR-specific (only allocated for SRGSPH)
    vec_t S;        // canonical momentum
    real e;         // canonical energy  
    real N;         // baryon number density
    real gamma_lor; // Lorentz factor
    real enthalpy;  // H
    real nu;        // baryon number in particle
};
```

**Option B: Create derived class**
```cpp
// include/srgsph/sr_particle.hpp
class SRParticle : public SPHParticle {
    vec_t S;        // canonical momentum
    real e;         // canonical energy
    real N;         // baryon number density
    real gamma_lor; // Lorentz factor
    real enthalpy;  // H
    real nu;        // baryon number in particle
};
```

**Recommendation**: Option A (extend SPHParticle) is simpler and follows existing architecture.

### Parameter Updates

**include/parameters.hpp:**
```cpp
enum struct SPHType {
    SSPH,
    DISPH,
    GSPH,
    GDISPH,
    SRGSPH,  // Add this
};

struct SPHParameters {
    // ... existing members ...
    
    struct SRGSPH {
        bool is_2nd_order;
        real c_speed;       // speed of light (default: 1.0)
        real c_shock;       // shock detection (default: 3.0)
        real c_cd;          // contact discontinuity (default: 1.0)
        real eta;           // smoothing parameter (default: 1.0)
        real c_smooth;      // smoothing length smoother (default: 2.0)
    } srgsph;
};
```

## Test Cases (Section 3 of Paper)

### 1D Tests

All with γ_c = 5/3, c = 1:

1. **Sod Problem** (Section 3.1.1)
   - Left:  (P, n, vx) = (1.0, 1.0, 0)
   - Right: (P, n, vx) = (0.1, 0.125, 0)
   - End time: t = 0.35
   - Tests: Equal baryon numbers, different baryon numbers

2. **Standard Blast Wave** (Section 3.1.2)
   - Left:  (P, n, vx) = (40/3, 10.0, 0)
   - Right: (P, n, vx) = (10^-6, 1.0, 0)
   - End time: t = 0.4

3. **Strong Blast Wave** (Section 3.1.3)
   - Left:  (P, n, vx) = (1000, 1.0, 0)
   - Right: (P, n, vx) = (0.01, 1.0, 0)
   - End time: t = 0.16
   - Test C_smooth sensitivity

4. **Ultra-Relativistic Shock** (Section 3.2)
   - Left:  (P, n, vx) = (1.0, 1.0, 0.9-0.999999999)
   - Right: (P, n, vx) = (1.0, 1.0, 0)
   - High Lorentz factor test

5. **Tangential Velocity Test** (Section 3.2.1)
   - Includes tangential velocity components
   - Tests 9 combinations: vt ∈ {0, 0.9, 0.99}

### 2D Tests

1. **2D Sod Problem** (Section 3.2.2)
   - ~60,000 particles
   - Periodic boundaries in y-direction
   - Triangular lattice

2. **Kelvin-Helmholtz Instability** (Section 3.3)
   - Two layers: vx = ±0.3
   - Initial perturbation: vy = A0 sin(2πx/λ)
   - ~250,000 particles
   - Tests growth rate vs. Bodo et al. (2004)

## Implementation Roadmap

### Phase 1: Core SR Infrastructure ✅
1. ✅ Create directory structure
2. ✅ Extend SPHParticle with SR variables
3. ✅ Implement primitive variable recovery
4. ✅ Add SRGSPH to parameters and solver

### Phase 2: Pre-Interaction Module
1. Implement volume-based particle density
2. Implement variable smoothing length calculation
3. Implement gradient calculation for MUSCL
4. Test with simple configurations

### Phase 3: Fluid Force Module
1. Implement SR Riemann solver (exact solver)
2. Implement equations of motion and energy
3. Implement MUSCL reconstruction
4. Implement monotonicity constraints
5. Test with 1D Sod problem

### Phase 4: Validation & Testing
1. 1D Sod problem (equal baryon numbers)
2. 1D Sod problem (different baryon numbers)
3. Standard blast wave
4. Strong blast wave
5. Ultra-relativistic shock
6. Tangential velocity tests
7. 2D Sod problem
8. Kelvin-Helmholtz instability

### Phase 5: Optimization & Documentation
1. Performance optimization
2. Documentation
3. Sample configurations
4. Comparison scripts

## Key Implementation Notes

### 1. Speed of Light Units

Paper uses c = 1 in code units. Need to handle conversion:
```cpp
// Internal: use c = 1
// I/O: may need physical units
```

### 2. Time Step (Section 3)

Paper uses simple Courant condition:
```cpp
Δt = C_CFL * min_i[h_i / c_s,i]
```
with C_CFL = 0.3

### 3. Kernel Function

Paper uses **Gaussian kernel** (Eq. 11):
```
W(x, h) = (1/(h√π))^d exp(-x²/h²)
```

But our code uses Cubic Spline / Wendland. Options:
- Add Gaussian kernel to kernel functions
- Test with existing kernels first

### 4. Volume Interpolation (Section 2.2.3)

The volume term V²_ij,interp requires interpolation of N^(-2)(x).

Paper references Inutsuka (2002) for details - may need to study original GSPH paper.

### 5. Baryon Number Conservation

Unlike non-relativistic SPH, baryon number ν in each particle can vary.
Need to track and conserve total baryon number.

## References

### Primary
- Kitajima, Inutsuka, Seno (2025) - arXiv:2510.18251v1 - **This paper**
- Inutsuka (2002) - J. Comp. Phys. 179, 238 - **Original GSPH**
- Pons et al. (2000) - J. Fluid Mech. 422, 125 - **Exact SR Riemann solver**

### Supporting
- Cha & Whitworth (2003) - MNRAS 340, 73 - **GSPH implementation**
- Murante et al. (2011) - MNRAS 417, 136 - **MUSCL for GSPH**
- van Leer (1979) - J. Comp. Phys. 32, 101 - **Limiter**
- Bodo et al. (2004) - Phys. Rev. E 70 - **KHI growth rate**

### In Our Codebase
- `src/gsph/` - Non-relativistic GSPH reference
- `docs/GDISPH_IMPLEMENTATION.md` - Recent similar implementation
- `sample/sr_sod/` - Existing SR test case (if any)

## Next Steps

1. ✅ Complete this design document
2. Review with team/supervisor
3. Start Phase 1 implementation
4. Iterate based on test results

---
**Document Status**: Draft v1.0  
**Author**: GitHub Copilot  
**Date**: 2025-11-15
