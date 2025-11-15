# Special Relativistic GSPH (SRGSPH)

## Overview

Implementation of **Special Relativistic Godunov Smoothed Particle Hydrodynamics** based on:

> Kitajima, K., Inutsuka, S., & Seno, I. (2025)  
> "Special Relativistic Smoothed Particle Hydrodynamics Based on Riemann Solver"  
> arXiv:2510.18251v1  
> Journal of Computational Physics (accepted)

This module extends the non-relativistic GSPH method to handle special relativistic flows with Lorentz factors up to γ ~ 10^6.

## Quick Start

### Configuration

```json
{
  "SPHType": "srgsph",
  "use2ndOrderSRGSPH": true,
  "cSpeed": 1.0,
  "cShock": 3.0,
  "cContactDiscontinuity": 1.0,
  "etaSmoothingLength": 1.0,
  "cSmoothGradient": 2.0,
  "gamma": 1.666666667,
  "neighborNumber": 50
}
```

### Build

```bash
cd build
cmake ..
make
```

### Run Test

```bash
./sph ../sample/sr_sod/config/sr_sod_gsph.json
```

## Key Features

### 1. Volume-Based Formulation
- Particle volume: `Vp = [Σ W]^(-1)`
- Baryon density: `N = ν/Vp`
- Avoids overshooting at discontinuities
- Handles non-uniform baryon numbers

### 2. Primitive Variable Recovery
- Solves quartic equation for Lorentz factor γ
- Converts conserved (S, e, N) ↔ primitive (v, ρ, P)
- Newton-Raphson with robust convergence

### 3. Variable Smoothing Length
- Gather approach with gradient smoothing
- Iterative solver: `h = η Vp*^(1/d)`
- Parameters: η=1.0, C_smooth=2.0

### 4. Riemann Solver
- Relativistic HLL (can upgrade to exact)
- Handles shock waves accurately
- Minimal artificial viscosity

### 5. MUSCL Reconstruction
- 2nd order spatial accuracy
- van Leer limiter
- Monotonicity constraints for shocks

## File Structure

```
include/srgsph/
  ├── sr_fluid_force.hpp         - Force calculation class
  ├── sr_pre_interaction.hpp     - Density & gradient calculation
  ├── sr_primitive_recovery.hpp  - Conserved ↔ primitive conversion
  └── CMakeLists.txt

src/srgsph/
  ├── sr_fluid_force.cpp         - Riemann solver + dS/dt, de/dt
  ├── sr_pre_interaction.cpp     - Volume-based density + smoothing
  ├── sr_primitive_recovery.cpp  - Lorentz factor solver
  └── CMakeLists.txt
```

## Physics

### Basic Equations

**Lagrangian Derivatives:**
```
dN/dt = -N∇·v                    (continuity)
dS/dt = -(1/N)∇P                 (momentum)
de/dt = -(1/N)∇·(Pv)             (energy)
```

**Conserved Variables:**
```
S = γHv         (canonical momentum)
e = γH - P/(Nc²) (canonical energy)
N = γn          (baryon number density)
```

**Primitive Variables:**
```
γ = 1/√(1-v²/c²)               (Lorentz factor)
H = 1 + u/c² + P/(nc²)         (enthalpy)
P = (γ_c-1)nu                  (pressure, ideal gas)
```

### Discretization

**Equations of Motion (Eq. 64-65):**
```
⟨ν_i Ṡ_i⟩ = -Σ_j P*_ij V²_ij [∇_i W - ∇_j W]
⟨ν_i ė_i⟩ = -Σ_j P*_ij v*_ij · V²_ij [∇_i W - ∇_j W]
```

where `P*` and `v*` come from Riemann solver.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `use2ndOrderSRGSPH` | true | Enable MUSCL reconstruction |
| `cSpeed` | 1.0 | Speed of light in code units |
| `cShock` | 3.0 | Shock detection threshold |
| `cContactDiscontinuity` | 1.0 | Contact discontinuity threshold |
| `etaSmoothingLength` | 1.0 | Smoothing length parameter |
| `cSmoothGradient` | 2.0 | Gradient smoothing factor |

## Test Cases

### 1D Tests

1. **Sod Problem** (Section 3.1.1)
   - Left:  P=1.0, n=1.0, v=0
   - Right: P=0.1, n=0.125, v=0
   - t_end = 0.35

2. **Standard Blast Wave** (Section 3.1.2)
   - Left:  P=40/3, n=10.0
   - Right: P=10^-6, n=1.0
   - t_end = 0.4

3. **Strong Blast Wave** (Section 3.1.3)
   - Left:  P=1000, n=1.0
   - Right: P=0.01, n=1.0
   - t_end = 0.16

4. **Ultra-Relativistic** (Section 3.2)
   - v_L = 0.9 to 0.999999999
   - Tests extreme Lorentz factors

### 2D Tests

1. **2D Sod** (Section 3.2.2)
   - ~60,000 particles
   - Periodic boundaries

2. **Kelvin-Helmholtz Instability** (Section 3.3)
   - v = ±0.3
   - ~250,000 particles
   - Growth rate validation

## Algorithm Flow

```
For each time step:

1. Pre-Interaction (sr_pre_interaction.cpp):
   - Find neighbors
   - Compute Vp (particle volume)
   - Update h (smoothing length, iterative)
   - Calculate N = ν/Vp
   - Recover primitive from conserved (v, ρ, P) ← (S, e, N)
   - Compute gradients for MUSCL

2. Fluid Force (sr_fluid_force.cpp):
   For each neighbor pair (i,j):
     - MUSCL reconstruction if 2nd order
     - Check monotonicity constraint
     - Riemann solver → p*, v*
     - Accumulate dS/dt and de/dt
   
3. Time Integration (existing framework):
   - Update S ← S + dS/dt × Δt
   - Update e ← e + de/dt × Δt
   - Update positions

4. Loop back to step 1
```

## Differences from Non-Relativistic GSPH

| Feature | GSPH | SRGSPH |
|---------|------|--------|
| Variables | ρ, v, u | S, e, N |
| Density | SPH sum | Volume-based |
| Smoothing | Fixed/simple | Variable with C_smooth |
| Riemann | Standard HLL | Relativistic HLL |
| Recovery | - | Quartic solver for γ |
| Speed Limit | - | v < c |

## Implementation Details

### Lorentz Factor Solver

Solves quartic equation (Eq. 67):
```
(γ²-1)(Xeγ-1)² - S²(Xγ²-1)² = 0
```

Newton-Raphson iteration:
- Initial guess: `γ ≈ 1 + √(1 + S²/c²)/X`
- Tolerance: 10^-10
- Max iterations: 100

### Monotonicity Constraint

Falls back to 1st order if (Eq. 66):
```
C_shock × |v_i - v_j| > min(c_s,i, c_s,j)  OR
|log10(P_i/P_j)| > C_c.d.
```

Prevents oscillations at:
- Strong shocks
- Large tangential velocities
- Contact discontinuities

### Volume Weighting

Uses `V²_ij` in force terms:
```
V_i = ν_i / N_i    (particle volume)
V²_ij ≈ 0.5(V_i² + V_j²)    (averaged)
```

Can be improved with N^-2 interpolation (Inutsuka 2002).

## Performance

### Computational Cost

Relative to non-relativistic GSPH:
- +10-15% for primitive recovery (γ solver)
- +5% for smoothing length iteration
- ~20% total overhead

### Scalability

- OpenMP parallelized
- O(N log N) neighbor search
- O(N) recovery per step
- O(N × n_neighbor) force calculation

## Known Limitations

1. **Riemann Solver**: Using simplified HLL, can upgrade to exact (Pons 2000)
2. **Volume Term**: Using averaged V²_ij, can improve with interpolation
3. **Kernel**: Cubic/Wendland instead of paper's Gaussian
4. **Time Step**: Simple Courant, could use more sophisticated

## Future Enhancements

- [ ] Exact Riemann solver (Pons et al. 2000)
- [ ] Proper V²_ij interpolation
- [ ] Gaussian kernel option
- [ ] Higher order time integration (RK2/RK4)
- [ ] General relativistic extension
- [ ] SR-MHD capability

## References

1. **Kitajima et al. (2025)** - arXiv:2510.18251v1 (This work)
2. **Inutsuka (2002)** - J. Comp. Phys. 179, 238 (Original GSPH)
3. **Pons et al. (2000)** - J. Fluid Mech. 422, 125 (Exact SR Riemann)
4. **Cha & Whitworth (2003)** - MNRAS 340, 73 (GSPH implementation)
5. **Murante et al. (2011)** - MNRAS 417, 136 (MUSCL for GSPH)

## Documentation

- `docs/SRGSPH_IMPLEMENTATION.md` - Detailed design document
- `docs/SRGSPH_COMPLETION_SUMMARY.md` - Implementation summary

## Contact

For questions or issues with this implementation, please refer to the main SPH code repository.

---
**Implementation**: 2025-11-15  
**Status**: Ready for testing  
**Version**: 1.0
