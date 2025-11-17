# SRGSPH: Special Relativistic Godunov Smoothed Particle Hydrodynamics

Implementation of the SRGSPH method from:

**Kitajima, K., Inutsuka, S., & Seno, I. (2025)**  
*"Special Relativistic Smoothed Particle Hydrodynamics Based on Riemann Solver"*  
[arXiv:2510.18251v1](https://arxiv.org/abs/2510.18251)

## Overview

SRGSPH is a novel numerical method for special relativistic fluid dynamics that combines:
- **Godunov SPH formulation**: Uses Riemann solver for accurate shock capturing
- **Volume-based particle definition**: Prevents smoothing length discontinuities at contact discontinuities
- **Convolution integrals**: Higher accuracy through proper SPH formulation
- **MUSCL reconstruction**: Second-order spatial accuracy with monotonicity preservation

### Key Features

1. **Special Relativistic Framework**
   - Baryon number density formulation (not mass density)
   - Explicit light speed (not c=1)
   - Canonical momentum and energy variables
   - Proper Lorentz transformation handling

2. **Volume-Based Approach** (Section 2.3-2.4)
   - Particle volume: Vₚ(x) = [Σⱼ W(x - xⱼ, h(x))]⁻¹
   - Smoothing length: h(x) = η Vₚ*(x)^(1/d)
   - Number density: N(x) = ν(x) / Vₚ(x)
   - Prevents artificial smoothing length jumps

3. **Riemann Solver Integration**
   - Exact Riemann solver for particle interactions
   - Minimal artificial viscosity
   - Accurate shock wave treatment

4. **MUSCL Reconstruction** (Section 2.5.2)
   - Second-order spatial accuracy
   - Shock detection and limiting (Eq. 66)
   - Monotonicity preservation

## Installation

```bash
cd /Users/guo/Downloads/sphcode/relativistic-riemann
python -m pip install -e .
```

## Quick Start

### Running Test Problems

The implementation includes three standard test problems from the paper:

1. **Sod Problem** (Section 3.1.1)
```bash
python tests/test_problems.py sod --t_end 0.35
```

2. **Standard Relativistic Blast Wave** (Section 3.1.2)
```bash
python tests/test_problems.py standard_blast --t_end 0.4
```

3. **Strong Relativistic Blast Wave** (Section 3.1.3)
```bash
python tests/test_problems.py strong_blast --t_end 0.16
```

### Visualization

Plot final state:
```bash
python tests/visualize_srgsph.py sod --mode plot
```

Create animation:
```bash
python tests/visualize_srgsph.py sod --mode animate --frames 50 --output sod_animation.gif
```

### Python API

```python
from srgsph.simulator import SRGSPH1D
import numpy as np

# Create simulator
sim = SRGSPH1D(
    gamma_c=5.0/3.0,      # Adiabatic index
    c=1.0,                 # Speed of light
    C_CFL=0.3,            # CFL constant
    use_variable_h=True,   # Variable smoothing length
    use_muscl=True         # MUSCL reconstruction
)

# Setup particles
n_particles = 1000
positions = np.linspace(0, 1, n_particles)
velocities = np.zeros(n_particles)
densities = np.ones(n_particles)
pressures = np.ones(n_particles)

sim.setup_particles(positions, velocities, densities, pressures)

# Run simulation
sim.run(t_end=0.5, output_freq=100)

# Get solution
solution = sim.get_solution()
print(f"Final time: {solution['time']}")
print(f"Particle positions: {solution['x']}")
print(f"Velocities: {solution['v']}")
print(f"Pressures: {solution['P']}")

# Save to file
sim.save_solution('output.dat')
```

## Implementation Details

### Code Structure

```
relativistic-riemann/
├── src/srgsph/
│   ├── __init__.py          # Package initialization
│   ├── particle.py          # Particle class with SR variables
│   ├── kernel.py            # Gaussian and Wendland kernels
│   ├── density.py           # Volume-based density calculation
│   └── simulator.py         # Main SRGSPH simulator
├── tests/
│   ├── test_problems.py     # Standard test problems
│   └── visualize_srgsph.py  # Visualization tools
├── kitajima_solver.py       # Exact Riemann solver
└── README_SRGSPH.md         # This file
```

### Physical Variables

The Particle class stores:

**Conserved quantities** (time-evolved):
- `S`: Canonical momentum per baryon = γHv
- `e`: Canonical energy per baryon = γH - P/(Nc²)

**Primitive variables** (recovered):
- `v`: Velocity
- `n`: Baryon number density (rest frame)
- `P`: Pressure
- `u`: Thermal energy per baryon

**Derived quantities**:
- `γ`: Lorentz factor = 1/√(1 - v²/c²)
- `H`: Enthalpy per baryon = 1 + u/c² + P/(nc²)
- `N`: Lab frame baryon density = γn
- `cs`: Sound speed

**SPH quantities**:
- `h`: Smoothing length
- `Vp`: Particle volume
- `nu`: Baryon number

### Equations of Motion

From Eqs. 64-65 with variable smoothing length:

```
⟨νᵢ Ṡᵢ⟩ = -Σⱼ P*ᵢⱼ V²ᵢⱼ [∇ᵢW(xᵢ-xⱼ, 2hᵢ) - ∇ⱼW(xᵢ-xⱼ, 2hⱼ)]

⟨νᵢ ėᵢ⟩ = -Σⱼ P*ᵢⱼ v*ᵢⱼ · V²ᵢⱼ [∇ᵢW(xᵢ-xⱼ, 2hᵢ) - ∇ⱼW(xᵢ-xⱼ, 2hⱼ)]
```

where P*ᵢⱼ and v*ᵢⱼ are from the Riemann solution.

### Time Integration

Euler method (Eqs. 70-72):

```
⟨νᵢSᵢ⟩ⁿ⁺¹ = ⟨νᵢSᵢ⟩ⁿ + ⟨νᵢṠᵢ⟩Δt
⟨νᵢeᵢ⟩ⁿ⁺¹ = ⟨νᵢeᵢ⟩ⁿ + ⟨νᵢėᵢ⟩Δt
xᵢⁿ⁺¹ = xᵢⁿ + ⟨vᵢ⟩Δt
```

Timestep from CFL condition (Eq. 73):
```
Δt = C_CFL × min_i[hᵢ / cs,ᵢ]
```

## Test Results

The implementation reproduces the test results from Section 3 of the paper:

### Sod Problem (Fig. 4)
- Sharp shock front capture
- Accurate contact discontinuity
- Minimal pressure oscillations

### Blast Wave Problems (Figs. 6-8)
- Correct shock speeds
- Accurate rarefaction waves
- Proper treatment of strong shocks

### Volume-Based Approach Benefits
- Smooth smoothing length profiles (Fig. 3)
- Better handling of contact discontinuities (Fig. 5, 7)
- Reduced overshooting/undershooting

## Configuration Parameters

### Simulation Parameters
- `gamma_c`: Adiabatic index (default: 5/3)
- `c`: Speed of light (default: 1.0 for natural units)
- `C_CFL`: CFL constant for timestep (default: 0.3)

### Volume Calculation
- `eta`: Smoothing length coefficient (default: 1.0)
- `C_smooth`: Smoothing coefficient > 1 for smooth h(x) (default: 2.0)

### MUSCL Reconstruction
- `C_shock`: Shock detection constant (default: 3.0)
- `C_cd`: Contact discontinuity detection (default: 1.0)

### Flags
- `use_variable_h`: Enable variable smoothing length (default: True)
- `use_muscl`: Enable MUSCL reconstruction (default: True)

## Comparison with Paper

This implementation follows the paper equations exactly:
- ✅ Volume-based particle definition (Eqs. 33-37)
- ✅ Variable smoothing length formulation (Eqs. 35-36)
- ✅ SRGSPH force equations (Eqs. 64-65)
- ✅ Primitive variable recovery (Eqs. 67-69)
- ✅ MUSCL with limiters (Eq. 66)
- ✅ Euler time integration (Eqs. 70-72)
- ✅ CFL timestep criterion (Eq. 73)

## Limitations

Current implementation:
- 1D only (2D/3D requires extension)
- Euler time integration (could use RK2/RK4)
- Simple neighbor search (could optimize with tree/cell methods)
- No periodic boundaries (easy to add)
- No external forces (gravity, etc.)

## Future Extensions

Possible improvements following the paper's discussion:
1. 2D/3D extension
2. Higher-order time integration
3. Adaptive particle refinement
4. Kelvin-Helmholtz instability tests (Section 3.3)
5. General relativistic extension (Section 4)

## References

Kitajima, K., Inutsuka, S., & Seno, I. (2025). Special Relativistic Smoothed Particle Hydrodynamics Based on Riemann Solver. *Journal of Computational Physics*, arXiv:2510.18251v1.

Related methods:
- Inutsuka (2002): Godunov SPH
- Pons et al. (2000): Exact Riemann solver for SR hydro
- Rosswog (2010): Previous SR SPH methods

## License

This implementation is provided for research and educational purposes. See LICENSE file for details.

## Contact

For questions or issues, please refer to the original paper or create an issue in the repository.
