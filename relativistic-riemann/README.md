# Relativistic Riemann Solver

A Python package for solving relativistic Riemann problems in special relativistic hydrodynamics.

## Overview

This package provides an exact Riemann solver for special relativistic hydrodynamics, based on the algorithm described in Martí and Müller, *J. Fluid Mech.* (1994). It can handle various test cases including shock tubes, blast waves, and high-velocity relativistic flows.

## Features

- **Exact Riemann solver** for special relativistic hydrodynamics
- **Multiple test cases**: SR Sod shock tube, blast waves, relativistic shocks, etc.
- **Visualization tools**: Create animations and static plots
- **Command-line interface**: Easy interactive solving
- **Pure Python**: Built with NumPy and Matplotlib

## Installation

Using `uv` (recommended):

```bash
uv pip install -e .
```

Or with development dependencies:

```bash
uv pip install -e ".[dev]"
```

## Quick Start

### As a Library

```python
from relativistic_riemann import RelativisiticRiemannSolver, test_case_sr_sod

# Get predefined test case
test = test_case_sr_sod()

# Create solver
solver = RelativisiticRiemannSolver(gamma=test['gamma'])
solver.set_initial_states(
    test['pl'], test['rhol'], test['vell'],
    test['pr'], test['rhor'], test['velr']
)

# Solve at t=0.4
x, pressure, density, velocity, internal_energy = solver.solve(t=0.4)

# Save solution
solver.save_solution('solution.dat', x, pressure, density, velocity, internal_energy, t=0.4)
```

### Command-Line Interface

Interactive mode:
```bash
relativistic-riemann
```

### Create Animations

Animate the SR Sod shock tube:
```bash
animate-riemann --test sod --fps 30 --duration 5
```

Create a static comparison plot:
```bash
animate-riemann --test sod --static
```

Available test cases:
- `sod`: SR Sod shock tube
- `blast`: Strong blast wave
- `relativistic`: Relativistic shock with high velocities (0.9c)
- `twoshock`: Two shock collision
- `all`: Run all test cases
# Kitajima Formulation: Baryon Number Density Solver

This implementation follows the formulation presented in:
**Kitajima et al. (2025)** - arXiv:2510.18251v1  
"Baryon number density formulation for special relativistic hydrodynamics"

## Key Differences from Standard Formulation

### Standard (Marti & Mueller 1994)
- Uses **mass density** ρ (lab frame)
- Speed of light **c = 1** (natural units)
- Conservative variables: D = γρ, S = ρhγ²v, E = ρhγ² - P

### Kitajima Formulation
- Uses **baryon number density** N = γn (lab frame)
  - n is rest frame baryon number density
  - Related to mass density via ρ = m_b n (m_b = baryon mass)
- **Explicit speed of light** c (configurable, not set to 1)
- **Canonical variables**:
  - S = γHv (canonical momentum)
  - e = γH - P/(Nc²) (canonical energy)
  - H = 1 + u/c² + P/(nc²) (specific enthalpy)

## Fundamental Equations

Continuity equation (baryon conservation):
```
∂N/∂t + ∇·(Nv) = 0
```

Momentum equation:
```
∂S/∂t + ∇·(Sv) = -∇P/N
```

Energy equation:
```
∂e/∂t + ∇·(ev) = -∇·(Pv)/N
```

Equation of state (ideal gas):
```
P = (γ_c - 1)nu
```
where γ_c is the adiabatic index and u is specific internal energy.

## Usage

### Python Script
```python
from kitajima_solver import KitajimaRiemannSolver

# Create solver with adiabatic index and speed of light
solver = KitajimaRiemannSolver(gamma_c=5.0/3.0, c=1.0)

# Set initial conditions: (P, n, v) for left and right states
solver.set_initial_states(
    Pl=1.0, nl=1.0, vl=0.0,    # Left state
    Pr=0.1, nr=0.125, vr=0.0   # Right state
)

# Solve at time t
x, P, n, N, v, u, gamma, S, e = solver.solve(t=0.4, n_points=400)
```

Returns:
- `x`: Spatial positions
- `P`: Pressure
- `n`: Rest frame baryon number density
- `N`: Lab frame baryon number density (N = γn)
- `v`: Velocity
- `u`: Specific internal energy
- `gamma`: Lorentz factor (γ = 1/√(1 - v²/c²))
- `S`: Canonical momentum (S = γHv)
- `e`: Canonical energy (e = γH - P/(Nc²))

### Animation
```bash
# SR Sod shock tube (c=1, natural units)
python animate_kitajima.py --test sod

# Static comparison plot at multiple times
python animate_kitajima.py --test sod --static

# With explicit speed of light (SI units: m/s)
python animate_kitajima.py --test sod --c 299792458

# Ultra-relativistic case (v ≈ 0.999c)
python animate_kitajima.py --test ultra --fps 30 --duration 5

# All test cases
python animate_kitajima.py --test all
```

## Test Cases

1. **SR Sod** (`--test sod`)
   - Classic special relativistic shock tube
   - P_L=1.0, n_L=1.0, v_L=0 | P_R=0.1, n_R=0.125, v_R=0

2. **Blast Wave** (`--test blast`)
   - Strong pressure jump
   - P_L=1000, n_L=1.0, v_L=0 | P_R=0.01, n_R=1.0, v_R=0

3. **Relativistic Shock** (`--test relativistic`)
   - Moving medium at v=0.9c
   - P_L=1.0, n_L=1.0, v_L=0.9c | P_R=1.0, n_R=1.0, v_R=0

4. **Ultra-Relativistic** (`--test ultra`)
   - Extreme case at v=0.999c
   - P_L=1.0, n_L=1.0, v_L=0.999c | P_R=1.0, n_R=1.0, v_R=0

## Physical Interpretation

### N (Lab Frame Baryon Number Density)
- Represents the baryon number per unit volume as measured in the lab frame
- Contracted by Lorentz factor: N = γn
- Conserved quantity in the continuity equation

### S (Canonical Momentum)
- S = γHv represents the momentum per baryon
- H = 1 + u/c² + P/(nc²) is the specific enthalpy
- More natural variable than standard momentum density

### e (Canonical Energy)
- e = γH - P/(Nc²) is the energy per baryon
- Includes contribution from enthalpy and pressure work
- Dimensionally consistent with explicit c

## Why Use Baryon Number Instead of Mass?

1. **Physical Clarity**: Baryon number is conserved in relativistic flows (no particle creation/annihilation)
2. **Natural Variables**: Canonical variables (S, e) arise naturally from variational principles
3. **Dimensional Analysis**: Explicit c makes units clear and enables SI calculations
4. **Numerical Stability**: Can improve conditioning in some regimes

## Files

- `kitajima_solver.py` - Core solver implementation
- `animate_kitajima.py` - Visualization and test cases
- `kitajima_sod.gif` - Example animation output
- `kitajima_sod_static.png` - Example static plot

## Comparison with Original Solver

Both formulations solve the same physics but with different variables:

| Aspect | Original (Marti & Mueller) | Kitajima |
|--------|---------------------------|----------|
| Primary density | ρ (mass) | n (baryon number) |
| Lab frame density | D = γρ | N = γn |
| Momentum | S = ρhγ²v | S = γHv |
| Energy | E = ρhγ² - P | e = γH - P/(Nc²) |
| Speed of light | c = 1 (implicit) | c explicit |
| EOS | P = (γ-1)ρε | P = (γ-1)nu |

For non-relativistic limits and ideal gases, both give identical physical solutions when properly converted.

## References

Kitajima et al. (2025). "Baryon number density formulation for special relativistic hydrodynamics."  
arXiv:2510.18251v1. https://arxiv.org/abs/2510.18251

Marti, J. M., & Müller, E. (1994). "The analytical solution of the Riemann problem in relativistic hydrodynamics."  
Journal of Fluid Mechanics, 258, 317-333.

## License

MIT License
