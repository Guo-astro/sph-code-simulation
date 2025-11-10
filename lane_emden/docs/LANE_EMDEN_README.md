# Lane-Emden Relaxation Module

## Overview

This module implements density relaxation for Lane-Emden polytropic spheres using analytical forces derived from the exact Lane-Emden equation solution.

## Folder Structure

```
sphcode/
├── data/
│   └── lane_emden/
│       ├── n1.5_3d.dat          # Numerical solution of Lane-Emden ODE
│       ├── n1.5_2d.dat          # 2D solution (for future use)
│       └── solution_*.png       # Diagnostic plots
├── scripts/
│   └── preprocessing/
│       └── generate_lane_emden.py  # Python ODE solver
├── include/
│   └── relaxation/
│       ├── lane_emden_data.hpp      # Data loader and interpolator
│       └── lane_emden_relaxation.hpp # Relaxation force calculator
└── src/
    ├── relaxation/
    │   ├── CMakeLists.txt
    │   ├── lane_emden_data.cpp
    │   └── lane_emden_relaxation.cpp
    └── sample/
        └── lane_emden.cpp           # Initial conditions generator
```

## Workflow

### 1. Generate Lane-Emden Solution (One-time)

```bash
python3 scripts/preprocessing/generate_lane_emden.py
```

This solves the Lane-Emden differential equation:
```
1/ξ² d/dξ(ξ² dθ/dξ) = -θ^n   (3D, n=1.5)
```

Output:
- `data/lane_emden/n1.5_3d.dat` - Solution data (ξ, θ, dθ/dξ)
- Diagnostic plots showing density profile

### 2. Configure Simulation

Edit `sample/lane_emden/lane_emden.json`:

```json
{
  "N": 30,                    // Grid resolution for particle placement
  "useGravity": true,         // Enable self-gravity
  "useRelaxation": true,      // Enable density relaxation
  "relaxationSteps": 100,     // Number of relaxation steps
  "endTime": 5.0,
  "neighborNumber": 50        // Increased for steep density gradients
}
```

### 3. Run Simulation

```bash
./build/sph lane_emden
```

The simulation will:
1. Load Lane-Emden solution from `data/lane_emden/n1.5_3d.dat`
2. Generate equal-mass particles using Fibonacci sphere distribution
3. Apply relaxation forces for first N steps to settle particles
4. Continue with normal SPH evolution

## Physics

### Lane-Emden Equation

For polytrope n=1.5 (γ=5/3):
- **Density**: ρ(r) = ρ_c × θ(ξ)^1.5
- **Pressure**: P(r) = K × ρ^γ
- **Scaling**: r = α × ξ, where α = R/ξ₁

### Hydrostatic Equilibrium

The analytical force from pressure gradient:
```
F_relax = -(K γ n / α) θ^(γ-1) dθ/dξ  (radial direction)
```

This drives particles toward hydrostatic balance.

### Equal-Mass Particle Placement

**Problem**: Uniform grid + 10^7:1 density contrast → neighbor overflow

**Solution**: Place particles at radii such that each encloses equal mass:
1. Compute cumulative mass profile M(r) from Lane-Emden
2. Place particle i at radius where M(r_i) = (i+0.5) × m_particle
3. Use Fibonacci sphere for uniform angular distribution

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| n | 1.5 | Polytropic index |
| γ | 5/3 | Adiabatic index |
| ξ₁ | 3.654 | First zero of θ (dimensionless surface) |
| ρ_c/ρ_s | 10^7 | Central/surface density contrast |
| M_total | 1.0 | Total mass |
| R | 1.0 | Sphere radius |

## Modules

### LaneEmdenData
**File**: `include/relaxation/lane_emden_data.hpp`

Loads and interpolates the numerical solution:
- `load_from_file()` - Read data file
- `get_theta(ξ)` - Interpolate θ at any ξ
- `dtheta_dxi(ξ)` - Interpolate dθ/dξ

### LaneEmdenRelaxation
**File**: `include/relaxation/lane_emden_relaxation.hpp`

Computes and applies relaxation forces:
- `initialize()` - Load parameters and data
- `compute_relaxation_force()` - Calculate force for one particle
- `apply_relaxation()` - Apply to all particles with velocity damping

## Usage Example

```cpp
#include "relaxation/lane_emden_relaxation.hpp"

// In solver initialization
LaneEmdenRelaxationParams params;
params.alpha_scaling = R / xi_1;
params.rho_center = M / (4π α³ ξ₁² |dθ/dξ|);
params.K = /* polytropic constant */;
params.gamma = 5.0/3.0;

auto relax = std::make_shared<LaneEmdenRelaxation>();
relax->initialize(params);

// During first N timesteps
if (step < relaxation_steps) {
    relax->apply_relaxation(sim, damping_factor=0.9);
}
```

## Expected Results

With proper relaxation:
- ✓ Particles settle into equilibrium within ~100 steps
- ✓ 90% of particles remain at r < R
- ✓ Energy oscillations < 1%
- ✓ No neighbor overflow (with equal-mass placement)

## Troubleshooting

**Neighbor overflow**: Increase `neighbor_list_size` in `defines.hpp`

**Particles escape**: 
- Increase `relaxationSteps`
- Check that `useGravity = true`
- Verify Lane-Emden data file exists

**Simulation crashes**:
- Check that `data/lane_emden/n1.5_3d.dat` exists
- Run Python script to generate solution
- Verify file format (header + data columns)

## References

- Chandrasekhar, S. (1939). *An Introduction to the Study of Stellar Structure*
- Lane-Emden equation: Standard polytrope model for self-gravitating spheres
- Equal-mass SPH: Prevents resolution bias in steep density gradients
