# Kernel-Convolved Gravity Implementation

## Overview
The SPH code supports kernel-convolved gravity in all dimensions (1D, 2D, 3D) for consistency with SPH pressure smoothing. Additionally, specialized gravity modes exist for different geometries.

## Key Files
- `include/gravity_force.hpp`: Declares kernel gravity functions
- `src/gravity_force.cpp`: Implements kernel gravity for all dimensions
- `include/parameters.hpp`: Configuration parameters
- `src/solver.cpp`: JSON config parsing
- `tests/test_kernel_gravity.cpp`: BDD-style Google Test suite

## Configuration Parameters
| Parameter | Default | Description |
|-----------|---------|-------------|
| `useKernelGravity1D` | true | Enable 1D kernel gravity for DIM=1 |
| `useKernelGravity2D` | true | Enable 2D radial kernel gravity for DIM=2 |
| `useKernelGravityPlanar2D` | false | Enable 1D planar gravity in y-direction for DIM=2 |
| `useKernelGravityCylinder3D` | false | Enable 2D radial gravity in xy-plane for DIM=3 |

## Geometry Modes

### 1D Slab (DIM=1)
Standard 1D self-gravitating slab geometry.
- Uses: `useKernelGravity1D: true`
- Sample: `polytropic_slab`

### 2D Disk/Cylindrical (DIM=2)
Standard 2D radial gravity (infinite cylinder cross-section).
- Uses: `useKernelGravity2D: true`
- Sample: `lane_emden`

### 2D Planar Slab (DIM=2)
2D simulation with 1D gravity (infinite slab in x, finite in y).
- Uses: `useKernelGravityPlanar2D: true`
- Sample: `polytropic_slab_2d`
- Gravity acts only in y-direction

### 3D Sphere (DIM=3)
Standard 3D spherical gravity.
- Uses: Wendland C4 or Hernquist-Katz softening
- Sample: `lane_emden`

### 3D Cylinder (DIM=3)
3D simulation with 2D radial gravity (infinite cylinder along z).
- Uses: `useKernelGravityCylinder3D: true`
- Sample: `lane_emden_cylinder`
- Gravity acts only in xy-plane (radially)

## Physical Formulas

### 1D (Slab Geometry)
```
g_i = -πG Σ_j m_j [2F((x_i - x_j)/h_ij) - 1]
```
Where F(q) is the cumulative cubic spline kernel function.

### 2D Radial (Cylindrical/Disk Geometry)
```
F_i = -2G Σ_j m_j × kernel_gravity_2d(r_ij, h_ij) × r̂_ij
kernel_gravity_2d(r, h) = F(r/h) / r
```
Where F(q) = ∫₀^q W(q')2πq'dq' is the enclosed mass fraction.

### 2D Planar (in 2D space)
Same formula as 1D slab, but applied to y-coordinate:
```
g_y = -πG Σ_j m_j [2F((y_i - y_j)/h_ij) - 1]
g_x = 0
```

### 3D Cylinder (in 3D space)
Same formula as 2D radial, but applied to xy-plane:
```
F = -2G Σ_j m_j × g(r_⊥, h_ij) × r̂_⊥
g_z = 0
```
Where r_⊥ = √(x² + y²) is the perpendicular distance from cylinder axis.

## Test Suites

Four test suites study energy conservation and grad-h effects:

| Location | Geometry | Dimension | Gravity |
|----------|----------|-----------|---------|
| `simulations/stability/gradh_study_1d/` | Planar slab | 1D | 1D planar |
| `simulations/stability/gradh_study_2d_planar/` | Planar slab | 2D | 1D planar (y-only) |
| `lane_emden/` (2D build) | Disk | 2D | 2D radial |
| `simulations/stability/gradh_study_3d_cylinder/` | Cylinder | 3D | 2D radial (xy-only) |

Each suite compares:
- GSPH with grad-h (stable)
- GSPH without grad-h (core collapse)
- SSPH with grad-h (stable)
- SSPH without grad-h (core collapse)

## Sample Types

### `polytropic_slab` (DIM=1)
1D polytropic slab using planar Lane-Emden equation.

### `polytropic_slab_2d` (DIM=2)
2D planar slab using planar Lane-Emden equation.
- Density varies in y, uniform in x
- Gravity in y-direction only

### `lane_emden` (DIM=2 or DIM=3)
Standard Lane-Emden disk (2D) or sphere (3D).

### `lane_emden_cylinder` (DIM=3)
3D cylinder using cylindrical Lane-Emden equation.
- Density varies radially in xy, uniform in z
- Gravity in xy-plane only

## Lane-Emden Equations

| Geometry | Equation | Surface ξ₁ (n=1.5) |
|----------|----------|-------------------|
| Planar | d²θ/dξ² = -θⁿ | ~2.75 |
| Cylindrical | (1/ξ)d/dξ(ξ dθ/dξ) = -θⁿ | ~2.65 |
| Spherical | (1/ξ²)d/dξ(ξ² dθ/dξ) = -θⁿ | ~3.65 |

## Running Tests

### Using Makefiles
```bash
# 2D Planar Study
cd build && rm -f CMakeCache.txt && cmake -DSPH_DIM=2 .. && make -j4 && cd ..
make gradh_2d_compare_all   # Run all 4 methods + plots + animation

# 3D Cylinder Study  
cd build && rm -f CMakeCache.txt && cmake -DSPH_DIM=3 .. && make -j4 && cd ..
make gradh_3d_compare_all   # Run all 4 methods + plots + animation
```

### Individual Targets
```bash
# 2D Planar
make gradh_2d_help              # Show available targets
make gradh_2d_compare_run       # Run simulations only
make gradh_2d_compare_viz       # Generate comparison plots
make gradh_2d_compare_animate   # Generate GIF animation

# 3D Cylinder
make gradh_3d_help              # Show available targets
make gradh_3d_compare_run       # Run simulations only
make gradh_3d_compare_viz       # Generate comparison plots
make gradh_3d_compare_animate   # Generate GIF animation
```

### Scripts
```
simulations/stability/gradh_study_2d_planar/scripts/
├── compare_methods.py      # Static comparison plots
└── animate_comparison.py   # 2x2 animated GIF

simulations/stability/gradh_study_3d_cylinder/scripts/
├── compare_methods.py      # Static comparison plots
└── animate_comparison.py   # 2x2 animated GIF (xy cross-section)
```

### Output Files
```
simulations/stability/gradh_study_2d_planar/results/comparison/
├── central_density_evolution.png
├── max_density_evolution.png
├── extent_evolution.png
├── velocity_evolution.png
├── energy_evolution.png
├── density_profiles.png
└── 4method_comparison.gif

simulations/stability/gradh_study_3d_cylinder/results/comparison/
├── central_density_evolution.png
├── max_density_evolution.png
├── radius_evolution.png
├── velocity_evolution.png
├── energy_evolution.png
├── density_profiles.png
└── 4method_comparison.gif
```
