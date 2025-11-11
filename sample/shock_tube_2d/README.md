# 2D Shock Tube Test

A 2D benchmark test for compressible SPH solvers based on the classic Riemann problem, extended to two dimensions. This test validates shock-capturing schemes with analytical solutions.

## Overview

The 2D shock tube extends the classic Sod shock tube to two dimensions, testing SPH methods' ability to handle strong shocks, contact discontinuities, and rarefaction waves in a planar geometry. The analytical solution is known and provides an exact benchmark.

**Key Features:**
- Strong shock propagation in 2D
- Exact analytical solution (1D + planar symmetry)
- Tests artificial viscosity and shock capturing
- Contact discontinuity tracking
- Rarefaction wave structure

## Physics

**Initial Conditions:**
- **Left state** (x < 0.5): ρ₁ = 1.0, P₁ = 1.0, u₁ = 0
- **Right state** (x ≥ 0.5): ρ₂ = 0.125, P₂ = 0.1, u₂ = 0
- **Adiabatic index**: γ = 1.4 (diatomic gas)
- **Planar symmetry**: No y-dependence

**Expected Wave Structure:**
1. **Leftward rarefaction fan** - smooth expansion
2. **Contact discontinuity** - density jump, continuous pressure/velocity
3. **Rightward shock** - sharp density, pressure, and velocity jump

**Analytical Solution:**
Exact Riemann solver provides profiles for:
- Density ρ(x,t)
- Velocity u(x,t)
- Pressure P(x,t)
- Internal energy e(x,t)

Solution is independent of y-coordinate (planar symmetry).

## Configuration

**Domain**: [0, 1] × [0, 0.2], periodic in y, reflective in x
**Resolution**: 200 × 40 particles (8,000 total)
**End time**: t = 0.2
**Output interval**: Δt = 0.01 (21 snapshots)
**Gamma**: γ = 1.4 (ideal gas)
**Kernels**: Wendland C4 or Cubic Spline

## Quick Start

### Single Method Runs
```bash
make shock_tube_2d_run                  # Default run (GSPH + Wendland)
make shock_tube_2d_gsph_wendland        # GSPH with Wendland kernel
make shock_tube_2d_gsph_cubic           # GSPH with Cubic Spline kernel
make shock_tube_2d_ssph_wendland        # SSPH with Wendland kernel
```

### Multi-Method Comparison (Wendland kernel)
```bash
make shock_tube_2d_compare_run          # Run all 5 methods: GSPH, SSPH, DISPH, GDISPH, GDISPH+Balsara
make shock_tube_2d_compare_viz          # Generate analytical overlays (1D cuts)
make shock_tube_2d_compare_2d           # Generate 2D spatial plots
make shock_tube_2d_compare_animate      # Create comparison animations (1D + 2D)
make shock_tube_2d_compare_all          # Run all + visualizations
```

### Kernel Comparison (Wendland vs Cubic Spline)
```bash
make shock_tube_2d_kernel_run           # Run all 10 simulations (5 methods × 2 kernels)
make shock_tube_2d_kernel_viz           # Generate analytical overlays for all configs
make shock_tube_2d_kernel_animate       # Create kernel comparison animations
make shock_tube_2d_kernel_all           # Run all + visualizations
```

### Cleanup
```bash
make shock_tube_2d_clean                # Remove all results
make shock_tube_2d_compare_clean        # Remove comparison results only
make shock_tube_2d_kernel_clean         # Remove kernel comparison results only
```

## Output Structure

```
sample/shock_tube_2d/results/
├── gsph_wendland/
│   ├── snapshot_*.csv         # Simulation data
│   ├── energy.dat             # Energy conservation history
│   ├── plots/                 # 1D analytical comparison plots
│   │   └── snapshot_*_comparison.png  # 21 analytical overlays (y-averaged)
│   └── plots_2d/              # 2D spatial visualizations
│       └── snapshot_*_2d.png  # 21 spatial plots
├── ssph_wendland/
│   ├── plots/                 # 1D analytical overlays
│   └── plots_2d/              # 2D spatial plots
├── disph_wendland/
├── gdisph_wendland/
├── gdisph_balsara_wendland/
├── gsph_cubic/                # (kernel comparison only)
├── ssph_cubic/                # (kernel comparison only)
├── disph_cubic/               # (kernel comparison only)
├── gdisph_cubic/              # (kernel comparison only)
├── gdisph_balsara_cubic/      # (kernel comparison only)
├── comparison/
│   ├── comparison_1d_*.png    # Multi-method 1D comparisons
│   └── comparison_animation.gif
├── animations/
│   ├── gsph_wendland_1d.gif   # 1D evolution with analytical overlay
│   ├── gsph_wendland_2d.gif   # 2D spatial evolution
│   └── ...
└── kernel_comparison/
    ├── kernel_comparison_*.png
    └── kernel_comparison.gif
```

## Visualization Types

### 1D Analytical Comparison Plots
Four-panel plots showing y-averaged quantities vs analytical solution:
1. **Density ρ(x)** - Shows shock, contact, and rarefaction
2. **Velocity u(x)** - Characteristic velocity profile
3. **Pressure P(x)** - Pressure jumps across waves
4. **Internal energy e(x)** - Thermal energy distribution

Each plot shows:
- Blue dots: SPH simulation (y-averaged)
- Red line: Exact analytical solution

### 2D Spatial Plots
Full 2D visualization of particle distribution:
- **Density field** - Color-coded particle plot
- **Velocity field** - Vector arrows showing flow
- **Pressure field** - Color-coded spatial distribution
- Shows planar symmetry and any numerical artifacts

### Animations
- **1D animations**: Evolution of 1D profiles with analytical overlay
- **2D animations**: Spatial evolution of particle distribution
- **Comparison animations**: Side-by-side method comparison

## Expected Results

A good SPH implementation should:
- Capture sharp shock front (x ≈ 0.85 at t=0.2)
- Maintain contact discontinuity (x ≈ 0.68 at t=0.2)
- Resolve rarefaction fan structure
- Preserve planar symmetry (no y-dependence)
- Match analytical solution within numerical error

**Key Metrics:**
- Shock position: x_shock ≈ 0.85 (exact: 0.8504)
- Contact position: x_contact ≈ 0.68 (exact: 0.6838)
- Shock compression: ρ_behind/ρ_ahead ≈ 3.03
- Post-shock velocity: u_behind ≈ 0.293

**Planar Symmetry Check:**
- Quantities should be independent of y
- Y-averaging should not change results
- Minimal scatter in y-direction

## Analytical Solution

The `shock_tube_2d_analytical.py` script computes the exact Riemann solution using an iterative solver.

**Method:**
- Exact Riemann solver (Toro 1999)
- Handles rarefaction, shock, and contact waves
- Arbitrary γ support (default: 1.4)
- Y-averaging for 2D comparison

**Output:**
Four-panel comparison plots showing:
1. Density ρ(x) - SPH vs exact
2. Velocity u(x) - SPH vs exact
3. Pressure P(x) - SPH vs exact
4. Internal energy e(x) - SPH vs exact

## Common Issues

**Excessive diffusion:**
- Shock/contact too wide
- → Increase resolution or reduce artificial viscosity

**Oscillations:**
- Unphysical wiggles near discontinuities
- → Adjust Balsara switch or limiter parameters

**Y-direction artifacts:**
- Waves breaking planar symmetry
- → Check particle distribution, boundary conditions

**Contact discontinuity:**
- Most challenging feature to resolve
- Expect some smearing even with good methods

## Comparison with 1D Shock Tube

**2D Advantages:**
- Tests planar symmetry preservation
- Validates 2D gradient calculations
- Checks for multi-dimensional artifacts

**2D Challenges:**
- Higher computational cost (8,000 vs 800 particles)
- Requires proper y-boundary conditions
- More complex visualization

## Directory Structure

```
shock_tube_2d/
├── config/
│   └── presets/
│       ├── gsph_wendland.json          # GSPH + Wendland config
│       ├── gsph_cubic.json             # GSPH + Cubic Spline config
│       ├── ssph_wendland.json
│       ├── disph_wendland.json
│       ├── gdisph_wendland.json
│       ├── gdisph_balsara_wendland.json
│       └── ...
├── scripts/
│   ├── shock_tube_2d_analytical.py    # Exact Riemann solver
│   ├── process_all_snapshots.py       # Batch 1D analytical overlays
│   ├── plot_2d_spatial.py             # 2D spatial visualization
│   ├── generate_all_animations.py     # Create all animations
│   ├── compare_shock_tube_2d.py       # Multi-method comparison
│   └── animate_comparison.py          # Comparison animations
├── results/                           # Simulation output (created on run)
├── shock_tube_2d.json                 # Active configuration (auto-generated)
├── Makefile.shock_tube_2d             # Build targets
└── README.md                          # This file
```

## References

1. **Sod, G. A. (1978)** "A Survey of Several Finite Difference Methods for Systems of Nonlinear Hyperbolic Conservation Laws", *Journal of Computational Physics*, 27(1), 1-31
2. **Toro, E. F. (1999)** "Riemann Solvers and Numerical Methods for Fluid Dynamics"
3. **Price, D. J. (2012)** "Smoothed Particle Hydrodynamics and Magnetohydrodynamics", *Journal of Computational Physics*, 231, 759-794
4. **Springel, V. (2010)** "Smoothed Particle Hydrodynamics in Astrophysics" (Section on shock tests)

## Notes

- The 2D shock tube should show perfect planar symmetry
- Y-averaging is used to compare with 1D analytical solution
- Small deviations from 1D behavior indicate numerical issues
- This test is particularly sensitive to gradient estimation in 2D
- Contact discontinuity resolution is a key discriminator between methods
