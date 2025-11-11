# Sedov Blast Wave Test

A 2D benchmark test for strong point explosions, featuring a self-similar analytical solution. The Sedov blast wave tests the ability of SPH methods to handle strong shocks and extreme gradients.

## Overview

The Sedov-Taylor blast wave describes the expansion of a strong shock from a point explosion. The solution is self-similar, meaning the spatial structure scales with time according to a power law. This makes it an excellent benchmark for validating shock-capturing schemes.

**Key Features:**
- Strong shock propagation
- Self-similar analytical solution
- Tests artificial viscosity and shock capturing
- Extreme density and pressure gradients

## Physics

**Initial Conditions:**
- Uniform density: ρ₀ = 1.0
- Low background pressure: p₀ = 10⁻⁶
- Point energy injection: E = 1.0 in central region

**Governing Equations:**
The solution follows from dimensional analysis:
```
R_shock(t) ∝ (E₀ t² / ρ₀)^(1/(ν+2))
```
where ν = 2 for cylindrical geometry (2D simulation).

**Analytical Solution:**
The Sedov solution provides exact profiles for:
- Density ρ(r,t)
- Velocity v(r,t)
- Pressure p(r,t)
- Internal energy e(r,t)

## Configuration

**Domain**: [-0.5, 0.5] × [-0.5, 0.5], non-periodic
**Resolution**: 100 × 100 particles (10,000 total)
**End time**: t = 0.1
**Output interval**: Δt = 0.005 (21 snapshots)
**Gamma**: γ = 1.4 (ideal gas)
**Kernel**: Wendland C2 or Cubic Spline

## Quick Start

### 🚀 ONE-SHOT: Complete Multi-Method Comparison
```bash
make sedov_compare_all          # Run all 4 SPH methods + analytical overlays + animation
                                # Includes: GSPH, SSPH, DISPH, GDISPH (Wendland kernel)
                                # Generates: 84 analytical comparison plots (21 snapshots × 4 methods)
```

### 🚀 ONE-SHOT: Complete Kernel Comparison
```bash
make sedov_kernel_all           # Run all 8 simulations (4 methods × 2 kernels) + visualizations
                                # Includes: All methods with Wendland and Cubic Spline kernels
                                # Generates: 168 analytical comparison plots (21 snapshots × 8 configs)
```

### Single Method Runs
```bash
make sedov_run                  # Default run (GSPH + Wendland)
make sedov_gsph_wendland        # GSPH with Wendland kernel
make sedov_gsph_cubic           # GSPH with Cubic Spline kernel
make sedov_ssph_wendland        # SSPH with Wendland kernel
make sedov_disph_wendland       # DISPH with Wendland kernel
make sedov_gdisph_wendland      # GDISPH with Wendland kernel
```

### Multi-Method Comparison (Wendland kernel)
```bash
make sedov_compare_run          # Run all 4 methods: GSPH, SSPH, DISPH, GDISPH
make sedov_compare_viz          # Generate analytical overlays for all methods
make sedov_compare_animate      # Create comparison animation
```

### Kernel Comparison (Wendland vs Cubic Spline)
```bash
make sedov_kernel_run           # Run all 8 simulations (4 methods × 2 kernels)
make sedov_kernel_viz           # Generate analytical overlays for all configs
make sedov_kernel_animate       # Create kernel comparison animation
```

### Analytical Comparison Plots
```bash
make sedov_plot                 # Single snapshot comparison
make sedov_plot_all             # All snapshots with analytical overlay
```

### Cleanup
```bash
make sedov_clean                # Remove all results
make sedov_compare_clean        # Remove comparison results only
make sedov_kernel_clean         # Remove kernel comparison results only
```

## What the ONE-SHOT Commands Do

### `make sedov_compare_all`
1. **Runs 4 simulations** with Wendland kernel:
   - GSPH, SSPH, DISPH, GDISPH
   - Each generates 21 snapshots (t = 0 to 0.1, Δt = 0.005)
2. **Generates analytical overlays**:
   - 4-panel comparison plots for each snapshot
   - Density, velocity, pressure (log), internal energy (log)
   - SPH particles (blue scatter) vs Sedov solution (red curves)
   - Total: 84 comparison plots (21 × 4)
3. **Creates comparison animation** (placeholder)
4. **Total runtime**: ~10-15 seconds (4 simulations + plotting)

### `make sedov_kernel_all`
1. **Runs 8 simulations** (4 methods × 2 kernels):
   - GSPH, SSPH, DISPH, GDISPH
   - Wendland C2 and Cubic Spline kernels
   - Each generates 21 snapshots
2. **Generates analytical overlays**:
   - Same 4-panel plots for all 8 configurations
   - Total: 168 comparison plots (21 × 8)
3. **Creates kernel comparison** (placeholder)
4. **Total runtime**: ~20-30 seconds (8 simulations + plotting)

### Output Structure
```
sample/sedov/results/
├── gsph_wendland/
│   ├── snapshot_*.csv         # Simulation data
│   ├── energy.dat             # Energy history
│   └── plots/
│       └── snapshot_*_comparison.png  # 21 analytical overlays
├── ssph_wendland/
│   └── plots/                 # 21 analytical overlays
├── disph_wendland/
│   └── plots/                 # 21 analytical overlays
├── gdisph_wendland/
│   └── plots/                 # 21 analytical overlays
├── gsph_cubic/                # (kernel comparison only)
├── ssph_cubic/                # (kernel comparison only)
├── disph_cubic/               # (kernel comparison only)
├── gdisph_cubic/              # (kernel comparison only)
└── comparison/
    └── comparison_animation.gif (placeholder)
```

### Manual plotting
```bash
# Single snapshot
python3 sample/sedov/scripts/sedov_analytical.py \
    sample/sedov/results/gsph_wendland/snapshot_0020.csv \
    -o comparison.png

# All snapshots
python3 sample/sedov/scripts/process_all_snapshots.py \
    sample/sedov/results/gsph_wendland
```

## Directory Structure

```
sedov/
├── config/
│   └── presets/
│       ├── gsph_wendland.json     # GSPH + Wendland config
│       └── gsph_cubic.json        # GSPH + Cubic Spline config
├── scripts/
│   ├── sedov_analytical.py        # Analytical solution calculator
│   └── process_all_snapshots.py   # Batch processing script
├── results/
│   ├── gsph_wendland/             # Simulation output
│   │   ├── snapshot_*.csv
│   │   └── plots/                 # Analytical comparison plots
│   └── gsph/                      # Cubic Spline output
├── sedov.json                      # Default configuration
├── Makefile.sedov                  # Build targets
└── README.md                       # This file
```

## Analytical Solution

The `sedov_analytical.py` script computes the exact Sedov-Taylor solution and overlays it on SPH results.

**Theory:**
- Based on Sedov (1959) similarity solution
- Implements Kamm & Timmes (2007) formulation
- Accounts for cylindrical geometry (ν = 2)
- Handles arbitrary γ (default: 1.4)

**Output:**
Four-panel comparison plots showing:
1. **Density profile**: ρ(r)
2. **Velocity profile**: v(r)
3. **Pressure profile**: p(r) (log scale)
4. **Internal energy**: e(r) (log scale)

Each plot shows:
- Blue dots: SPH simulation particles
- Red line: Analytical solution

## Expected Results

A good SPH implementation should:
- Capture the shock front position accurately
- Maintain proper shock compression ratios
- Preserve spherical symmetry
- Show smooth profiles behind the shock
- Match analytical solution within numerical error

**Key Metrics:**
- Shock radius: R_shock(t) ∝ t^(2/4) = t^0.5
- Density jump at shock: ρ/ρ₀ = (γ+1)/(γ-1) ≈ 6 for γ=1.4
- Velocity at shock: v_shock = (2/(γ+1)) × dR/dt

## Common Issues

**Excessive diffusion:**
- Shock front too wide
- Reduced compression ratios
- → Increase resolution or reduce artificial viscosity

**Oscillations:**
- Unphysical wiggles behind shock
- → Adjust Balsara switch or kernel choice

**Asymmetry:**
- Non-circular shock front
- → Check particle distribution, increase neighbor number

## References

1. **Sedov, L. I. (1959)** "Similarity and Dimensional Methods in Mechanics"
2. **Taylor, G. I. (1950)** "The Formation of a Blast Wave by a Very Intense Explosion"
3. **Kamm, J. R. & Timmes, F. X. (2007)** "On Efficient Generation of Numerically Robust Sedov Solutions"
4. **Springel, V. (2010)** "Smoothed Particle Hydrodynamics in Astrophysics" (Section on shock tests)

## Notes

- The Sedov test is scale-free: results should match for any consistent choice of E₀, ρ₀
- Early times (t < 0.01) may deviate from analytical solution due to finite blast radius
- Late times show excellent agreement if shock capturing is properly implemented
- This test is particularly sensitive to artificial viscosity parameters
