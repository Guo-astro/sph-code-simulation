# Hydrostatic Test - Multi-Method & Kernel Comparison

## Overview
Comprehensive benchmark infrastructure for testing hydrostatic equilibrium across 5 SPH methods and 2 kernel types.

## Test Configuration
- **Dimensions**: 2D (64×64 = 4096 particles)
- **Domain**: [-0.5, 0.5] × [-0.5, 0.5], periodic boundaries
- **Duration**: t = 0 → 8.0
- **Output**: Every 0.2 time units (41 snapshots)
- **Gamma**: 5/3 (ideal gas)
- **Neighbor Number**: 32

## SPH Methods (5)
1. **GSPH** - Godunov SPH (2nd order Riemann solver)
2. **SSPH** - Standard SPH (traditional formulation)
3. **DISPH** - Density-Independent SPH (pressure-energy formulation)
4. **GDISPH** - Godunov DISPH (1st order Riemann + DISPH)
5. **GDISPH+Balsara** - GDISPH with Balsara switch (density-based AV limiter)

## Kernel Types (2)
- **Cubic Spline** (M₄) - Classic SPH kernel, C¹ continuous
- **Wendland** (C²) - Modern smooth kernel, C² continuous, compact support

## Directory Structure
```
sample/hydrostatic/
├── config/presets/          # 10 preset configurations (5 methods × 2 kernels)
├── scripts/                 # Visualization scripts
├── results/                 # Simulation output directories
└── Makefile.hydrostatic     # Build targets
```

## Quick Start

### Single Method Run
```bash
# Run with default (SSPH)
make hydrostatic_run

# Run specific method
make hydrostatic_run HYDROSTATIC_PRESET=hydrostatic_gsph
```

### Multi-Method Comparison (Cubic Spline)
```bash
# Complete workflow: run 5 methods + generate plots + animation
make hydrostatic_compare_all

# Individual steps
make hydrostatic_compare_run      # Run all 5 methods
make hydrostatic_compare_viz      # Generate comparison plots
make hydrostatic_compare_animate  # Generate animation
```

### Kernel Comparison (10 runs total)
```bash
# Complete workflow: run 10 simulations + visualizations
make hydrostatic_kernel_all

# Individual steps
make hydrostatic_kernel_run       # Run 5 methods × 2 kernels
make hydrostatic_kernel_viz       # Generate side-by-side plots
make hydrostatic_kernel_animate   # Generate comparison animation
```

## Available Presets

### Cubic Spline Kernel
- `hydrostatic_gsph.json`
- `hydrostatic_ssph.json`
- `hydrostatic_disph.json`
- `hydrostatic_gdisph.json`
- `hydrostatic_gdisph_balsara.json`

### Wendland Kernel
- `hydrostatic_gsph_wendland.json`
- `hydrostatic_ssph_wendland.json`
- `hydrostatic_disph_wendland.json`
- `hydrostatic_gdisph_wendland.json`
- `hydrostatic_gdisph_balsara_wendland.json`

## Output Structure
```
results/
├── gsph/                    # GSPH with cubic spline
├── ssph/                    # SSPH with cubic spline
├── disph/                   # DISPH with cubic spline
├── gdisph/                  # GDISPH with cubic spline
├── gdisph_balsara/          # GDISPH+Balsara with cubic spline
├── gsph_wendland/           # GSPH with Wendland
├── ssph_wendland/           # SSPH with Wendland
├── disph_wendland/          # DISPH with Wendland
├── gdisph_wendland/         # GDISPH with Wendland
├── gdisph_balsara_wendland/ # GDISPH+Balsara with Wendland
├── comparison/              # 5-method comparison plots & animation
└── kernel_comparison/       # Kernel comparison plots & animation
```

## Visualization

### Method Comparison Plots
- Hot colormap with log₁₀(density)
- Black background for enhanced contrast
- 5-panel layout showing all methods
- Multiple time snapshots + final state

### Kernel Comparison
- Side-by-side: Wendland (left) vs Cubic Spline (right)
- 5 rows (methods) × 2 columns (kernels)
- Synchronized time evolution
- Direct visual comparison of kernel effects

## Scientific Goals
- Test hydrostatic equilibrium maintenance
- Compare dissipation characteristics across methods
- Evaluate kernel smoothness effects
- Measure pressure field stability
- Assess energy conservation

## Cleanup
```bash
make hydrostatic_clean              # Clean all results
make hydrostatic_compare_clean      # Clean comparison results only
make hydrostatic_kernel_clean       # Clean kernel comparison only
```

## Requirements
- DIM=2 build configuration
- Python 3 with numpy, matplotlib
- Optional: tqdm (progress bars), pillow (animations)

## Performance
- Typical runtime: ~1-3 seconds per simulation (4096 particles, 41 snapshots)
- Total workflow time: ~20-30 seconds for 10 runs
- Plot generation: ~10-30 seconds
- Animation generation: ~1-2 minutes

## Notes
- Hydrostatic test validates pressure gradient accuracy
- Methods with better gradient operators maintain equilibrium longer
- Kernel smoothness affects long-term stability
- Balsara switch reduces artificial viscosity in shear-free regions
