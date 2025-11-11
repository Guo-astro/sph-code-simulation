# KHI (Kelvin-Helmholtz Instability) - Multi-Method & Kernel Comparison

## Overview
Comprehensive benchmark infrastructure for testing Kelvin-Helmholtz instability evolution across 5 SPH methods and 2 kernel types.

## Test Configuration
- **Dimensions**: 2D (256×256 = 65,536 particles)
- **Domain**: [0, 1] × [0, 1], periodic boundaries
- **Duration**: t = 0 → 3.0
- **Output**: Every 0.1 time units (31 snapshots)
- **Gamma**: 5/3 (ideal gas)
- **Neighbor Number**: 32
- **Special Features**: Time-dependent artificial viscosity enabled

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
sample/khi/
├── config/presets/          # 10 preset configurations (5 methods × 2 kernels)
├── scripts/                 # Visualization scripts
├── results/                 # Simulation output directories
└── Makefile.khi             # Build targets
```

## Quick Start

### Single Method Run
```bash
# Run with default (SSPH)
make khi_run

# Run specific method
make khi_run KHI_PRESET=khi_gsph
```

### Multi-Method Comparison (Cubic Spline)
```bash
# Complete workflow: run 5 methods + generate plots + animation
make khi_compare_all

# Individual steps
make khi_compare_run      # Run all 5 methods
make khi_compare_viz      # Generate comparison plots
make khi_compare_animate  # Generate animation
```

### Kernel Comparison (10 runs total)
```bash
# Complete workflow: run 10 simulations + visualizations
make khi_kernel_all

# Individual steps
make khi_kernel_run       # Run 5 methods × 2 kernels
make khi_kernel_viz       # Generate side-by-side plots
make khi_kernel_animate   # Generate comparison animation
```

## Available Presets

### Cubic Spline Kernel
- `khi_gsph.json`
- `khi_ssph.json`
- `khi_disph.json`
- `khi_gdisph.json`
- `khi_gdisph_balsara.json`

### Wendland Kernel
- `khi_gsph_wendland.json`
- `khi_ssph_wendland.json`
- `khi_disph_wendland.json`
- `khi_gdisph_wendland.json`
- `khi_gdisph_balsara_wendland.json`

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
- Captures vortex roll-up evolution

### Kernel Comparison
- Side-by-side: Wendland (left) vs Cubic Spline (right)
- 5 rows (methods) × 2 columns (kernels)
- Synchronized time evolution
- Direct visual comparison of instability growth

## Scientific Goals
- Test shear flow instability evolution
- Compare vortex resolution across methods
- Evaluate mixing and turbulence development
- Measure kinetic energy dissipation
- Assess numerical diffusion effects
- Compare kernel impact on small-scale structure

## KHI Physics
The Kelvin-Helmholtz instability develops when two fluids with different velocities interact:
- Initial setup: Two layers with velocity shear
- Perturbation grows exponentially (linear phase)
- Vortices roll up (nonlinear phase)
- Eventually turbulent cascade develops
- Tests: gradient accuracy, dissipation control, small-scale capturing

## Cleanup
```bash
make khi_clean              # Clean all results
make khi_compare_clean      # Clean comparison results only
make khi_kernel_clean       # Clean kernel comparison only
```

## Requirements
- DIM=2 build configuration
- Python 3 with numpy, matplotlib
- Optional: tqdm (progress bars), pillow (animations)
- **Note**: Higher resolution (256×256) requires more memory than other tests

## Performance
- Typical runtime: ~30-90 seconds per simulation (65,536 particles, 31 snapshots)
- Total workflow time: ~5-15 minutes for 10 runs
- Plot generation: ~30-60 seconds
- Animation generation: ~2-5 minutes
- **Memory**: ~500MB-1GB per simulation

## Expected Behavior

### Method Characteristics
- **SSPH**: Smooth but diffusive, vortices may be damped
- **DISPH**: Better kinetic energy conservation, sharper features
- **GSPH**: Excellent shock/contact capturing, minimal diffusion
- **GDISPH**: Combines DISPH stability with Riemann accuracy
- **GDISPH+Balsara**: Reduced AV in vortices, cleaner small-scale structure

### Kernel Effects
- **Cubic Spline**: Standard resolution, moderate smoothness
- **Wendland**: Smoother gradients, may affect instability growth rate
- Kernel choice affects: vortex core size, mixing layer thickness, energy dissipation

## Diagnostics to Check
1. Vortex roll-up timing (linear growth rate)
2. Number and size of vortices formed
3. Symmetry preservation
4. Kinetic energy history
5. Enstrophy evolution
6. Mixing layer thickness growth

## Notes
- KHI is sensitive to numerical dissipation
- Time-dependent AV helps reduce excess viscosity
- Higher resolution captures more vortex scales
- This test separates low-diffusion from high-diffusion methods
- Kernel smoothness impacts sub-grid turbulence modeling
