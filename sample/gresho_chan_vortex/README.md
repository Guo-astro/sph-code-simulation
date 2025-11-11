# Gresho-Chan Vortex Test

A 2D hydrodynamic test case featuring a rotating vortex with discontinuous angular velocity profile. This test evaluates the ability of SPH methods to maintain vortex stability and preserve angular momentum.

## Overview

The Gresho-Chan vortex consists of three regions with different angular velocities:
- Inner region (r < 0.2): Solid body rotation
- Middle region (0.2 ≤ r < 0.4): Differential rotation
- Outer region (r ≥ 0.4): Zero velocity

This test is particularly challenging for SPH methods due to the velocity discontinuities and the need to maintain pressure equilibrium.

## Configuration

**Domain**: [-0.5, 0.5] × [-0.5, 0.5] with periodic boundaries
**Resolution**: 128 × 128 particles (16,384 total)
**End time**: t = 3.0
**Neighbor number**: 50
**Gamma**: 5/3 (ideal gas)

## SPH Methods

The test includes 5 SPH formulations:

1. **GSPH** - Gather SPH (traditional)
2. **SSPH** - Scatter SPH  
3. **DISPH** - Density-Independent SPH
4. **GDISPH** - Gather-Density-Independent SPH
5. **GDISPH+Balsara** - GDISPH with Balsara viscosity switch

## Kernel Comparison

Each method is run with two different smoothing kernels:
- **Wendland C2** kernel
- **Cubic Spline** kernel

This results in 10 total simulations (5 methods × 2 kernels).

## Quick Start

### Run all simulations
```bash
make gresho_all
```

### Run specific kernel
```bash
make gresho_wendland    # All methods with Wendland kernel
make gresho_cubic       # All methods with Cubic Spline kernel
```

### Run specific method
```bash
make gresho_gsph_wendland
make gresho_ssph_cubic
# ... etc
```

### Generate visualizations
```bash
make gresho_kernel_viz        # Comparison plots
make gresho_kernel_animate    # Animation GIF
```

## Directory Structure

```
gresho_chan_vortex/
├── config/                          # Configuration files
│   ├── gsph_wendland.json
│   ├── gsph_cubic.json
│   ├── ssph_wendland.json
│   ├── ssph_cubic.json
│   ├── disph_wendland.json
│   ├── disph_cubic.json
│   ├── gdisph_wendland.json
│   ├── gdisph_cubic.json
│   ├── gdisph_balsara_wendland.json
│   └── gdisph_balsara_cubic.json
├── scripts/                         # Visualization scripts
│   ├── compare_kernels.py          # Generate comparison plots
│   └── animate_kernel_comparison.py # Generate animation
├── results/                         # Simulation outputs
│   ├── gsph_wendland/
│   ├── gsph/                       # Cubic Spline results
│   ├── ssph_wendland/
│   ├── ssph/
│   ├── disph_wendland/
│   ├── disph/
│   ├── gdisph_wendland/
│   ├── gdisph/
│   ├── gdisph_balsara_wendland/
│   ├── gdisph_balsara/
│   └── kernel_comparison/          # Comparison plots & animations
├── Makefile.gresho                  # Build targets
└── README.md                        # This file
```

## Expected Results

A well-performing SPH method should:
- Maintain the vortex structure over time
- Preserve angular momentum
- Keep pressure equilibrium across discontinuities
- Minimize artificial viscosity dissipation

The kernel comparison helps assess:
- Kernel influence on vortex stability
- Angular momentum conservation differences
- Numerical dissipation characteristics

## Visualization

The kernel comparison tools generate:
- **Comparison plots**: Side-by-side snapshots at key times
- **Animation**: Full evolution comparison (GIF format)
- **Final state**: High-resolution comparison at t=3.0

All visualizations use a logarithmic density colormap with consistent scaling for direct comparison.

## Physics

The vortex is initialized with:
- Density: ρ = 1.0 everywhere
- Pressure: Set to maintain hydrostatic equilibrium
- Angular velocity profile with discontinuities at r = 0.2 and r = 0.4

The test duration (t = 3.0) is approximately 1.5 rotation periods for the inner region.

## References

- Gresho, P. M., & Chan, S. T. (1990). "On the theory of semi-implicit projection methods for viscous incompressible flow"
- Hopkins, P. F. (2015). "A new class of accurate, mesh-free hydrodynamic simulation methods"
