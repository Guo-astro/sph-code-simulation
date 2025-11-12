# 1D Vacuum Test

One-dimensional vacuum formation test from the DISPH paper (Yuasa & Mori 2024, Section 4.2.2). Tests SPH methods' ability to handle vacuum regions and rarefaction waves.

## Important: Requires 1D Build

⚠️ **This simulation requires `DIM=1` compilation.** The vacuum test is a 1D problem.

Before running:

1. **Reconfigure build for 1D**:
   ```bash
   sed -i '' 's/#define DIM [0-9]/#define DIM 1/' include/defines.hpp
   cd build && make -j8 && cd ..
   ```

2. **Run vacuum test**:
   ```bash
   make vacuum_run
   ```

3. **Visualize results**:
   ```bash
   make vacuum_plot_all
   ```

4. **Multi-method comparison**:
   ```bash
   make vacuum_compare_all
   ```

5. **Kernel comparison**:
   ```bash
   make vacuum_kernel_all
   ```

## Physics

**Vacuum Formation Problem:**
- Left state (x ≤ 0): ρ = 1.0, P = 0.4, v = -2.0
- Right state (x > 0): ρ = 1.0, P = 0.4, v = 2.0
- γ = 1.4 (diatomic gas)

**Expected behavior:**
- Left fluid moves left at v = -2.0
- Right fluid moves right at v = +2.0
- Vacuum region forms at x = 0
- Two symmetric rarefaction waves propagate outward

**Challenge:**
This test is particularly challenging because:
- True vacuum (ρ → 0) develops at the center
- SPH struggles with very low particle density regions
- Godunov methods can have accuracy issues in vacuum regime
- Tests the robustness of pressure gradient calculations

## Configuration

**Domain:** [-0.5, 1.5] with periodic boundaries
**Particles:** 400 (left region) + 400 (right region) = 800 total
**Neighbor number:** Nngb = 5.2
**End time:** t = 0.14154
**Output interval:** Δt = 0.01

Presets are located in `config/presets/`:
- `vacuum_gsph_wendland.json` - GSPH with Wendland C4 kernel
- `vacuum_gsph_cubic.json` - GSPH with Cubic Spline kernel
- `vacuum_ssph_wendland.json` - SSPH with Wendland C4 kernel
- `vacuum_ssph_cubic.json` - SSPH with Cubic Spline kernel
- `vacuum_disph_wendland.json` - DISPH with Wendland C4 kernel
- `vacuum_disph_cubic.json` - DISPH with Cubic Spline kernel
- `vacuum_gdisph_wendland.json` - GDISPH with Wendland C4 kernel
- `vacuum_gdisph_cubic.json` - GDISPH with Cubic Spline kernel
- `vacuum_gdisph_balsara_wendland.json` - GDISPH+Balsara with Wendland C4 kernel
- `vacuum_gdisph_balsara_cubic.json` - GDISPH+Balsara with Cubic Spline kernel

## Visualization

The visualization scripts generate:
- **Density profile** - Should show vacuum region (ρ → 0) at center
- **Velocity profile** - Should show symmetric pattern (left: v < 0, right: v > 0)
- **Pressure profile** - Should decrease toward vacuum region
- **Internal energy profile** - Should track rarefaction wave structure

Analytical solution is computed using exact Riemann solver for the vacuum case.

## Quick Start

```bash
# Initial run (GSPH + Wendland)
make vacuum_run

# Generate all analytical comparison plots
make vacuum_plot_all

# Run all 5 SPH methods (Wendland kernel)
make vacuum_compare_run

# Generate comparison plots
make vacuum_compare_viz

# Create animations
make vacuum_compare_animate

# Complete workflow
make vacuum_compare_all

# Run all 10 variants (5 methods × 2 kernels)
make vacuum_kernel_all
```

## Output Structure

```
sample/vacuum/results/
├── gsph_wendland/
│   ├── snapshot_0000.csv
│   ├── snapshot_0001.csv
│   ├── ...
│   └── plots/
│       ├── snapshot_0000_comparison.png
│       └── ...
├── ssph_wendland/
├── disph_wendland/
├── gdisph_wendland/
├── gdisph_balsara_wendland/
├── gsph_cubic/
├── ssph_cubic/
├── disph_cubic/
├── gdisph_cubic/
├── gdisph_balsara_cubic/
├── comparison/
│   └── vacuum_multi_method_comparison.png
└── animations/
    ├── gsph_wendland_vacuum.gif
    └── ...
```

## Expected Results

According to the DISPH paper:
- **GSPH Case 3** shows moderate accuracy but may have issues in vacuum regime
- **SSPH** requires artificial viscosity tuning
- **DISPH** shows good performance but still needs artificial viscosity
- **GDISPH Case 1 & 2** should handle vacuum regions better than GSPH
- All Godunov methods may show ~50% overestimation error in vacuum regime (known limitation)

## Makefile Targets

```bash
make vacuum_help                    # Show all available targets
make vacuum_run                     # Run default (GSPH + Wendland)
make vacuum_plot_all               # Generate all analytical plots
make vacuum_clean                  # Clean results

# Multi-method comparison (5 methods, Wendland kernel)
make vacuum_compare_run            # Run all 5 methods
make vacuum_compare_viz            # Generate comparison plots
make vacuum_compare_animate        # Create animations
make vacuum_compare_all            # Complete comparison workflow
make vacuum_compare_clean          # Clean comparison results

# Kernel comparison (5 methods × 2 kernels = 10 runs)
make vacuum_kernel_run             # Run all 10 simulations
make vacuum_kernel_viz             # Generate kernel comparison plots
make vacuum_kernel_animate         # Create animations
make vacuum_kernel_all             # Complete kernel comparison
make vacuum_kernel_clean           # Clean kernel results
```

## References

- Yuasa, T., & Mori, M. (2024). "Novel Hydrodynamic Schemes Capturing Shocks and Contact Discontinuities and Comparison Study with Existing Methods." arXiv:2312.03224v3
- Toro, E. F. (2009). "Riemann Solvers and Numerical Methods for Fluid Dynamics"
- Saitoh, T. R., & Makino, J. (2013). "A density-independent formulation of smoothed particle hydrodynamics"

## Common Issues

**Issue:** Particles have NaN values or simulation crashes
- **Cause:** Vacuum region causes division by zero
- **Solution:** Check that SPH code handles ρ → 0 gracefully with minimum density floor

**Issue:** Poor resolution in vacuum region
- **Cause:** Too few particles in low-density areas
- **Solution:** This is inherent to SPH; expected behavior documented in paper

**Issue:** Large errors compared to analytical solution
- **Cause:** Godunov methods known to have issues in vacuum regime
- **Solution:** Expected behavior; GDISPH Case 1 should have best performance (~50% error vs ~75% for GSPH)
