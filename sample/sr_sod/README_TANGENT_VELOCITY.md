# Special Relativistic Tangent Velocity Shock Tube Test

This directory contains preset configurations for the **Riemann Problem 5: One-Dimensional Shock Tube Problem With Tangential Velocity** from Kitajima et al. (2025) arXiv:2510.18251v1 Section 3.1.5.

## Physical Background

In relativistic hydrodynamics, the velocity component along the partition (`v^t`) affects the solution, unlike non-relativistic cases. When the initial tangential velocity is large, errors in numerical solutions become larger and extremely high resolution is required (Pons et al. 2000, Mignone et al. 2005).

The Lorentz factor includes tangential velocity:
```
γ = 1 / √(1 - (v_x² + v_t²) / c²)
```

## Initial Conditions

This test uses the same pressure conditions as the Strong Relativistic Blast Wave (Section 3.1.3), but with tangential velocity:

| Side  | Pressure (P) | Density (n) | Normal Velocity (v^x) | Tangent Velocity (v^t) |
|-------|--------------|-------------|----------------------|----------------------|
| Left  | 1000         | 1.0         | 0                    | {0, 0.9, 0.99}      |
| Right | 0.01         | 1.0         | 0                    | {0, 0.9, 0.99}      |

With γ = 5/3 (ratio of specific heats).

## Available Preset Files

### Standard Resolution (1600 particles per side)

| File | v^t_L | v^t_R | γ_L | γ_R |
|------|-------|-------|-----|-----|
| `sr_tangent_vt0_vt0.json` | 0 | 0 | 1.00 | 1.00 |
| `sr_tangent_vt09_vt0.json` | 0.9 | 0 | 2.29 | 1.00 |
| `sr_tangent_vt09_vt09.json` | 0.9 | 0.9 | 2.29 | 2.29 |
| `sr_tangent_vt09_vt099.json` | 0.9 | 0.99 | 2.29 | 7.09 |
| `sr_tangent_vt099_vt099.json` | 0.99 | 0.99 | 7.09 | 7.09 |

### High Resolution (51200 particles per side)

| File | v^t_L | v^t_R | Notes |
|------|-------|-------|-------|
| `sr_tangent_vt09_vt09_hires.json` | 0.9 | 0.9 | Resolves rarefaction wave (Fig. 9) |

## Running the Test

```bash
# Build with DIM=2 (required for tangent velocity)
cd /path/to/sphcode/build
cmake -DSPH_DIM=2 ..
make -j4

# Run a specific test
./sph sample/sr_sod/config/presets/sr_tangent_vt09_vt09.json

# Results will be in sample/sr_sod/results/tangent_vt09_vt09/
```

## Expected Behavior

From Kitajima et al. (2025):
- Deviation from analytical solution increases when initial tangential velocity is large on the rarefaction wave side
- Increased resolution on the expansion wave side effectively resolves this issue
- The tangential velocity v^t near the contact discontinuity behaves worse in SPH-based methods compared to grid-based methods at low resolution (Liptai & Price 2019)

## Notes

- The simulation ends when the shock front reaches x = 0.2 (not a fixed end time)
- The shock velocity depends on tangent velocity
- Requires DIM >= 2 compilation (particles arranged in 1D along x-axis but have y-velocity component)

## References

- Kitajima et al. (2025): arXiv:2510.18251v1
- Pons et al. (2000): JFM 422, 125
- Rezzolla & Zanotti (2002): JFM 449, 395
- Mignone et al. (2005): ApJS 160, 199
- Zhang & MacFadyen (2006): ApJS 164, 255
- Liptai & Price (2019): MNRAS 485, 819
