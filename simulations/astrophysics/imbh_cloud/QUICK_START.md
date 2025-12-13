# IMBH-Cloud Quick Start Guide

## TL;DR

Simulate IMBH (10⁵ M_☉) tidal disruption of molecular cloud (10⁴ M_☉, 5 pc):

```bash
# 1. Build (set DIM=3 in include/defines.hpp)
make clean && make

# 2. Run
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_b3pc_run

# 3. Visualize
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_visualize
```

## Physics Setup

- **Cloud**: Lane-Emden n=1.5, M=10⁴ M_☉, R=5 pc
- **IMBH**: M=10⁵ M_☉ at distance b=3-6 pc
- **Thermal**: Koyama & Inutsuka equilibrium (τ_cool << τ_tidal)
- **SPH**: GDISPH method (N = 2×10⁵ particles)

## Resolution

| Purpose | N particles | Runtime | h/R_cloud |
|---------|-------------|---------|-----------|
| Quick test | 5×10⁴ | ~2 hr | 0.014 |
| Production | 2×10⁵ | ~12 hr | 0.007 |
| High-res | 10⁶ | ~3 days | 0.003 |

**Recommended**: N = 2×10⁵ (resolves ~100 particles per Jeans mass)

## Files Structure

```
simulations/astrophysics/imbh_cloud/
├── README.md                        ← Read this first
├── docs/
│   ├── RESEARCH_SETUP.md           ← Physics & resolution analysis
│   └── IMPLEMENTATION_SUMMARY.md   ← Code architecture
├── config/presets/
│   ├── imbh_cloud_b3pc_gdisph.json ← b = 3 pc preset
│   └── imbh_cloud_b6pc_gdisph.json ← b = 6 pc preset
├── scripts/
│   └── visualize_disruption.py     ← Plot generator
├── Makefile.imbh_cloud             ← Build automation
└── results/                         ← Output (created at runtime)
```

## Makefile Targets

```bash
# See all options
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_help

# Single runs
imbh_b3pc_run          # b = 3 pc
imbh_b6pc_run          # b = 6 pc

# Parameter scan
imbh_scan_run          # b = 3, 4, 5, 6 pc (4 runs)

# Analysis
imbh_visualize         # Generate plots
imbh_clean             # Remove results
```

## Expected Output

### Snapshots
- `results/b3pc_gdisph/snapshot_NNNN.csv` (every 0.02 Myr)
- `results/b3pc_gdisph/energy.dat` (energy evolution)

### Plots
- `results/plots/density_evolution.png` (max/mean density vs time)
- `results/plots/axis_ratios.png` (cloud deformation a:b:c)
- `results/plots/phase_diagram.png` (n-T distribution)

### Animations
- `results/animations/disruption_animation.gif` (XY and XZ projections)

## Integration Needed

Before first compile, add to `src/CMakeLists.txt`:
```cmake
add_subdirectory(external_forces)
```

And register sample in `src/solver.cpp`:
```cpp
} else if (sample_type == "imbh_cloud") {
    make_imbh_cloud();
}
```

## Key Physics Diagnostics

1. **Tidal radius**: r_t ~ R_cloud (M_BH/M_cloud)^(1/3) ~ 10 pc
2. **Hill radius**: r_H ~ b (M_cloud/3M_BH)^(1/3) ~ 0.7 pc
3. **Thermal time**: τ_cool ~ 10³-10⁴ yr << τ_tidal ~ 10⁵ yr ✓
4. **Energy conservation**: ΔE/E < 1% (check energy.dat)

## Parameter Customization

Edit JSON preset to change:
```json
{
  "sample": {
    "N": 200000,          // Particle count
    "M_cloud": 10000,     // Cloud mass [M_☉]
    "R_cloud": 5.0,       // Cloud radius [pc]
    "M_BH": 100000,       // BH mass [M_☉]
    "b": 3.0,             // Impact parameter [pc]
    "v_rel": 0.0,         // Relative velocity
    "epsilon_BH": 0.05    // BH softening [pc]
  }
}
```

## Validation Checklist

- [ ] Energy conserved (< 1% drift)
- [ ] Thermal equilibrium maintained (n-T follows K&I curve)
- [ ] Tidal elongation along x-axis (a > b > c)
- [ ] Pancake compression in yz-plane (c << a)

## References

- **Physics**: `docs/RESEARCH_SETUP.md`
- **Usage**: `README.md`
- **Code**: `docs/IMPLEMENTATION_SUMMARY.md`

---

**Need help?** Check README.md for detailed explanations.
