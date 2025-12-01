# IMBH-Molecular Cloud Tidal Disruption Simulation

Simulates tidal disruption of molecular clouds by intermediate-mass black holes (IMBHs).

## Overview

This test case studies the tidal interaction between:
- **IMBH**: M_BH = 10⁵ M_☉ (intermediate-mass black hole)
- **Molecular Cloud**: M_cloud = 10⁴ M_☉, R_cloud = 5 pc
- **Impact Parameter**: b = 3-6 pc (pericenter distance)

### Physics Included

1. **Self-gravity**: Cloud's internal gravity (Lane-Emden n=1.5 polytrope)
2. **External force**: Point-mass BH gravity with Plummer softening
3. **Thermal equilibrium**: Koyama & Inutsuka (2000) cooling/heating
4. **SPH method**: GDISPH (best for tidal disruption dynamics)

## Key Features

- **3D Lane-Emden sphere** with polytropic index n = 5/3 (γ = 5/3)
- **Thermal timescale << dynamical timescale** (cloud stays in thermal equilibrium)
- **Resolution**: N = 2×10⁵ particles (resolves ~100 particles per Jeans mass)
- **Parameter scans**: Impact parameters b = 3, 4, 5, 6 pc

## Quick Start

### 1. Build Code (3D Configuration)

```bash
# Edit include/defines.hpp and set:
#define DIM 3

# Rebuild
make clean && make
```

### 2. Run Single Simulation

```bash
# Impact parameter b = 3 pc
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_b3pc_run

# Impact parameter b = 6 pc
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_b6pc_run
```

### 3. Run Parameter Scan

```bash
# Run all impact parameters: b = 3, 4, 5, 6 pc
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_scan_run
```

### 4. Generate Visualizations

```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_visualize
```

## Expected Physics

### Tidal Disruption Sequence

1. **Initial Approach** (t = 0 - 0.5 t_tidal):
   - Cloud experiences increasing tidal gradient
   - Elongation along BH direction (x-axis)
   - Compression perpendicular to orbit (yz-plane)

2. **Maximum Compression** (t ~ 0.5 - 1.0 Myr):
   - Peak density at pericenter passage
   - Pancake-like compression
   - Potential shock heating (if relative velocity is high)

3. **Disruption** (t > 1.0 Myr):
   - Leading and trailing tidal tails form
   - Separation into bound vs. unbound material
   - Remnant cloud mass determined by self-gravity

4. **Thermal Response**:
   - Compression increases density → T_eq decreases (CNM branch)
   - Cooling maintains thermal equilibrium
   - No runaway heating (unlike adiabatic case)

### Key Timescales

- **Tidal timescale**: t_tidal ~ √(R³/(GM_BH)) ~ 10⁴-10⁵ yr
- **Crossing time**: t_cross ~ R/v_rel ~ 10⁵-10⁶ yr  
- **Thermal time**: t_thermal ~ 10³-10⁴ yr << t_tidal ✓

**Assumption**: Thermal equilibrium maintained throughout.

## Resolution Analysis

### Jeans Mass Criterion

For CNM (n_H ~ 100-1000 cm⁻³, T ~ 10-50 K):
- Jeans mass: M_J ~ 1-10 M_☉
- **Requirement**: N_J ≥ 50-100 particles per Jeans mass
- **For M_cloud = 10⁴ M_☉**: N ≥ 5×10⁴ (minimum), N = 2×10⁵ (recommended)

### Spatial Resolution

- Smoothing length: h ~ R_cloud/√N ~ 0.01 R_cloud (for N = 2×10⁵)
- **Tidal radius**: r_t ~ R_cloud (M_BH/M_cloud)^(1/3) ~ 5 R_cloud
- **Resolution**: h/r_t ~ 0.002 ✓ (well-resolved)

### Performance

| N particles | Runtime (approx) | Memory | Recommended Use |
|-------------|------------------|--------|-----------------|
| 5×10⁴      | ~2 hours         | ~1 GB  | Quick tests, parameter scans |
| 2×10⁵      | ~12 hours        | ~4 GB  | Production runs |
| 10⁶        | ~3 days          | ~20 GB | High-resolution studies |

*Approximate for 1 Myr simulation on 8-core workstation*

## Configuration Files

### Presets (`config/presets/`)

- `imbh_cloud_b3pc_gdisph.json` - Impact parameter b = 3 pc
- `imbh_cloud_b6pc_gdisph.json` - Impact parameter b = 6 pc

### Key Parameters

```json
{
  "sample": {
    "type": "imbh_cloud",
    "N": 200000,              // Particle count
    "M_cloud": 10000.0,       // Cloud mass [M_☉]
    "R_cloud": 5.0,           // Cloud radius [pc]
    "M_BH": 100000.0,         // BH mass [M_☉]
    "b": 3.0,                 // Impact parameter [pc]
    "v_rel": 0.0,             // Relative velocity [code units]
    "epsilon_BH": 0.05        // BH softening [pc]
  },
  "numerical": {
    "sph_type": "gdisph",     // GDISPH method
    "kernel": "wendland",     // Wendland C4 kernel
    "neighbor_number": 50
  },
  "thermal": {
    "enable_cooling": true,
    "cooling_type": "koyama_inutsuka",
    "N_H_column": 1e20,
    "thermal_relax_time": 0.01
  }
}
```

## Output Structure

```
sample/imbh_cloud/results/
├── b3pc_gdisph/
│   ├── snapshot_0000.csv
│   ├── snapshot_0001.csv
│   ├── ...
│   └── energy.dat
├── b4pc_gdisph/
├── b5pc_gdisph/
├── b6pc_gdisph/
├── plots/
│   ├── density_evolution.png
│   ├── axis_ratios.png
│   └── phase_diagram.png
└── animations/
    └── disruption_animation.gif
```

## Diagnostics

### 1. Morphology

- **Axis ratios** (a:b:c) from moment of inertia tensor
- **Elongation**: e = 1 - c/a (0 = spherical, 1 = highly elongated)

### 2. Energetics

- **Kinetic energy**: E_kin = Σ (1/2) m v²
- **Thermal energy**: E_thermal = Σ m u
- **Gravitational**: E_grav = -Σ m Φ (cloud + BH)
- **Total energy**: Should be conserved (< 1% drift)

### 3. Mass Budget

- **Bound fraction**: Particles with E_total < 0
- **Unbound fraction**: Particles with E_total > 0
- **Accretion**: Mass within r < r_acc ~ 0.01 pc

### 4. Thermal State

- **n-T phase diagram**: Should follow K&I equilibrium curve
- **Cooling time**: τ_cool should remain << t_tidal

## Code Units

- **Length**: 1 code unit = 1 parsec (pc)
- **Mass**: 1 code unit = 1 M_☉
- **Time**: 1 code unit ~ 0.978 Myr
- **Velocity**: 1 code unit ~ 1.02 km/s
- **Density**: ρ [code] → n_H = ρ × 40.5 cm⁻³

## Makefile Targets

```bash
# Single runs
imbh_b3pc_run          # Run b = 3 pc
imbh_b6pc_run          # Run b = 6 pc

# Parameter scans
imbh_scan_run          # Run b = 3, 4, 5, 6 pc

# Visualization
imbh_visualize         # Generate all plots

# Cleanup
imbh_clean             # Remove all results

# Help
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_help
```

## Implementation Status

### ✅ Completed

- [x] Research setup documentation
- [x] External force module (point-mass BH)
- [x] Initial conditions (imbh_cloud.cpp)
- [x] Configuration presets
- [x] Makefile
- [x] Visualization scripts
- [x] README

### ⏳ Pending Integration

- [ ] Add external force to CMakeLists.txt
- [ ] Extend SPHParameters for external_force section
- [ ] Integrate imbh_cloud into solver.cpp sample dispatcher
- [ ] Unit tests for external force module
- [ ] Validate energy conservation with BH force
- [ ] Test thermal equilibrium tracking

## References

1. **Koyama & Inutsuka (2000)**, ApJ 532, 980  
   *Molecular Cloud Formation in Shock-Compressed Layers*

2. **Guillochon & Loeb (2015)**, ApJ 811, 20  
   *Tidal Disruption Events: Physics and Observables*

3. **Burkert & Naab (2013)**, MNRAS 434, 36  
   *Cloud-Black Hole Interactions in Galactic Nuclei*

4. **Hopkins (2013)**, MNRAS 428, 2840  
   *GDISPH Formulation*

## Next Steps

1. **Complete integration**: Add to CMakeLists.txt, extend parameters
2. **Test suite**: Energy conservation, thermal equilibrium validation
3. **Convergence study**: N = 5×10⁴, 2×10⁵, 10⁶ particles
4. **Parameter exploration**: b = 2-10 pc, M_BH = 10⁴-10⁶ M_☉
5. **Physical extensions**:
   - Magnetic fields (MHD)
   - Star formation triggered by compression
   - Multiple clouds + BH binary

## Contact

For questions about this setup, refer to `docs/RESEARCH_SETUP.md` for detailed physics and resolution analysis.
