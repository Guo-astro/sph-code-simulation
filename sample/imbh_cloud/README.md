# IMBH-Molecular Cloud Tidal Disruption Simulation

## Overview

This directory contains simulations of intermediate-mass black hole (IMBH) tidal disruption of molecular clouds using 3D Lane-Emden polytropic spheres.

**Physics:**
- Lane-Emden n=1.5 polytrope (γ=5/3)
- GDISPH method with Wendland C4 kernel
- Koyama & Inutsuka (2000) thermal equilibrium
- Impact parameters: b = 3-6 pc
- Relative velocities: v = 0-20 km/s
- Black hole mass: M_BH ~ 10^5 M☉

## Quick Start

All commands must be run from the **repository root** (`/Users/guo/Downloads/sphcode/`):

```bash
cd /Users/guo/Downloads/sphcode
```

### 1. Lane-Emden Relaxation (Required First Step)

**Run 20,000-step relaxation (~10k particles):**
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax
```

**Production run (200k particles, 100k steps - HIGH RESOLUTION):**
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_200k
```
**Warning:** This takes 4-8 hours and requires ~5-10 GB disk space.

**Resume from checkpoint:**
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_resume \
     SNAPSHOT=sample/imbh_cloud/results/lane_emden_50k_relax/snapshot_0021.csv \
     STEPS=10000 \
     FREQ=1000
```

**Visualize existing results:**
```bash
# For 10k particle relaxation
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_viz

# For 200k particle relaxation
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_200k_viz
```

### 2. Hydrostatic Self-Gravity Test (NEW! - Required Before IMBH Runs)

**Purpose:** Verify self-gravity implementation by testing equilibrium stability over multiple crossing times.

**Quick test (2k particles, ~5-10 minutes):**
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_2k
```

**Production test (200k particles, ~8-16 hours):**
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_200k
```

**What it tests:**
- ✓ Self-gravity correctly balances pressure forces
- ✓ Density profile remains constant over 10 crossing times
- ✓ Velocities stay << 1% sound speed (no spurious motion)
- ✓ Energy conserved to < 1% (numerical stability)
- ✓ Test duration >> tidal disruption timescale

**Timescales (in code units):**
- Crossing time: t_cross ~ 1.0 (sound crossing time)
- Dynamical time: t_dyn = 1.0 (free-fall time)
- Tidal time: t_tidal ~ 3.2 (for M_BH/M_cloud = 10)
- Test duration: 10.0 code units (>> all relevant timescales)

**Visualize results:**
```bash
# For 2k test
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_2k_viz

# For 200k test
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_200k_viz
```

**Pass criteria:**
- RMS density error < 5% from analytic Lane-Emden profile
- Max velocity < 1% of sound speed
- Median velocity < 0.1% of sound speed
- Energy drift < 1% over entire test duration

### 3. IMBH Disruption Simulations

**Single run (b=3pc, v=10km/s):**
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_b3pc_run
```

**Impact parameter scan (b=3,4,5,6 pc):**
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_scan_impact
```

**Velocity scan (v=0,5,10,15,20 km/s):**
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_scan_velocity
```

### 3. Visualization

**Create relaxation animation:**
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_viz
```

**View animation:**
```bash
open sample/imbh_cloud/results/lane_emden_50k_relax/relaxation_animation.gif
```

## Directory Structure

```
sample/imbh_cloud/
├── README.md                    # This file
├── Makefile.imbh_cloud         # Build and run automation
├── config/
│   └── presets/
│       ├── lane_emden_2k_relax.json           # 2k relaxation config
│       ├── lane_emden_200k_relax.json         # 200k relaxation config
│       ├── lane_emden_2k_hydrostatic.json     # 2k hydrostatic test (NEW!)
│       ├── lane_emden_200k_hydrostatic.json   # 200k hydrostatic test (NEW!)
│       ├── imbh_cloud_b3pc_gdisph.json        # b=3pc preset
│       ├── imbh_cloud_b4pc_gdisph.json        # b=4pc preset
│       └── ...                                 # Other presets
├── results/
│   ├── lane_emden_2k_relax/        # 2k relaxation outputs
│   ├── lane_emden_200k_relax/      # 200k relaxation outputs
│   ├── lane_emden_2k_hydrostatic/  # 2k hydrostatic test (NEW!)
│   ├── lane_emden_200k_hydrostatic/# 200k hydrostatic test (NEW!)
│   │   ├── snapshot_*.csv          # Particle data (50 snapshots)
│   │   ├── energy.dat              # Energy evolution
│   │   ├── hydrostatic_animation.gif    # 6-panel animation
│   │   ├── hydrostatic_summary.png      # Final state analysis
│   │   └── run.log                 # Simulation log
│   └── disruption_runs/            # IMBH disruption results
└── scripts/
    ├── visualize_relaxation.py     # Relaxation animation & plotting
    ├── visualize_hydrostatic.py    # Hydrostatic test visualization (NEW!)
    └── ...                         # Other analysis scripts
```

## Workflow: Complete IMBH Simulation Pipeline

### Step 1: Lane-Emden Relaxation (Create equilibrium initial conditions)

```bash
# From repository root
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_2k
```

**Output:**
- 400+ snapshots with relaxation steps
- Location: `sample/imbh_cloud/results/lane_emden_2k_relax/`

### Step 2: Hydrostatic Self-Gravity Test (Verify gravity implementation)

```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_2k
```

**Purpose:**
- Verify self-gravity balances pressure correctly
- Ensure no spurious forces or numerical instabilities
- Test stability over multiple crossing times (>> tidal timescale)

**Output:**
- 50 snapshots over 10 crossing times
- Location: `sample/imbh_cloud/results/lane_emden_2k_hydrostatic/`
- Animation shows density profile stability and velocity evolution

**Pass Criteria:**
- ✓ Density RMS error < 5% from Lane-Emden analytic
- ✓ Velocities < 1% sound speed
- ✓ Energy conserved to < 1%

### Step 3: IMBH Tidal Disruption (Main science simulation)

```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_b3pc_run
```

**Uses verified self-gravity from Step 2**

### Production Workflow (200k particles)

```bash
# Step 1: High-resolution relaxation (4-8 hours)
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_200k

# Step 2: Production hydrostatic test (8-16 hours)
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_200k

# Step 3: IMBH disruption simulations (ready for production runs)
# ... proceed with science simulations
```

## Workflow: Resume from Checkpoint

```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_resume \
     SNAPSHOT=sample/imbh_cloud/results/lane_emden_50k_relax/snapshot_0021.csv \
     STEPS=10000 \
     FREQ=1000 \
     OUTDIR=sample/imbh_cloud/results/resumed_relax
```

**Parameters:**
- `SNAPSHOT` - Path to checkpoint file (required)
- `STEPS` - Additional relaxation steps (default: 10000)
- `FREQ` - Output frequency (default: 500)
- `OUTDIR` - Output directory (default: results/resumed_relax)

### Step 3: Create Visualization

```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_viz
```

**Generated files:**
- `relaxation_animation.gif` - 6-panel animation with Lane-Emden analytic overlay
- `relaxation_summary.png` - Final state summary
- `velocity_decay.png` - Velocity evolution analysis

## Lane-Emden Parameters

**Polytrope Configuration:**
- n = 1.5 (polytropic index)
- γ = 5/3 (adiabatic index)
- ξ₁ = 3.6537540101 (first zero)
- α = R/ξ₁ = 0.273691 (scaling factor)

**Initial Conditions:**
- N_shells = 22 → ~2,130 particles (testing/quick runs)
- N_shells = 100 → 200,000 particles (production runs)
- ρ_c = 1.43009692 [code units]
- K = 0.424209 (polytropic constant)
- R = 1.0, M_total = 1.0 [code units]

**Particle Scaling:**
- Formula: N_particles = N³/5 (for 3D Lane-Emden)
- N=22: ~2,130 particles (fast testing)
- N=63: ~50,000 particles (medium resolution)
- N=100: 200,000 particles (production quality)
- N=158: ~1,000,000 particles (very high resolution)

**Numerical Settings:**
- SPH method: GDISPH
- Kernel: Wendland C4
- Neighbor count: 50
- Riemann solver: HLL
- Artificial viscosity: α=1.0 with Balsara switch

## Snapshot Format

CSV files contain metadata header followed by particle data:

**Header sections (61 lines):**
1. Simulation metadata (code version, dimensions, SPH method)
2. Physical parameters (gamma, kernel type)
3. Relaxation metadata (step number, Lane-Emden parameters)
4. Column descriptions
5. Column names

**Data columns:**
```
id, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, 
acc_x, acc_y, acc_z, mass, dens, pres, ene, 
sml, sound, alpha, balsara, gradh, phi, neighbor, is_ghost
```

## Resume Functionality

The resume feature preserves:
- ✓ Lane-Emden parameters (α, ρ_c, K, R, M)
- ✓ Particle positions, velocities, masses
- ✓ Snapshot counter (no file overwrites)
- ✓ Relaxation step number

**Example resume chain:**
```bash
# Initial run: 20,000 steps
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax

# Resume #1: +10,000 steps (total: 30,000)
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_resume \
     SNAPSHOT=sample/imbh_cloud/results/lane_emden_50k_relax/snapshot_0021.csv \
     STEPS=10000

# Resume #2: +10,000 steps (total: 40,000)
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_resume \
     SNAPSHOT=sample/imbh_cloud/results/resumed_relax/snapshot_0011.csv \
     STEPS=10000
```

## Cleanup

**Remove relaxation results only:**
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_clean
```

**Remove all IMBH results:**
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_clean
```

## Help

**View all available targets:**
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_help
```

## Troubleshooting

### "No rule to make target 'imbh_relax'"

**Cause:** Running from wrong directory or missing `-f` flag.

**Solution:** Always use full command from repo root:
```bash
cd /Users/guo/Downloads/sphcode
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax
```

### "Error loading snapshot: could not convert string 'id' to float64"

**Cause:** Visualization script not handling CSV metadata header.

**Solution:** Already fixed in `scripts/visualize_relaxation.py` (dynamic header detection).

### "Snapshot counter overwrites existing files"

**Cause:** Resume not setting counter based on snapshot step number.

**Solution:** Already fixed in `src/solver.cpp` lines 930-937.

### Dimension mismatch error

**Cause:** Binary compiled with wrong dimension (DIM=2 instead of DIM=3).

**Solution:**
```bash
cd build
cmake -DSPH_DIM=3 ..
make -j8
```

## Configuration Files

**Main relaxation config:**
`config/presets/lane_emden_50k_relax.json`

**Key parameters to adjust:**
```json
{
  "N": 22,                        // N_shells (22³ ≈ 10,648 particles)
  "relaxationSteps": 20000,       // Total relaxation steps
  "relaxationOutputFreq": 1000,   // Snapshot frequency
  "alpha_scaling": 0.273691,      // Lane-Emden α parameter
  "rho_center": 1.43009692,       // Central density
  "K": 0.424209,                  // Polytropic constant
  "neighborNumber": 50,           // SPH neighbors
  "kernel": "wendland"            // Kernel function
}
```

## Quality Checks

**After relaxation, verify:**

1. **Energy conservation:**
   ```bash
   tail -1 sample/imbh_cloud/results/lane_emden_50k_relax/energy.dat
   ```
   Drift should be < 0.01%

2. **Velocity decay:**
   - View `velocity_decay.png`
   - Final velocity should approach zero

3. **Density profile:**
   - Check `relaxation_animation.gif`
   - SPH density should match Lane-Emden analytic curve

4. **Central density:**
   ```bash
   python3 -c "import numpy as np; data = np.loadtxt('sample/imbh_cloud/results/lane_emden_50k_relax/snapshot_0021.csv', delimiter=',', skiprows=61); r = np.sqrt(np.sum(data[:, 1:4]**2, axis=1)); rho = data[:, 11]; central = rho[r < 0.05]; print(f'Central density: {central.mean():.3f} (analytic: 1.43)')"
   ```

## Expected Physics (IMBH Tidal Disruption)

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

### Key Timescales

- **Tidal timescale**: t_tidal ~ √(R³/(GM_BH)) ~ 10⁴-10⁵ yr
- **Crossing time**: t_cross ~ R/v_rel ~ 10⁵-10⁶ yr  
- **Thermal time**: t_thermal ~ 10³-10⁴ yr << t_tidal ✓

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


## Resolution Analysis

| N (shells) | Formula    | Total Particles | Use Case              |
|------------|------------|----------------|-----------------------|
| 22         | 22³/5      | ~2,130         | Quick testing         |
| 40         | 40³/5      | ~12,800        | Medium resolution     |
| 63         | 63³/5      | ~50,000        | Good quality          |
| 100        | 100³/5     | 200,000        | **Production (recommended)** |
| 158        | 158³/5     | ~1,000,000     | Very high resolution  |

**Recommendation:** 
- **Testing**: N=22 (minutes)
- **Production**: N=100 (hours, 200k particles - matches your requirement!)
- **Publication**: N≥100 (ensures N_J > 50-100 particles per Jeans mass)

## Physical Units (Code Units)

- **Length**: [R] = 1.0 (normalized)
- **Mass**: [M] = 1.0 (normalized)
- **Density**: [ρ_c] = 1.43009692 (Lane-Emden central density)
- **Time**: [t] = √(R³/(GM)) (dynamical time)

For physical molecular cloud:
- Scale to R = 5 pc, M = 10⁴ M_☉
- ρ_physical = ρ_code × (M/(4/3 π R³))

## Citation

If using this code for research, please cite:

- **Lane-Emden solution:** Chandrasekhar (1939)
- **GDISPH method:** Hopkins (2013), MNRAS 428, 2840
- **Thermal equilibrium:** Koyama & Inutsuka (2000), ApJ 532, 980

## Support

For issues or questions:
1. Check `make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_help`
2. Review this README troubleshooting section
3. Check simulation logs in `results/*/run.log`

---

**Last updated:** 2025-12-01  
**Author:** Guo  
**Repository:** sph-code-simulation

