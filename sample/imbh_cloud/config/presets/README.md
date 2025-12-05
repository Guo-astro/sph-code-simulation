# IMBH-Cloud Simulation Presets

Configuration presets for Lane-Emden polytropic sphere simulations with IMBH tidal disruption.

## Directory Structure

```
config/presets/
├── README.md                    # This file
├── relaxation/                  # Step 1: Generate relaxed initial conditions
│   ├── scenarios/
│   │   ├── 10k/                 # ~10,648 particles (N=22)
│   │   │   ├── gsph.json
│   │   │   └── gsph.json
│   │   ├── 61k/                 # ~64,000 particles (N=40) - 5% accuracy target
│   │   │   ├── gsph.json
│   │   │   └── gsph.json
│   │   └── 200k/                # ~200,000 particles (N=100) - production
│   │       ├── gsph.json
│   │       └── gsph.json
│   └── README.md
├── hydrostatic/                 # Step 2: Test hydrostatic equilibrium
│   ├── scenarios/
│   │   ├── 10k/
│   │   ├── 61k/
│   │   └── 200k/
│   ├── PARTICLE_COUNT_DERIVATION.md
│   └── README.md
└── simulation/                  # Step 3: IMBH encounter simulations
    └── scenarios/
        ├── b1p5pc_optimal/      # Impact parameter 1.5 pc
        ├── b3pc_cool/           # Impact parameter 3 pc with cooling
        ├── b3pc_nocool/         # Impact parameter 3 pc adiabatic
        └── b6pc_nocool/         # Impact parameter 6 pc adiabatic
```

## Complete Workflow

### Overview

The simulation workflow consists of three stages:

1. **Relaxation**: Generate relaxed initial conditions from Lane-Emden polytrope
2. **Hydrostatic Test**: Verify equilibrium stability with self-gravity
3. **IMBH Simulation**: Run tidal disruption encounter (optional)

### Step 1: Relaxation (Generate Initial Conditions)

Generate a relaxed Lane-Emden n=1.5 polytropic sphere using analytical force relaxation.

```bash
# Quick test (~2-5 min, 10k particles)
make -f sample/imbh_cloud/Makefile.relaxation relax_oneshot SIZE=10k METHOD=gsph

# High resolution for 5% accuracy (~30-60 min, 61k particles)
make -f sample/imbh_cloud/Makefile.relaxation relax_oneshot SIZE=61k METHOD=gsph

# Production run (~4-8 hours, 200k particles)
make -f sample/imbh_cloud/Makefile.relaxation relax_oneshot SIZE=200k METHOD=gsph
```

**Output**: `sample/imbh_cloud/results/relaxation/<SIZE>/<METHOD>/snapshot_NNNN.csv`

**Physics**:
- Lane-Emden polytrope with n=1.5 (γ=5/3)
- ξ₁ = 3.6537540101 (first zero of θ)
- Central density ρ_c = 1.43 [code units]
- Analytical force relaxation (gravity OFF during relaxation)

### Step 2: Hydrostatic Equilibrium Test

Verify that the relaxed sphere maintains hydrostatic equilibrium with self-gravity enabled.

```bash
# Test 10k particles
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_oneshot SIZE=10k METHOD=gsph

# Test 61k particles (5% accuracy)
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_oneshot SIZE=61k METHOD=gsph

# Test 200k particles (production)
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_oneshot SIZE=200k METHOD=gsph
```

**Input**: Loads relaxed snapshot via `resumeFromSnapshot` parameter

**Output**: `sample/imbh_cloud/results/hydrostatic/<SIZE>/<METHOD>/`

**Stability Criteria**:
- Density RMS error < 5%
- Velocity |v| < 1% sound speed
- Energy drift < 1%

### Step 3: IMBH Encounter Simulation (Optional)

Run tidal disruption simulation with intermediate-mass black hole.

```bash
# Standard b=3pc adiabatic
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_oneshot SCENARIO=b3pc_nocool METHOD=gsph

# Close encounter b=1.5pc
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_oneshot SCENARIO=b1p5pc_optimal METHOD=gsph

# With cooling
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_oneshot SCENARIO=b3pc_cool METHOD=gsph
```

## Size Options

| SIZE | N | Particles | Relaxation Steps | Expected Accuracy | Runtime |
|------|---|-----------|------------------|-------------------|---------|
| 10k  | 22 | ~10,648  | 400,000          | ~16% on ρ_c       | 2-5 min |
| 61k  | 40 | ~64,000  | 200,000          | ~5% on ρ_c        | 30-60 min |
| 200k | 100| ~200,000 | 100,000          | ~2% on ρ_c        | 4-8 hours |

**Accuracy Derivation**: See `hydrostatic/PARTICLE_COUNT_DERIVATION.md`

SPH error scales as: ε ~ N^(-2/3)

## SPH Methods

| Method | Description | Best For |
|--------|-------------|----------|
| `gsph` | Godunov DISPH (Riemann + density-independent) | General use, recommended |
| `gsph` | Godunov SPH (pure Riemann solver) | Comparison, no AV parameters |

## Quick Reference Commands

```bash
# Show help for each stage
make -f sample/imbh_cloud/Makefile.relaxation relax_help
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_help
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_help

# Visualization only (after simulation)
make -f sample/imbh_cloud/Makefile.relaxation relax_oneshot_viz SIZE=10k METHOD=gsph
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_oneshot_viz SIZE=10k METHOD=gsph

# Clean results
make -f sample/imbh_cloud/Makefile.relaxation relax_oneshot_clean SIZE=10k METHOD=gsph
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_oneshot_clean SIZE=10k METHOD=gsph

# List available configs
make -f sample/imbh_cloud/Makefile.relaxation relax_list_configs
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_list_configs
```

## Config File Parameters

### Relaxation Config Keys

| Parameter | Description | Example |
|-----------|-------------|---------|
| `N` | Grid resolution (N³ particles) | 22, 40, 100 |
| `sample` | Initial condition generator | `"lane_emden"` |
| `relaxationSteps` | Number of relaxation steps | 400000 |
| `relaxationOnly` | Stop after relaxation | `true` |
| `useGravity` | Gravity during relaxation | `false` |

### Hydrostatic Config Keys

| Parameter | Description | Example |
|-----------|-------------|---------|
| `resumeFromSnapshot` | Path to relaxed IC | `"results/relaxation/10k/GDISPH/snapshot_0020.csv"` |
| `resetTimeOnResume` | Reset time to 0 | `true` |
| `useGravity` | Enable self-gravity | `true` |
| `G` | Gravitational constant | 1.0 |
| `endTime` | Simulation end time | 100.0 |

## Lane-Emden Physics Reference

For polytropic index n=1.5 (γ=5/3):

| Parameter | Symbol | Value |
|-----------|--------|-------|
| First zero | ξ₁ | 3.6537540101 |
| Density ratio | ρ_c/ρ̄ | 5.99071 |
| Central density | ρ_c | 1.43 [code] |
| Expected K | K = P/ρ^γ | 0.4242 |

## Troubleshooting

### "Initial condition not found"
Run relaxation first:
```bash
make -f sample/imbh_cloud/Makefile.relaxation relax_oneshot SIZE=10k METHOD=gsph
```

### "DIM=3 required"
Rebuild with 3D:
```bash
cd build && cmake -DDIM=3 .. && make -j
```

### Central density error > expected
- Check K(r) profile for constancy (isentropic IC)
- Increase particle count for better resolution
- Verify `lane_emden.cpp` uses analytic density for internal energy

## References

- Lane-Emden equation: Chandrasekhar (1939)
- SPH convergence: Price (2012)
- GDISPH method: Saitoh & Makino (2013)
