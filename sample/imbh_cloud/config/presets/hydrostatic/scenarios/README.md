# Hydrostatic Equilibrium Test Scenarios

This directory contains pre-configured JSON files for hydrostatic equilibrium tests.

## Purpose

Verify that relaxed Lane-Emden spheres maintain hydrostatic equilibrium with self-gravity enabled. This is a crucial validation step before using the initial conditions in IMBH encounter simulations.

## Directory Structure (SSOT)

```
scenarios/
├── 10k/                        # Testing: ~10,648 particles
│   ├── gdisph.json             # GDISPH method
│   └── gsph.json               # GSPH method
└── 200k/                       # Production: ~200,000 particles
    ├── gdisph.json
    └── gsph.json
```

Each JSON file is the **Single Source of Truth (SSOT)** for that size/method combination.

## Prerequisites

**Requires relaxed initial conditions from Makefile.relaxation:**
```bash
make -f sample/imbh_cloud/Makefile.relaxation relax_oneshot SIZE=10k METHOD=gdisph
```

## Test Parameters

- **End time**: t = 100 code units (~few dynamical times)
- **Output interval**: dt = 1.0 (~100 snapshots)
- **Runtime**: ~5-10 minutes (10k), ~1-2 hours (200k)

## SPH Methods

### gdisph
- Godunov DISPH with pure Riemann solver (`avAlpha=0`)
- Balsara switch enabled
- Density-independent formulation

### gsph
- Godunov SPH with pure Riemann solver
- No artificial viscosity at all
- All dissipation from Riemann solver

## Usage

```bash
# Run hydrostatic test
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_oneshot SIZE=10k METHOD=gdisph
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_oneshot SIZE=10k METHOD=gsph

# Visualization
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_oneshot_viz SIZE=10k METHOD=gdisph

# Stability report
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_stability_report SIZE=10k METHOD=gdisph
```

## Stability Criteria

A hydrostatic test **passes** if:
- **Density RMS change** < 5%
- **Max velocity** < 1% of sound speed
- **Energy drift** < 1%

## Important: Snapshot Naming

The config files reference `snapshot_final.csv` as the initial condition. After running relaxation, you may need to update this path to match the actual final snapshot number (e.g., `snapshot_0024.csv`).

Check the actual file:
```bash
ls sample/imbh_cloud/results/relaxation/10k/GSPH/snapshot_*.csv | tail -1
```
