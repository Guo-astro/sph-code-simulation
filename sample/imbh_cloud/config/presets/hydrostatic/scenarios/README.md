# Hydrostatic Equilibrium Test Scenarios

This directory contains pre-configured JSON files for hydrostatic equilibrium tests.

## Purpose

Verify that relaxed Lane-Emden spheres maintain hydrostatic equilibrium with self-gravity enabled. This is a crucial validation step before using the initial conditions in IMBH encounter simulations.

## Directory Structure

```
scenarios/
├── 10k/                        # Testing: ~10,648 particles
│   ├── gdisph_short.json       # GDISPH, t=100 (quick test)
│   ├── gdisph_long.json        # GDISPH, t=1000 (stability)
│   ├── gsph_short.json         # GSPH, t=100
│   └── gsph_long.json          # GSPH, t=1000
└── 200k/                       # Production: ~200,000 particles
    ├── gdisph_short.json
    ├── gdisph_long.json
    ├── gsph_short.json
    └── gsph_long.json
```

## Prerequisites

**Requires relaxed initial conditions from Makefile.relaxation:**
```bash
make -f sample/imbh_cloud/Makefile.relaxation relax_10k
```

## Duration Options

### short (Quick Test)
- **End time**: t = 100 code units
- **Output interval**: dt = 1.0 (~100 snapshots)
- **Runtime**: ~5-10 minutes (10k), ~1-2 hours (200k)
- **Purpose**: Quick verification that IC is stable

### long (Stability Test)
- **End time**: t = 1000 code units
- **Output interval**: dt = 10.0 (~100 snapshots)
- **Runtime**: ~30-60 minutes (10k), ~8-12 hours (200k)
- **Purpose**: Extended stability verification

## SPH Methods

### gdisph
- Godunov DISPH with pure Riemann solver (`avAlpha=0`)
- Balsara switch enabled
- Recommended for production simulations

### gsph
- Godunov SPH (pure Riemann solver)
- No artificial viscosity parameters
- Good for comparison studies

## Usage

```bash
# Quick test (10k, short)
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_10k

# Long stability test
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_10k_long

# With method selection
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_oneshot SIZE=10k METHOD=gdisph DURATION=short
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_oneshot SIZE=10k METHOD=gsph DURATION=long

# Generate stability report
make -f sample/imbh_cloud/Makefile.hydrostatic hydro_stability_report SIZE=10k METHOD=gdisph DURATION=short
```

## Output Structure

```
results/hydrostatic/
├── 10k/
│   ├── GDISPH_short/
│   │   ├── snapshot_0000.csv
│   │   └── ...
│   ├── GDISPH_long/
│   ├── GSPH_short/
│   └── GSPH_long/
└── 200k/
    └── ...
```

## Stability Criteria

| Criterion | Threshold | Description |
|-----------|-----------|-------------|
| Energy drift | < 1% | Total energy should be conserved |
| Density RMS | < 5% | Density profile should remain stable |
| Velocity ratio | < 1% c_s | Velocities should remain small |

### Interpretation

- **✅ PASS**: All criteria met → IC ready for IMBH simulations
- **⚠️ MARGINAL**: Minor deviations → May need more relaxation
- **❌ FAIL**: Large deviations → Sphere is unstable, check relaxation

## Physics Notes

| Timescale | Value | Description |
|-----------|-------|-------------|
| t_cross | ~1.0 | Sound crossing time (R/c_s) |
| t_dyn | ~1.0 | Dynamical time (√(R³/GM)) |

A simulation of t=100 tests ~100 dynamical times.
A simulation of t=1000 tests ~1000 dynamical times.

## Next Steps

After successful hydrostatic test:
```bash
# Run IMBH encounter simulation
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_oneshot SCENARIO=b3pc_nocool METHOD=gdisph
```
