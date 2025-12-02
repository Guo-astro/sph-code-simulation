# Lane-Emden Relaxation Scenarios

This directory contains pre-configured JSON files for different relaxation scenarios.

## Directory Structure

```
scenarios/
├── 10k/                    # Testing: ~10,648 particles
│   ├── gdisph.json         # GDISPH method (recommended)
│   └── gsph.json           # GSPH method (pure Riemann)
└── 200k/                   # Production: ~200,000 particles
    ├── gdisph.json
    └── gsph.json
```

## Size Options

### 10k (Testing)
- **Particles**: ~10,648 (N=22, formula: N³/5)
- **Steps**: 400,000 relaxation steps
- **Output**: Every 1,000 steps (~400 snapshots)
- **Runtime**: ~2-5 minutes
- **Use case**: Quick testing, parameter tuning

### 200k (Production)
- **Particles**: ~200,000 (N=100)
- **Steps**: 100,000 relaxation steps
- **Output**: Every 2,000 steps (~50 snapshots)
- **Runtime**: ~4-8 hours
- **Use case**: Production IMBH encounter simulations

## SPH Methods

### gdisph (Recommended)
- Godunov DISPH (Riemann solver + density-independent formulation)
- Uses `avAlpha=1.0` for smooth relaxation
- Balsara switch enabled
- Best for: stable relaxation with good energy conservation

### gsph
- Godunov SPH (pure Riemann solver)
- No artificial viscosity parameters
- Can be noisier during relaxation
- Best for: comparison studies

## Usage

```bash
# Using Makefile (recommended)
make -f sample/imbh_cloud/Makefile.relaxation relax_10k
make -f sample/imbh_cloud/Makefile.relaxation relax_200k

# Or with method selection
make -f sample/imbh_cloud/Makefile.relaxation relax_oneshot SIZE=10k METHOD=gdisph
make -f sample/imbh_cloud/Makefile.relaxation relax_oneshot SIZE=10k METHOD=gsph

# Direct execution
./build/sph sample/imbh_cloud/config/presets/relaxation/scenarios/10k/gdisph.json
```

## Output Structure

```
results/relaxation/
├── 10k/
│   ├── GDISPH/
│   │   ├── snapshot_0000.csv
│   │   ├── snapshot_0001.csv
│   │   └── ...
│   │   └── snapshot_0400.csv  ← Use this as IC
│   └── GSPH/
│       └── ...
└── 200k/
    ├── GDISPH/
    └── GSPH/
```

## Lane-Emden Parameters

All scenarios use these physical parameters:

| Parameter | Value | Description |
|-----------|-------|-------------|
| Polytropic index | n = 1.5 | γ = 5/3 for monoatomic gas |
| ξ₁ | 3.6537540101 | First zero of Lane-Emden solution |
| Central density | ρ_c = 1.43 | Code units |
| Total mass | M = 1.0 | Code units |
| Radius | R = 1.0 | Code units |

## Quality Criteria

A successful relaxation should show:
- Velocities decaying to near zero (< 0.01 code units)
- Density profile matching analytic Lane-Emden solution
- Energy approximately conserved

## Next Steps

After successful relaxation:
1. Run hydrostatic equilibrium test:
   ```bash
   make -f sample/imbh_cloud/Makefile.hydrostatic hydro_10k
   ```
2. Use final snapshot as IC for IMBH encounter:
   ```bash
   make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_oneshot
   ```
