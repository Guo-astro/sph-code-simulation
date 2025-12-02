# IMBH Cloud Encounter Scenarios

This directory contains pre-configured JSON files for different IMBH-cloud encounter scenarios. Each scenario represents a different set of physical parameters (impact parameter, cooling).

## Directory Structure

```
scenarios/
├── b3pc_nocool/           # Standard case: b=3pc, adiabatic
│   ├── gdisph.json        # Godunov DISPH method
│   └── gsph.json          # Godunov SPH method
├── b3pc_cool/             # With cooling: b=3pc, radiative cooling enabled
│   ├── gdisph.json
│   └── gsph.json
├── b6pc_nocool/           # Wide orbit: b=6pc, adiabatic
│   ├── gdisph.json
│   └── gsph.json
└── b1p5pc_optimal/        # Close encounter: b=1.5pc, adiabatic
    ├── gdisph.json
    └── gsph.json
```

## Scenarios Overview

### b3pc_nocool (Standard Case)
- **Impact parameter**: 3 pc
- **Cooling**: Disabled (adiabatic)
- **Description**: Standard reference case for tidal disruption studies

### b3pc_cool (With Cooling)
- **Impact parameter**: 3 pc
- **Cooling**: Enabled (radiative cooling)
- **Description**: Tests effect of radiative cooling on disruption dynamics

### b6pc_nocool (Wide Orbit)
- **Impact parameter**: 6 pc
- **Cooling**: Disabled (adiabatic)
- **Description**: Weaker tidal interaction, cloud passes further from BH

### b1p5pc_optimal (Close Encounter)
- **Impact parameter**: 1.5 pc
- **Cooling**: Disabled (adiabatic)
- **Description**: Stronger tidal disruption, closer to optimal for stream formation

## SPH Methods

### gdisph (Godunov DISPH)
- Godunov-type Riemann solver combined with density-independent SPH
- Uses `avAlpha=0.0` (pure Riemann, no artificial viscosity)
- Balsara switch enabled for shock detection
- Best for: discontinuities, contact surfaces, tidal streams

### gsph (Godunov SPH)
- Pure Godunov SPH with Riemann solver
- No artificial viscosity parameters (dissipation fully handled by Riemann solver)
- Best for: strong shocks, blast waves, clean shock tube problems

## Usage

### Using Makefile (Recommended)

```bash
# Default scenario (b3pc_nocool with gdisph)
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_oneshot

# Specify scenario and method
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_oneshot SCENARIO=b3pc_nocool METHOD=gdisph
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_oneshot SCENARIO=b3pc_nocool METHOD=gsph

# Or use convenience targets
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_b3pc_gdisph
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_b3pc_gsph
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_b1p5pc_gdisph
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_cool_gdisph
```

### Direct Execution

```bash
./build/sph sample/imbh_cloud/config/presets/simulation/scenarios/b3pc_nocool/gdisph.json
```

## Output Structure

Results are organized by scenario and method:

```
results/
├── b3pc_nocool/
│   ├── GDISPH/
│   │   ├── snapshot_0000.csv
│   │   ├── snapshot_0001.csv
│   │   └── ...
│   └── GSPH/
│       └── ...
├── b3pc_cool/
│   └── ...
└── ...
```

Animations are similarly organized:

```
animations/
├── b3pc_nocool/
│   ├── GDISPH/
│   │   ├── tidal.gif
│   │   ├── shock_diagnostics.gif
│   │   ├── shock_diagnostics_com.gif
│   │   ├── density.gif
│   │   ├── temperature.gif
│   │   ├── mach.gif
│   │   └── shock.gif
│   └── GSPH/
└── ...
```

## Physical Parameters

All scenarios share these base parameters:

| Parameter | Value | Description |
|-----------|-------|-------------|
| IMBH mass | 10⁵ M☉ | 100 code units |
| Cloud mass | ~1000 M☉ | 10,648 particles |
| Cloud radius | ~1 pc | Relaxed Lane-Emden n=1.5 |
| Approach velocity | 10 km/s | Code units |
| Unit system | IMBH_ENCOUNTER | [L]=1pc, [M]=1000M☉, [V]=1km/s |
| Simulation time | 4.0 | ~3.9 Myr |

## Creating New Scenarios

To add a new scenario:

1. Create a new directory: `scenarios/<scenario_name>/`
2. Copy an existing config as template
3. Modify the `initialCondition.transform.translate` to set the impact parameter
4. Modify `thermal.cooling` as needed
5. Update `outputDirectory` and `checkpointPath` to match the new scenario name
