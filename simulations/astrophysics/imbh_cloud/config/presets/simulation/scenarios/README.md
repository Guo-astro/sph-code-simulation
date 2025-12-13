# IMBH Cloud Encounter Scenarios

Pre-configured JSON files for IMBH-cloud tidal disruption simulations.

## Directory Structure

```
scenarios/
├── Mc1e3_Mbh1e5_b1p5_v10/          # Impact parameter 1.5 pc (deep penetration)
│   ├── adiabatic_61k_gsph.json
│   ├── adiabatic_61k_gdisph.json
│   ├── radiative_61k_gsph.json
│   └── radiative_61k_gdisph.json
├── Mc1e3_Mbh1e5_b3_v10/            # Impact parameter 3.0 pc (moderate)
│   ├── adiabatic_61k_gsph.json
│   ├── adiabatic_61k_gdisph.json
│   ├── radiative_61k_gsph.json
│   └── radiative_61k_gdisph.json
├── Mc1e3_Mbh1e5_b6_v10/            # Impact parameter 6.0 pc (grazing)
│   ├── adiabatic_61k_gsph.json
│   ├── adiabatic_61k_gdisph.json
│   ├── radiative_61k_gsph.json
│   └── radiative_61k_gdisph.json
├── README.md                        # This file
└── SIMULATION_MATRIX.md             # Complete physics documentation
```

## Naming Convention

### Folder Names: `Mc{cloud}_Mbh{bh}_b{impact}_v{vel}/`

| Component | Meaning | Example |
|-----------|---------|---------|
| `Mc1e3` | Cloud mass = 10³ M☉ | 1.0 code units |
| `Mbh1e5` | BH mass = 10⁵ M☉ | 100.0 code units |
| `b1p5` | Impact parameter = 1.5 pc | 1.5 code units |
| `v10` | Approach velocity = 10 km/s | 10.0 code units |

Note: `p` replaces decimal point (e.g., `b1p5` = b=1.5 pc)

### Config Files: `{thermal}_{resolution}_{method}.json`

| Component | Options | Description |
|-----------|---------|-------------|
| thermal | `adiabatic`, `radiative` | Thermal physics model |
| resolution | `61k` | Particle count (~64,000) |
| method | `gsph`, `gdisph` | SPH formulation |

## Scenario Physics

### Impact Parameter Comparison

| Scenario | b (pc) | b/r_tidal | Disruption Strength |
|----------|--------|-----------|---------------------|
| `b1p5` | 1.5 | 0.41 | Strong (deep penetration) |
| `b3` | 3.0 | 0.83 | Moderate |
| `b6` | 6.0 | 1.65 | Weak (grazing encounter) |

Tidal radius: r_tidal ≈ 3.63 pc (where tidal force = self-gravity)

### Thermal Physics Options

| Type | Cooling | Description |
|------|---------|-------------|
| `adiabatic` | None | Pure γ=5/3 gas dynamics |
| `radiative` | Koyama-Inutsuka | ISM-like radiative cooling |

## SPH Methods

### GSPH (Godunov SPH)
- Pure Riemann solver for flux computation
- No artificial viscosity parameters
- Best for: strong shocks, blast waves

### GDISPH (Godunov DISPH)
- Riemann solver + density-independent formulation
- `avAlpha=0.0` (pure Riemann dissipation)
- Balsara switch for shock detection
- Best for: contact discontinuities, tidal streams

## Usage

### Quick Start

```bash
cd /Users/guo/Downloads/sphcode

# Copy desired preset to config.json
cp simulations/astrophysics/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph.json config.json

# Run simulation
./build/sph
```

### Available Configurations

**b=1.5 pc (Deep Penetration)**
```bash
cp .../Mc1e3_Mbh1e5_b1p5_v10/adiabatic_61k_gsph.json config.json
cp .../Mc1e3_Mbh1e5_b1p5_v10/adiabatic_61k_gdisph.json config.json
cp .../Mc1e3_Mbh1e5_b1p5_v10/radiative_61k_gsph.json config.json
cp .../Mc1e3_Mbh1e5_b1p5_v10/radiative_61k_gdisph.json config.json
```

**b=3.0 pc (Moderate)**
```bash
cp .../Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph.json config.json
cp .../Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gdisph.json config.json
cp .../Mc1e3_Mbh1e5_b3_v10/radiative_61k_gsph.json config.json
cp .../Mc1e3_Mbh1e5_b3_v10/radiative_61k_gdisph.json config.json
```

**b=6.0 pc (Grazing)**
```bash
cp .../Mc1e3_Mbh1e5_b6_v10/adiabatic_61k_gsph.json config.json
cp .../Mc1e3_Mbh1e5_b6_v10/adiabatic_61k_gdisph.json config.json
cp .../Mc1e3_Mbh1e5_b6_v10/radiative_61k_gsph.json config.json
cp .../Mc1e3_Mbh1e5_b6_v10/radiative_61k_gdisph.json config.json
```

## Output Structure

Results organized by scenario, thermal model, and method:

```
results/
├── Mc1e3_Mbh1e5_b1p5_v10/
│   ├── adiabatic/
│   │   ├── GSPH/
│   │   │   ├── snapshot_0000.csv
│   │   │   └── ...
│   │   └── GDISPH/
│   └── radiative/
│       ├── GSPH/
│       └── GDISPH/
├── Mc1e3_Mbh1e5_b3_v10/
│   └── ...
└── Mc1e3_Mbh1e5_b6_v10/
    └── ...
```

## Physical Parameters

All scenarios share these base parameters:

| Parameter | Value | Code Units | Physical |
|-----------|-------|------------|----------|
| Cloud mass | M_c | 1.0 | 10³ M☉ |
| BH mass | M_BH | 100.0 | 10⁵ M☉ |
| Approach velocity | v∞ | 10.0 | 10 km/s |
| Cloud radius | R_c | ~0.9 | ~0.9 pc |
| Tidal radius | r_tidal | ~3.63 | ~3.63 pc |
| Particle count | N | 64,000 | - |
| Simulation time | t_end | 4.0 | ~3.9 Myr |

### Unit System (IMBH_ENCOUNTER)

| Dimension | Code Unit | Physical Value |
|-----------|-----------|----------------|
| Length | 1 | 1 pc |
| Mass | 1 | 1000 M☉ |
| Velocity | 1 | 1 km/s |
| Time | 1 | 0.978 Myr |
| G | 1 | (derived) |

## Initial Condition

All configs use the same relaxed hydrostatic equilibrium IC:

```
simulations/astrophysics/imbh_cloud/results/hydrostatic/61k/GSPH_kernel_gravity/snapshot_0000.csv
```

Properties:
- Lane-Emden n=3/2 polytrope (γ=5/3)
- 64,000 particles
- Kernel-convolved self-gravity (Wendland C4)
- Verified hydrostatic equilibrium

## Key Config Parameters

```json
{
  "resumeFromSnapshot": "simulations/astrophysics/imbh_cloud/results/hydrostatic/61k/GSPH_kernel_gravity/snapshot_0000.csv",
  "gravitySofteningType": "wendland_c4",
  "theta": 0.5,
  "externalForces": {
    "pointMass": {
      "mass": 100.0,
      "position": [-20.0, <b>, 0.0],
      "velocity": [10.0, 0.0, 0.0]
    }
  },
  "thermal": {
    "cooling": "none" | "koyama_inutsuka"
  }
}
```

## Creating New Scenarios

1. Create directory: `Mc{X}_Mbh{Y}_b{Z}_v{W}/`
2. Copy existing config as template
3. Update `externalForces.pointMass.position[1]` for impact parameter
4. Update `thermal.cooling` as needed
5. Update `outputDirectory` and `checkpointPath`
6. Update `SIMULATION_MATRIX.md` with new scenario

## See Also

- `SIMULATION_MATRIX.md` - Complete physics derivations and parameter table
- `../../README.md` - Full workflow documentation
