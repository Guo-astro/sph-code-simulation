# IMBH-Cloud Encounter Simulation Matrix

## Overview

This document defines the simulation parameter space for studying tidal disruption of molecular clouds by intermediate-mass black holes (IMBHs). The simulations use a 5/3-polytrope (Lane-Emden n=3/2) cloud as the initial condition.

## Unit System

The IMBH_ENCOUNTER unit system uses:
- **Length unit**: 1 pc (parsec) = 3.086×10¹⁸ cm
- **Mass unit**: 1000 M☉ (solar masses) = 1.989×10³⁶ g  
- **Velocity unit**: 1 km/s = 10⁵ cm/s
- **Time unit**: ~0.978 Myr (derived: L/V)
- **G**: 1.0 in code units (self-consistent)

## Initial Cloud Properties (from Hydrostatic Equilibrium)

The cloud is a γ=5/3 polytrope (Lane-Emden n=3/2) in hydrostatic equilibrium.

### Measured from 61k Relaxed Snapshot

| Property | Code Units | Physical Units | Notes |
|----------|------------|----------------|-------|
| Total Mass (M_cloud) | 1.0 | 1000 M☉ | Sum of all particles |
| Central Density (ρ_c) | ~1.4 | ~1.4 × ρ_unit | From snapshot max(dens) |
| Cloud Radius (R_cloud) | ~0.9 | ~0.9 pc | 90% containment |
| Central Sound Speed | ~0.94 | ~0.94 km/s | At max density |
| Central Pressure | ~0.74 | | P = K ρ^γ |
| Polytropic K | ~0.30 | | K = P_c / ρ_c^γ |
| Particle Count | 64000 | | N=40 per dimension |
| Particle Mass | 1.5625×10⁻⁵ | ~0.0156 M☉ | M_cloud / N_particles |

### Lane-Emden n=3/2 Theory

For a γ=5/3 polytrope (n=3/2):
- First zero: ξ₁ = 3.65375
- Mass parameter: ω_n = -ξ²θ'(ξ₁) = 2.71406
- R = α × ξ₁ where α² = K(n+1)ρ_c^(1/n-1)/(4πG)

## Black Hole Parameters

| Property | Code Units | Physical Units |
|----------|------------|----------------|
| IMBH Mass (M_BH) | 100 | 10⁵ M☉ |
| Position | (0, 0, 0) | Origin |
| Softening | 0.05 | 0.05 pc |

### Tidal Radius

The tidal radius where BH tidal force equals cloud self-gravity:

```
r_tidal = R_cloud × (M_BH / M_cloud)^(1/3)
        ≈ 0.9 × (100/1)^(1/3)
        ≈ 0.9 × 4.64
        ≈ 4.18 pc (in code units)
        ≈ 3.63 pc (using exact R_cloud from snapshot)
```

## Simulation Matrix

### Parameter Space

| Parameter | Values | Description |
|-----------|--------|-------------|
| Impact parameter (b) | 1.5, 3.0, 6.0 pc | Cloud center y-offset |
| Approach velocity (v) | 10 km/s | Fixed |
| Thermal physics | Adiabatic, Radiative | Cooling on/off |
| SPH Method | GSPH, GDISPH | Godunov methods |
| Resolution | 61k (64000 particles) | Standard |

### Scenario Matrix

| Folder Name | b (pc) | b/r_tidal | Regime | Thermal |
|-------------|--------|-----------|--------|---------|
| `Mc1e3_Mbh1e5_b1p5_v10/adiabatic_*` | 1.5 | 0.41 | Deep disruption | Adiabatic |
| `Mc1e3_Mbh1e5_b1p5_v10/radiative_*` | 1.5 | 0.41 | Deep disruption | Radiative |
| `Mc1e3_Mbh1e5_b3_v10/adiabatic_*` | 3.0 | 0.83 | Partial disruption | Adiabatic |
| `Mc1e3_Mbh1e5_b3_v10/radiative_*` | 3.0 | 0.83 | Partial disruption | Radiative |
| `Mc1e3_Mbh1e5_b6_v10/adiabatic_*` | 6.0 | 1.65 | Weak interaction | Adiabatic |
| `Mc1e3_Mbh1e5_b6_v10/radiative_*` | 6.0 | 1.65 | Weak interaction | Radiative |

### Naming Convention

```
Mc{cloud_mass}_Mbh{bh_mass}_b{impact_param}_v{velocity}/
    {thermal}_{resolution}_{sph_method}.json
```

Example: `Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph.json`
- Cloud mass: 10³ M☉
- BH mass: 10⁵ M☉
- Impact parameter: 3 pc
- Velocity: 10 km/s
- Thermal: adiabatic (no cooling)
- Resolution: 61k particles
- SPH method: GSPH

## Initial Condition Setup

The cloud is placed at:
- **Position**: (-20, b, 0) pc, where b = impact parameter
- **Velocity**: (10, 0, 0) km/s (toward BH at origin)

The cloud approaches from negative x, passes closest to BH at x≈0, then continues along +x.

## Thermal Physics Options

### Adiabatic
- `cooling: false, heating: false`
- Pure γ=5/3 polytropic evolution
- Energy equation: dε/dt = -(P/ρ)∇·v

### Radiative (Koyama-Inutsuka Cooling)
- `cooling: true, heating: true, coolingType: "koyama_inutsuka"`
- Radiative cooling/heating appropriate for ISM conditions
- Allows thermal instability and fragmentation

## Gravity Treatment

All simulations use:
- **Kernel-convolved self-gravity**: `gravitySofteningType: "wendland_c4"`
- **Tree opening angle**: `theta: 0.5`
- **External point mass BH**: Plummer softening = 0.05 pc

## Timescales

| Timescale | Expression | Value (code units) | Physical |
|-----------|------------|-------------------|----------|
| Free-fall | √(3π/32Gρ_c) | ~0.46 | ~0.45 Myr |
| Sound crossing | R_cloud/c_s | ~0.96 | ~0.94 Myr |
| Encounter | 2×20/10 | 4.0 | ~3.9 Myr |
| Dynamical | √(R³/GM) | ~0.85 | ~0.83 Myr |

## Directory Structure

```
simulations/astrophysics/imbh_cloud/
├── config/presets/simulation/scenarios/
│   ├── Mc1e3_Mbh1e5_b1p5_v10/       # b=1.5 pc scenarios
│   │   ├── adiabatic_61k_gsph.json
│   │   ├── adiabatic_61k_gdisph.json
│   │   ├── radiative_61k_gsph.json
│   │   └── radiative_61k_gdisph.json
│   ├── Mc1e3_Mbh1e5_b3_v10/         # b=3.0 pc scenarios
│   │   └── ...
│   ├── Mc1e3_Mbh1e5_b6_v10/         # b=6.0 pc scenarios
│   │   └── ...
│   └── ic_relaxation/               # Initial condition generation
│       └── 61k/
├── results/
│   ├── hydrostatic/61k/GSPH_kernel_gravity/  # Relaxed IC snapshots
│   └── Mc1e3_Mbh1e5_*/               # Simulation outputs
└── scripts/                          # Analysis tools
```

## Running Simulations

### Generate Initial Condition (if needed)

```bash
make -f simulations/astrophysics/imbh_cloud/Makefile.hydrostatic hydro_oneshot SIZE=61k METHOD=gsph GRAVITY_TYPE=kernel
```

### Run Encounter Simulation

```bash
# Example: b=3pc, adiabatic, GSPH
./build/sph simulations/astrophysics/imbh_cloud/config/presets/simulation/scenarios/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph.json
```

## Expected Physics

### b=1.5 pc (Deep Disruption, b/r_tidal ≈ 0.41)
- Strong tidal forces dominate
- Complete or near-complete disruption
- Long tidal streams form
- Significant mass loss from cloud

### b=3.0 pc (Partial Disruption, b/r_tidal ≈ 0.83)
- Tidal forces comparable to self-gravity
- Partial stripping of outer layers
- Core may survive
- Tidal tails develop

### b=6.0 pc (Weak Interaction, b/r_tidal ≈ 1.65)
- Cloud mostly outside tidal radius
- Weak tidal distortion
- Possible oscillations/breathing modes
- Cloud survives largely intact

## Notes on Radiative Cooling

When Koyama-Inutsuka cooling is enabled:
- Cooling can trigger thermal instability
- Dense regions cool faster → fragmentation
- May form clumps/filaments in tidal streams
- Significantly different evolution than adiabatic case
- Important for star formation in tidal debris
