# IMBH-Cloud Encounter Simulation Parameters

## Overview

This document provides the complete physics-based parameter matrix for IMBH (Intermediate Mass Black Hole) - molecular cloud encounter simulations. All parameters are derived from first principles using a γ=5/3 Lane-Emden polytropic cloud model.

---

## Unit System

The simulation uses the **IMBH_ENCOUNTER** unit system optimized for galactic-scale cloud-BH interactions:

| Quantity | Code Unit | Physical Value | Symbol |
|----------|-----------|----------------|--------|
| Length | 1 | 1 pc = 3.086×10¹⁸ cm | L |
| Mass | 1 | 1000 M☉ = 1.988×10³⁶ g | M |
| Velocity | 1 | 1 km/s = 10⁵ cm/s | V |
| Time | 1 | L/V = 977.8 kyr ≈ 1 Myr | T |
| Density | 1 | M/L³ = 6.77×10⁻²⁰ g/cm³ | ρ |
| G | 1 | (Normalized) | G |

**Note:** G is set to 1 in code units, which provides natural scaling for self-gravitating systems.

---

## Initial Condition: Lane-Emden γ=5/3 Polytrope

The molecular cloud is modeled as a self-gravitating polytropic sphere in hydrostatic equilibrium.

### Lane-Emden Solution (n = 3/2)

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Polytropic index | n | 3/2 |
| Adiabatic index | γ | 5/3 = 1.6667 |
| First zero | ξ₁ | 3.65375 |
| Derivative at ξ₁ | θ'(ξ₁) | -0.20330 |
| Dimensionless mass | ωₙ | 2.714 |

### Cloud Physical Properties

| Property | Code Value | Physical Value | Notes |
|----------|------------|----------------|-------|
| **Total mass** | M_cloud = 1.0 | 1000 M☉ | Normalized |
| **Radius** | R_cloud ≈ 1.13 | 1.13 pc | Surface where ρ→0 |
| **Central density** | ρ_c = 1.0 | 6.77×10⁻²⁰ g/cm³ | Lane-Emden normalized |
| **Polytropic constant** | K = 1.0 | (code units) | P = K ρ^γ |

### Particle Resolution

| Resolution | N_particles | Particle Mass (code) | Particle Mass (M☉) |
|------------|-------------|----------------------|---------------------|
| 10k | 10,000 | 1.0×10⁻⁴ | 0.1 |
| **61k** | 64,000 | 1.5625×10⁻⁵ | 0.0156 |
| 200k | 200,000 | 5.0×10⁻⁶ | 0.005 |

---

## IMBH Properties

| Property | Code Value | Physical Value |
|----------|------------|----------------|
| **BH mass** | M_BH = 100 | 10⁵ M☉ |
| **Position** | (0, 0, 0) | Origin |
| **Softening** | ε = 0.05 | 0.05 pc |
| **Mass ratio** | M_BH/M_cloud | 100 |

### Tidal Radius

The tidal radius where BH gravity equals cloud self-gravity:

```
r_tidal = R_cloud × (M_BH / 3M_cloud)^(1/3) = 3.63 pc
```

---

## Encounter Geometry

### Initial Setup

| Parameter | Value | Physical Meaning |
|-----------|-------|------------------|
| Cloud start position | (-20, b, 0) pc | 20 pc from BH along -x |
| Cloud velocity | (10, 0, 0) km/s | Moving toward BH |
| Impact parameter | b | Perpendicular offset |

### Impact Parameter Scenarios

| Scenario | b (pc) | b/r_tidal | b/R_cloud | Tidal Regime |
|----------|--------|-----------|-----------|--------------|
| **b1p5pc** | 1.5 | 0.41 | 1.33 | **Strong disruption** |
| **b3pc** | 3.0 | 0.83 | 2.66 | **Partial disruption** |
| **b6pc** | 6.0 | 1.65 | 5.33 | **Tidal stripping** |

### Regime Definitions

- **Strong disruption** (b/r_t < 0.5): Cloud passes deep within tidal radius, severe mass loss
- **Partial disruption** (0.5 < b/r_t < 1.0): Significant tidal stripping, cloud survives partially
- **Tidal stripping** (1.0 < b/r_t < 2.0): Outer layers stripped, core survives
- **Weak interaction** (b/r_t > 2.0): Minimal disruption, small perturbation

---

## Time Scales

| Time Scale | Code Value | Physical Value | Definition |
|------------|------------|----------------|------------|
| Free-fall time | t_ff = 0.54 | 531 kyr | √(3π/32Gρ_c) |
| Sound crossing | t_cross = 0.87 | 853 kyr | R_cloud/c_s |
| Approach time | t_approach = 2.0 | 1.96 Myr | |x₀|/v∞ |
| **Simulation end** | t_end = 4.0 | 3.9 Myr | Full encounter + aftermath |

---

## Thermal Physics Options

### Adiabatic (No Cooling)

- Pure pressure-gravity dynamics
- γ = 5/3 polytropic EOS maintained
- Energy equation: dε/dt = -P∇·v (compression/expansion only)
- Config: `"thermal": {"cooling": false, "heating": false}`

### Radiative (Koyama-Inutsuka Cooling)

- Realistic ISM thermal balance
- Net cooling function: Λ(T) - Γ
- Thermal instability captured
- Config: `"thermal": {"cooling": true, "heating": true, "coolingType": "koyama_inutsuka"}`

---

## Simulation Parameter Matrix

### Full Matrix (Impact × Thermal × Resolution × Method)

| Impact | Thermal | Resolution | SPH Method | Scenario ID |
|--------|---------|------------|------------|-------------|
| b1.5pc | adiabatic | 61k | GSPH | `Mc1e3_Mbh1e5_b1p5_v10_adiabatic_61k_gsph` |
| b1.5pc | adiabatic | 61k | GDISPH | `Mc1e3_Mbh1e5_b1p5_v10_adiabatic_61k_gdisph` |
| b1.5pc | radiative | 61k | GSPH | `Mc1e3_Mbh1e5_b1p5_v10_radiative_61k_gsph` |
| b1.5pc | radiative | 61k | GDISPH | `Mc1e3_Mbh1e5_b1p5_v10_radiative_61k_gdisph` |
| b3pc | adiabatic | 61k | GSPH | `Mc1e3_Mbh1e5_b3_v10_adiabatic_61k_gsph` |
| b3pc | adiabatic | 61k | GDISPH | `Mc1e3_Mbh1e5_b3_v10_adiabatic_61k_gdisph` |
| b3pc | radiative | 61k | GSPH | `Mc1e3_Mbh1e5_b3_v10_radiative_61k_gsph` |
| b3pc | radiative | 61k | GDISPH | `Mc1e3_Mbh1e5_b3_v10_radiative_61k_gdisph` |
| b6pc | adiabatic | 61k | GSPH | `Mc1e3_Mbh1e5_b6_v10_adiabatic_61k_gsph` |
| b6pc | adiabatic | 61k | GDISPH | `Mc1e3_Mbh1e5_b6_v10_adiabatic_61k_gdisph` |
| b6pc | radiative | 61k | GSPH | `Mc1e3_Mbh1e5_b6_v10_radiative_61k_gsph` |
| b6pc | radiative | 61k | GDISPH | `Mc1e3_Mbh1e5_b6_v10_radiative_61k_gdisph` |

### Naming Convention

```
Mc{cloud_mass}_Mbh{bh_mass}_b{impact_param}_v{velocity}_{thermal}_{resolution}_{method}
```

Where:
- `Mc1e3` = Cloud mass 10³ M☉
- `Mbh1e5` = BH mass 10⁵ M☉  
- `b{X}` = Impact parameter X pc (e.g., b1p5 = 1.5 pc, b3 = 3 pc)
- `v10` = Approach velocity 10 km/s
- `adiabatic` or `radiative` = Thermal physics
- `61k` = Particle resolution
- `gsph` or `gdisph` = SPH method

---

## Directory Structure (Recommended)

```
simulations/astrophysics/imbh_cloud/
├── config/presets/
│   └── simulation/
│       └── scenarios/
│           ├── Mc1e3_Mbh1e5_b1p5_v10/
│           │   ├── adiabatic_61k_gsph.json
│           │   ├── adiabatic_61k_gdisph.json
│           │   ├── radiative_61k_gsph.json
│           │   └── radiative_61k_gdisph.json
│           ├── Mc1e3_Mbh1e5_b3_v10/
│           │   └── ...
│           ├── Mc1e3_Mbh1e5_b6_v10/
│           │   └── ...
│           └── ic_relaxation/
│               ├── 10k/
│               ├── 61k/
│               └── 200k/
├── results/
│   ├── Mc1e3_Mbh1e5_b1p5_v10/
│   │   ├── adiabatic_61k_gsph/
│   │   ├── adiabatic_61k_gdisph/
│   │   └── ...
│   └── ...
└── docs/
    └── SIMULATION_PARAMETERS.md (this file)
```

---

## SPH Method Comparison

### GSPH (Godunov SPH)

- Riemann solver: HLL
- Artificial viscosity: None (dissipation via Riemann solver)
- Best for: Strong shocks, clean discontinuities

### GDISPH (Godunov Density-Independent SPH)

- Riemann solver: HLL
- Artificial viscosity: α=0 (pure Riemann)
- Balsara switch: Enabled
- Best for: Contact discontinuities, mixing, tidal streams

---

## Gravity Settings

All encounter simulations use **kernel-convolved self-gravity** for consistency:

```json
{
  "useGravity": true,
  "G": 1.0,
  "theta": 0.5,
  "gravitySofteningType": "wendland_c4"
}
```

The `wendland_c4` kernel softening matches the SPH kernel, providing consistent force resolution.

---

## Initial Condition (IC) Files

| Resolution | IC File Path | Relaxation Steps |
|------------|--------------|------------------|
| 10k | `results/ic_relaxation/10k/gsph/snapshot_final.csv` | 50,000 |
| 61k | `results/ic_relaxation/61k/gsph/snapshot_final.csv` | 200,000 |
| 200k | `results/ic_relaxation/200k/gsph/snapshot_final.csv` | 500,000 |

---

## Quick Reference: Key Numbers

| Parameter | Value |
|-----------|-------|
| Cloud mass | 1000 M☉ |
| Cloud radius | 1.13 pc |
| BH mass | 10⁵ M☉ |
| Tidal radius | 3.63 pc |
| Approach velocity | 10 km/s |
| Simulation duration | ~4 Myr |
| 1 code time | ~1 Myr |

---

## References

1. Lane-Emden equation: Chandrasekhar, S. (1939), "An Introduction to the Study of Stellar Structure"
2. Tidal disruption physics: Hills, J.G. (1975), Nature 254, 295
3. Koyama-Inutsuka cooling: Koyama, H. & Inutsuka, S. (2002), ApJ 564, L97

---

*Last updated: 2025-12-10*
