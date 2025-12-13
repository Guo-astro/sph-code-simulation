# IMBH-Cloud Simulation Category System

## Overview

This document defines the **systematic category naming** for IMBH-cloud tidal disruption simulations. The system is designed for easy-to-use Makefile targets and clear physics traceability.

## Category Naming Convention

```
CAT{physics}_{resolution}_{scenario}
```

Example: `CAT1_B_rp05` = Category 1 (Adiabatic), Resolution B (200k), Pericenter r_p=0.5 pc

---

## Physics Categories

| Category | Code | Physics | Key Features |
|----------|------|---------|--------------|
| **CAT1** | `adiabatic` | Adiabatic (γ=5/3) | Pure gas dynamics, no cooling/heating |
| **CAT3** | `radiative_sg` | Radiative + Self-gravity | Full physics, potential fragmentation |
| **CAT_OKA** | `oka_exact` | Oka et al. (2017) parameters | Hyperbolic orbit matching observations |

### CAT1: Adiabatic Dynamics
- **Physics**: Pure hydrodynamics with γ=5/3 polytropic EOS
- **Purpose**: Baseline for understanding pure tidal dynamics
- **Key Observable**: Tidal stream morphology, shock heating
- **Runtime**: Fastest (no cooling subcycling)
- **Self-gravity**: Cloud self-gravity enabled

### CAT3: Full Physics (Radiative + Self-Gravity)
- **Physics**: Inoue & Inutsuka (2008) ISM cooling + cloud self-gravity
- **Cooling function**: Λ/Γ = 10^7 exp(-114800/(T+1000)) + 0.014 √T exp(-92/T)
- **Purpose**: Most physically complete simulation
- **Key Observable**: Fragmentation, bound clump formation, density enhancement
- **Runtime**: Slowest (gravity + cooling)
- **Note**: Requires sufficient resolution to resolve Jeans length

### CAT_OKA: Oka et al. (2017) Exact Parameters
- **Physics**: Radiative cooling + self-gravity with Oka hyperbolic orbit
- **Orbital parameters**: Hyperbolic trajectory with pericentre = 1.69 pc, e = 1.24
- **Purpose**: Direct comparison with CO-0.40-0.22 observations
- **Key Observable**: P-V diagram morphology, velocity width ~100 km/s

---

## Resolution Levels

| Level | Code | Particles | Smoothing Length | Use Case |
|-------|------|-----------|------------------|----------|
| **A** | `61k` | 61,000 | ~0.03 pc | Quick tests, parameter scans |
| **B** | `200k` | 200,000 | ~0.02 pc | Production, convergence study |
| **C** | `1M` | 1,000,000 | ~0.01 pc | High-resolution, publication |

### Resolution Guidelines

1. **Jeans Resolution** (for CAT3):
   ```
   h_max < λ_J / 4 = 0.004 pc (for n=10^4 cm^-3, T=60K)
   ```
   - Level A: Marginal (may miss fragmentation)
   - Level B: Adequate for most purposes
   - Level C: Well-resolved, suitable for convergence

2. **Pericenter Resolution** (for close encounters):
   - Level A: ~3-5 particles across r_peri = 0.3 pc
   - Level B: ~10 particles across r_peri = 0.3 pc
   - Level C: ~30 particles across r_peri = 0.3 pc

---

## Scenario Parameters: Impact Parameter / Pericenter Study

### Close Encounter Scenarios (Pericenter 0.3-1.0 pc)

These scenarios probe the **deep disruption regime** where pericenter < cloud radius (R_cloud = 1.13 pc).

| Scenario | r_peri (pc) | r_peri/R_cloud | b (pc) | Eccentricity | Regime |
|----------|-------------|----------------|--------|--------------|--------|
| `rp03` | 0.3 | 0.27 | 2.08 | 1.043 | **Extreme disruption** |
| `rp05` | 0.5 | 0.44 | 2.70 | 1.071 | **Very close encounter** |
| `rp07` | 0.7 | 0.62 | 3.21 | 1.100 | **Close encounter** |
| `rp10` | 1.0 | 0.88 | 3.88 | 1.142 | **Moderate close** |

### Standard Scenarios (Pericenter > R_cloud)

| Scenario | r_peri (pc) | r_peri/R_cloud | b (pc) | Regime |
|----------|-------------|----------------|--------|--------|
| `rp15` | 1.5 | 1.33 | 4.83 | Grazing disruption |
| `rp20` | 2.0 | 1.77 | 5.64 | Partial disruption |
| `rp30` | 3.0 | 2.65 | 7.02 | Weak interaction |

### Orbital Parameters

All hyperbolic scenarios use:
- **Cloud mass**: M_cloud = 1000 M☉
- **BH mass**: M_BH = 10^5 M☉
- **Velocity at infinity**: v_∞ = 8 km/s
- **Semi-major axis**: a = -7.03 pc (hyperbolic)
- **Cloud radius**: R_cloud = 1.13 pc
- **Tidal radius**: r_tidal = 3.64 pc (at r = R_cloud)

### Orbital Mechanics Summary

For hyperbolic orbit with v_∞ = 8 km/s:
```
GM_BH = 449.8 pc·(km/s)²
ε = v_∞²/2 = 32 (km/s)²  [specific orbital energy]
a = -GM/(2ε) = -7.03 pc  [semi-major axis, negative for hyperbola]
e = 1 + r_peri/|a|       [eccentricity]
h = √(GM·r_peri·(1+e))   [specific angular momentum]
b = h/v_∞                [impact parameter]
```

### Initial Conditions for Close Encounters

Starting position: r_0 = 20 pc from BH

| Scenario | Position (pc) | Velocity (km/s) | v_peri (km/s) |
|----------|--------------|-----------------|---------------|
| `rp03` | (19.89, -2.08, 0) | (-7.92, 0.83, 0) | 54.7 |
| `rp05` | (19.82, -2.70, 0) | (-7.85, 1.08, 0) | 42.4 |
| `rp07` | (19.74, -3.21, 0) | (-7.78, 1.28, 0) | 35.8 |
| `rp10` | (19.62, -3.88, 0) | (-7.68, 1.55, 0) | 30.0 |
| `rp15` | (19.41, -4.83, 0) | (-7.52, 1.93, 0) | 24.5 |
| `rp20` | (19.19, -5.64, 0) | (-7.35, 2.25, 0) | 21.2 |
| `rp30` | (18.73, -7.02, 0) | (-7.02, 2.81, 0) | 17.3 |

---

## Complete Simulation Matrix

### Priority 1: Close Encounter Study (New)

| Config Name | Category | Resolution | Pericenter | Particles |
|-------------|----------|------------|------------|-----------|
| `CAT1_B_rp03` | CAT1 | B (200k) | 0.3 pc | 200,000 |
| `CAT1_B_rp05` | CAT1 | B (200k) | 0.5 pc | 200,000 |
| `CAT1_B_rp07` | CAT1 | B (200k) | 0.7 pc | 200,000 |
| `CAT1_B_rp10` | CAT1 | B (200k) | 1.0 pc | 200,000 |
| `CAT3_B_rp03` | CAT3 | B (200k) | 0.3 pc | 200,000 |
| `CAT3_B_rp05` | CAT3 | B (200k) | 0.5 pc | 200,000 |
| `CAT3_B_rp07` | CAT3 | B (200k) | 0.7 pc | 200,000 |
| `CAT3_B_rp10` | CAT3 | B (200k) | 1.0 pc | 200,000 |

### Priority 2: High-Resolution Close Encounters

| Config Name | Category | Resolution | Pericenter | Particles |
|-------------|----------|------------|------------|-----------|
| `CAT3_C_rp03` | CAT3 | C (1M) | 0.3 pc | 1,000,000 |
| `CAT3_C_rp05` | CAT3 | C (1M) | 0.5 pc | 1,000,000 |

### Legacy Scenarios (Retained for Comparison)

| Config Name | Category | Resolution | Impact b (pc) | Notes |
|-------------|----------|------------|---------------|-------|
| `CAT1_A_b1p5` | CAT1 | A | 1.5 | Legacy format |
| `CAT1_B_b3` | CAT1 | B | 3.0 | Legacy format |
| `CAT_OKA_B_oka` | CAT_OKA | B | 5.17 | Oka exact |

---

## Physical Regimes

### Deep Disruption (r_peri < R_cloud)

When pericenter is smaller than cloud radius, the cloud **passes through** the BH position:

```
        r_peri = 0.3 pc                    r_peri = 1.0 pc
        ════════════════                   ════════════════

              ┌─────┐                           ┌─────┐
            ╱ │     │ ╲                       ╱ │     │ ╲
           │  │ BH  │  │                     │  │     │  │
           │  │  ●──┼──│── r_peri           │  ●─────┼──│── r_peri
           │  │     │  │                     │  BH   │  │
            ╲ │     │ ╱                       ╲ │     │ ╱
              └─────┘                           └─────┘
           Cloud envelope                    Cloud envelope
           passes through BH                 grazes BH position

           Extreme tidal                     Strong tidal
           disruption                        disruption
```

**Expected Outcomes:**
- **r_peri = 0.3 pc**: Complete disruption, cloud torn into multiple streams
- **r_peri = 0.5 pc**: Severe disruption, extended tidal tails
- **r_peri = 0.7 pc**: Strong disruption, possible core survival
- **r_peri = 1.0 pc**: Significant disruption, core likely survives

### Tidal Strength Parameter

The key physics parameter is the **tidal strength**:

$$\frac{r_\mathrm{peri}}{r_\mathrm{tidal}} = \frac{r_\mathrm{peri}}{R_\mathrm{cloud}} \left(\frac{M_\mathrm{cloud}}{3 M_\mathrm{BH}}\right)^{-1/3}$$

| Scenario | r_peri/R_cloud | r_peri/r_tidal | Tidal Regime |
|----------|---------------|----------------|--------------|
| rp03 | 0.27 | 0.08 | **Extreme** (>> 1 crossing) |
| rp05 | 0.44 | 0.14 | **Very strong** |
| rp07 | 0.62 | 0.19 | **Strong** |
| rp10 | 0.88 | 0.27 | **Moderate-strong** |
| rp15 | 1.33 | 0.41 | Moderate |
| rp20 | 1.77 | 0.55 | Weak-moderate |

---

## Makefile Usage

### One-Command Execution

```bash
# Run specific configuration
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud CAT1_B_rp05
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud CAT3_B_rp03

# Run close encounter parameter study
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud CAT1_B_close_encounters
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud CAT3_B_close_encounters

# Compare pericenter variation
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud compare_rperi_CAT1
```

### Individual Targets

```bash
# Simulation only (no animation)
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud CAT3_B_rp05_run

# Animation only (after simulation)
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud CAT3_B_rp05_anim

# Clean specific run
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud CAT3_B_rp05_clean
```

---

## Output Directory Structure

```
simulations/astrophysics/imbh_cloud/
├── results/
│   ├── CAT1/                          # Adiabatic simulations
│   │   ├── A_61k/
│   │   │   ├── rp03/                  # Pericenter 0.3 pc
│   │   │   ├── rp05/                  # Pericenter 0.5 pc
│   │   │   ├── rp07/                  # Pericenter 0.7 pc
│   │   │   ├── rp10/                  # Pericenter 1.0 pc
│   │   │   └── b1p5/                  # Legacy: impact b=1.5 pc
│   │   ├── B_200k/
│   │   └── C_1M/
│   ├── CAT3/                          # Full physics simulations
│   │   ├── A_61k/
│   │   ├── B_200k/
│   │   └── C_1M/
│   └── CAT_OKA/                       # Oka et al. exact
│       ├── A_61k/
│       ├── B_200k/
│       └── C_1M/
└── animations/
    └── (same structure as results)
```

---

## Recommended Workflow

### 1. Close Encounter Parameter Scan (Priority)

```bash
# Quick scan at low resolution
for rp in rp03 rp05 rp07 rp10; do
    make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud CAT1_A_${rp}
done

# Production runs at medium resolution
for rp in rp03 rp05 rp07 rp10; do
    make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud CAT3_B_${rp}
done
```

### 2. Physics Comparison (Fixed Pericenter)

```bash
# Compare CAT1 vs CAT3 at r_peri = 0.5 pc
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud CAT1_B_rp05
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud CAT3_B_rp05
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud compare_CAT_rp05
```

### 3. Resolution Convergence

```bash
# Test convergence for CAT3 at r_peri = 0.5 pc
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud CAT3_A_rp05
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud CAT3_B_rp05
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud CAT3_C_rp05
```

### 4. Publication Run

```bash
# High-resolution extreme disruption
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud CAT3_C_rp03
```

---

## Initial Conditions

All simulations share the same relaxed Lane-Emden sphere:

| Resolution | IC File | Source |
|------------|---------|--------|
| A (61k) | `results/IC/61k/snapshot_0000.csv` | Relaxed polytrope |
| B (200k) | `results/IC/200k/snapshot_0000.csv` | Relaxed polytrope |
| C (1M) | `results/IC/1M/snapshot_0000.csv` | Relaxed polytrope |

---

## Runtime Estimates

| Config | Resolution | r_peri | Est. Runtime | Notes |
|--------|------------|--------|--------------|-------|
| CAT1_A_rp05 | 61k | 0.5 pc | ~15 min | Quick test |
| CAT1_B_rp05 | 200k | 0.5 pc | ~2 hours | Production |
| CAT3_B_rp05 | 200k | 0.5 pc | ~4 hours | Full physics |
| CAT3_B_rp03 | 200k | 0.3 pc | ~6 hours | More timesteps |
| CAT3_C_rp03 | 1M | 0.3 pc | ~48 hours | Publication |

**Note**: Close encounters (small r_peri) require more timesteps due to higher velocities near pericenter.

---

## Physical Questions Addressed

| Question | Recommended Config | Notes |
|----------|-------------------|-------|
| Minimum pericenter for core survival? | CAT3_B_rp03/05/07/10 | Compare outcomes |
| Stream morphology vs r_peri? | CAT1_B_rp03-rp10 | Pure kinematics |
| Does cooling affect disruption threshold? | CAT1 vs CAT3 at rp05 | Physics comparison |
| Velocity at pericenter? | All rp configs | v_peri ∝ 1/√r_peri |
| P-V diagram vs observations? | CAT_OKA_B_oka | Compare with Oka |
| Maximum compression ratio? | CAT3_B_rp03 | Smallest r_peri |

---

## Numerical Considerations for Close Encounters

### Softening Length

The BH softening ε = 0.05 pc prevents numerical singularities:
- For rp03 (r_peri = 0.3 pc): r_peri/ε = 6 (well-resolved)
- Gravitational force is Plummer-softened: F ∝ 1/(r² + ε²)

### Timestep Constraints

Near pericenter, velocities are high:
- v_peri(rp03) ≈ 55 km/s
- v_peri(rp10) ≈ 30 km/s

CFL condition requires smaller timesteps for close encounters.

### Mass Conservation Check

For extreme disruption, monitor:
- Total mass in domain
- Mass crossing softening sphere
- Numerical artifacts near BH

---

## References

1. Oka, T., et al. (2017). Nature Astronomy - CO-0.40-0.22 observations
2. Theoretical derivation: `docs/imbh-sim/theoretical_derivation_imbh_cloud_scattering.md`
