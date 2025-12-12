# IMBH-Cloud Simulation Category System

## Overview

This document defines the **systematic category naming** for IMBH-cloud tidal disruption simulations. The system is designed for easy-to-use Makefile targets and clear physics traceability.

## Category Naming Convention

```
CAT{physics}_{resolution}_{scenario}
```

Example: `CAT2_B_b3` = Category 2 (radiative), Resolution B (200k), Impact parameter b=3pc

---

## Physics Categories

| Category | Code | Physics | Key Features |
|----------|------|---------|--------------|
| **CAT1** | `adiabatic` | Adiabatic (γ=5/3) | Pure gas dynamics, no cooling/heating |
| **CAT2** | `radiative` | Inoue-Inutsuka cooling | ISM thermal equilibrium, density enhancement |
| **CAT3** | `radiative_sg` | Radiative + Self-gravity | Full physics, potential fragmentation |

### CAT1: Adiabatic Dynamics
- **Physics**: Pure hydrodynamics with γ=5/3 polytropic EOS
- **Purpose**: Baseline for understanding pure tidal dynamics
- **Key Observable**: Tidal stream morphology, shock heating
- **Runtime**: Fastest (no cooling subcycling)

### CAT2: Radiative Cooling/Heating  
- **Physics**: Inoue & Inutsuka (2008) ISM cooling function
- **Cooling function**: Λ/Γ = 10^7 exp(-114800/(T+1000)) + 0.014 √T exp(-92/T)
- **Purpose**: Realistic ISM thermal physics
- **Key Observable**: Density enhancement, thermal instability
- **Runtime**: Moderate (cooling timestep constraint)

### CAT3: Full Physics (Radiative + Self-Gravity)
- **Physics**: Cooling + cloud self-gravity enabled
- **Purpose**: Most physically complete simulation
- **Key Observable**: Fragmentation, bound clump formation
- **Runtime**: Slowest (gravity + cooling)
- **Note**: Requires sufficient resolution to resolve Jeans length

---

## Resolution Levels

| Level | Code | Particles | Smoothing Length | Use Case |
|-------|------|-----------|------------------|----------|
| **A** | `60k` | 64,000 | ~0.03 pc | Quick tests, parameter scans |
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

2. **Tidal Stream Resolution**:
   - Level A: ~10 particles across stream cross-section
   - Level B: ~20 particles across stream
   - Level C: ~50 particles across stream

---

## Scenario Parameters

### Impact Parameter Variations

| Scenario | b (pc) | b/r_tidal | Regime |
|----------|--------|-----------|--------|
| `b1p5` | 1.5 | 0.41 | Deep disruption |
| `b3` | 3.0 | 0.83 | Partial disruption |
| `b6` | 6.0 | 1.65 | Weak interaction |

All scenarios use:
- Cloud mass: M_cloud = 1000 M☉
- BH mass: M_BH = 10^5 M☉
- Approach velocity: v = 10 km/s

---

## Complete Simulation Matrix

### Full 27-Configuration Matrix (per impact parameter)

| Config Name | Category | Resolution | Thermal | Self-Gravity | Particles |
|-------------|----------|------------|---------|--------------|-----------|
| `CAT1_A_b3` | CAT1 | A (60k) | Adiabatic | Cloud only | 64,000 |
| `CAT1_B_b3` | CAT1 | B (200k) | Adiabatic | Cloud only | 200,000 |
| `CAT1_C_b3` | CAT1 | C (1M) | Adiabatic | Cloud only | 1,000,000 |
| `CAT2_A_b3` | CAT2 | A (60k) | Radiative | Cloud only | 64,000 |
| `CAT2_B_b3` | CAT2 | B (200k) | Radiative | Cloud only | 200,000 |
| `CAT2_C_b3` | CAT2 | C (1M) | Radiative | Cloud only | 1,000,000 |
| `CAT3_A_b3` | CAT3 | A (60k) | Radiative | Yes | 64,000 |
| `CAT3_B_b3` | CAT3 | B (200k) | Radiative | Yes | 200,000 |
| `CAT3_C_b3` | CAT3 | C (1M) | Radiative | Yes | 1,000,000 |

For b=1.5pc and b=6pc: replace `_b3` with `_b1p5` or `_b6`.

---

## Makefile Usage

### One-Command Execution

```bash
# Run specific configuration
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT1_A_b3
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT2_B_b3
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT3_C_b3

# Run all resolutions for a category/scenario
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT1_all_b3

# Run all categories for a resolution/scenario
make -f sample/imbh_cloud/Makefile.imbh_cloud all_A_b3

# Run entire category across all scenarios
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT2_B_all

# Generate comparison animations
make -f sample/imbh_cloud/Makefile.imbh_cloud compare_CAT_b3     # Compare CAT1 vs CAT2 vs CAT3
make -f sample/imbh_cloud/Makefile.imbh_cloud compare_RES_CAT2   # Compare A vs B vs C
```

### Individual Targets

```bash
# Simulation only (no animation)
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT2_B_b3_run

# Animation only (after simulation)
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT2_B_b3_anim

# Clean specific run
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT2_B_b3_clean
```

---

## Output Directory Structure

```
sample/imbh_cloud/
├── results/
│   ├── CAT1/                          # Adiabatic simulations
│   │   ├── A_60k/                     # Low resolution
│   │   │   ├── b1p5/                  # Impact parameter 1.5 pc
│   │   │   │   ├── snapshot_*.csv
│   │   │   │   ├── energy.csv
│   │   │   │   └── checkpoints/
│   │   │   ├── b3/                    # Impact parameter 3 pc
│   │   │   └── b6/                    # Impact parameter 6 pc
│   │   ├── B_200k/                    # Mid resolution
│   │   └── C_1M/                      # High resolution
│   ├── CAT2/                          # Radiative simulations
│   │   ├── A_60k/
│   │   ├── B_200k/
│   │   └── C_1M/
│   └── CAT3/                          # Full physics simulations
│       ├── A_60k/
│       ├── B_200k/
│       └── C_1M/
└── animations/
    └── (same structure as results)
```

---

## Recommended Workflow

### 1. Quick Parameter Exploration (Level A)
```bash
# Test all impact parameters with adiabatic physics
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT1_A_b1p5
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT1_A_b3
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT1_A_b6
```

### 2. Physics Comparison (Level B, fixed scenario)
```bash
# Compare thermal physics at b=3pc
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT1_B_b3
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT2_B_b3
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT3_B_b3
make -f sample/imbh_cloud/Makefile.imbh_cloud compare_CAT_b3
```

### 3. Resolution Convergence (fixed physics/scenario)
```bash
# Test convergence for CAT2 at b=3pc
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT2_A_b3
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT2_B_b3
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT2_C_b3
make -f sample/imbh_cloud/Makefile.imbh_cloud compare_RES_CAT2_b3
```

### 4. Publication Run (Level C)
```bash
# High-resolution with full physics
make -f sample/imbh_cloud/Makefile.imbh_cloud CAT3_C_b3
```

---

## Initial Conditions

All simulations share the same relaxed Lane-Emden sphere, rescaled to the appropriate particle count:

| Resolution | IC File | Source |
|------------|---------|--------|
| A (60k) | `hydrostatic/60k/GSPH/snapshot_0000.csv` | Relaxed polytrope |
| B (200k) | `hydrostatic/200k/GSPH/snapshot_0000.csv` | Relaxed polytrope |
| C (1M) | `hydrostatic/1M/GSPH/snapshot_0000.csv` | Relaxed polytrope |

To generate ICs:
```bash
make -f sample/imbh_cloud/Makefile.relaxation relax_60k
make -f sample/imbh_cloud/Makefile.relaxation relax_200k
make -f sample/imbh_cloud/Makefile.relaxation relax_1M
```

---

## Runtime Estimates

| Config | Resolution | Est. Runtime | Wall Clock (8 cores) |
|--------|------------|--------------|----------------------|
| CAT1_A | 60k | 1 hour | ~10 min |
| CAT1_B | 200k | 8 hours | ~1 hour |
| CAT1_C | 1M | 100 hours | ~12 hours |
| CAT2_A | 60k | 2 hours | ~15 min |
| CAT2_B | 200k | 16 hours | ~2 hours |
| CAT2_C | 1M | 200 hours | ~24 hours |
| CAT3_A | 60k | 3 hours | ~20 min |
| CAT3_B | 200k | 24 hours | ~3 hours |
| CAT3_C | 1M | 300 hours | ~36 hours |

---

## Physical Questions Addressed

| Question | Category | Resolution |
|----------|----------|------------|
| Does tidal tail form? | CAT1_A | Any |
| Effect of cooling on compression? | CAT1 vs CAT2 | B+ |
| Does cloud fragment? | CAT3 | B+ |
| Shock Mach number distribution? | CAT2 | A+ |
| P-V diagram comparison with Oka? | CAT2/CAT3 | B+ |
| Convergence of mass-loss rate? | All | A,B,C |

