# IMBH-Cloud Simulation Category System

A systematic approach to running IMBH-molecular cloud tidal disruption simulations with varying physics, resolution, and impact parameters.

## Quick Start

```bash
# Run a single low-resolution adiabatic simulation
make -f simulations/astrophysics/imbh_cloud/Makefile.categories CAT1_A_b3

# Run all physics categories at medium resolution
make -f simulations/astrophysics/imbh_cloud/Makefile.categories all_B_b3

# Compare physics (adiabatic vs radiative vs full)
make -f simulations/astrophysics/imbh_cloud/Makefile.categories compare_physics_b3
```

## Category System Overview

### Naming Convention: `CAT{physics}_{resolution}_{scenario}`

| Component | Options | Description |
|-----------|---------|-------------|
| **CAT** | 1, 2, 3 | Physics category |
| **Resolution** | A, B, C | Particle count |
| **Scenario** | b1p5, b3, b6 | Impact parameter |

Example: `CAT2_B_b3` = Radiative physics, 200k particles, b=3pc

## Physics Categories

| Category | Name | Physics | Use Case |
|----------|------|---------|----------|
| **CAT1** | Adiabatic | No cooling, γ=5/3 EOS | Baseline dynamics |
| **CAT2** | Radiative | Inoue-Inutsuka cooling | ISM thermal physics |
| **CAT3** | Full Physics | Radiative + enhanced self-gravity | Fragmentation studies |

### Physics Details

**CAT1 (Adiabatic)**
- Pure hydrodynamics with γ=5/3 polytropic equation of state
- No cooling or heating
- Fastest runtime
- Best for: Understanding pure tidal dynamics, shock heating

**CAT2 (Radiative)**
- Inoue & Inutsuka (2008) ISM cooling function:
  ```
  Λ/Γ = 10^7 exp(-114800/(T+1000)) + 0.014 √T exp(-92/T)
  ```
- Thermal equilibrium at T ≈ 60-80 K
- Best for: Realistic ISM physics, density enhancement

**CAT3 (Full Physics)**
- Radiative cooling + enhanced self-gravity (θ=0.3)
- Lower Barnes-Hut opening angle for accurate fragmentation
- Slowest runtime
- Best for: Fragmentation, bound clump formation

## Resolution Levels

| Level | Particles | Smoothing Length | Runtime (est.) | Use Case |
|-------|-----------|------------------|----------------|----------|
| **A** | 64,000 | ~0.03 pc | ~10 min | Quick tests |
| **B** | 200,000 | ~0.02 pc | ~1-3 hours | Production |
| **C** | 1,000,000 | ~0.01 pc | ~12-36 hours | Publication |

## Impact Parameter Scenarios

| Scenario | b (pc) | b/r_tidal | Regime |
|----------|--------|-----------|--------|
| **b1p5** | 1.5 | 0.41 | Deep disruption |
| **b3** | 3.0 | 0.83 | Partial disruption |
| **b6** | 6.0 | 1.65 | Weak interaction |

## Makefile Targets

### Single Simulation
```bash
make -f simulations/astrophysics/imbh_cloud/Makefile.categories CAT1_A_b3        # Full (sim + anim)
make -f simulations/astrophysics/imbh_cloud/Makefile.categories CAT1_A_b3_run    # Simulation only
make -f simulations/astrophysics/imbh_cloud/Makefile.categories CAT1_A_b3_anim   # Animation only
```

### Batch Runs
```bash
# All resolutions for one category/scenario
make -f simulations/astrophysics/imbh_cloud/Makefile.categories CAT1_all_b3

# All categories for one resolution/scenario
make -f simulations/astrophysics/imbh_cloud/Makefile.categories all_B_b3
```

### Comparison Workflows
```bash
# Physics comparison (CAT1 vs CAT2 vs CAT3) at Resolution B
make -f simulations/astrophysics/imbh_cloud/Makefile.categories compare_physics_b3

# Resolution convergence (A vs B vs C) for CAT2
make -f simulations/astrophysics/imbh_cloud/Makefile.categories compare_resolution_CAT2_b3
```

### Cleanup
```bash
make -f simulations/astrophysics/imbh_cloud/Makefile.categories clean_CAT1
make -f simulations/astrophysics/imbh_cloud/Makefile.categories clean_all
```

## Directory Structure

```
simulations/astrophysics/imbh_cloud/
├── config/presets/categories/
│   ├── CATEGORY_SYSTEM.md           # Full documentation
│   ├── CAT1/                        # Adiabatic configs
│   │   ├── A_60k/
│   │   │   ├── b1p5.json
│   │   │   ├── b3.json
│   │   │   └── b6.json
│   │   ├── B_200k/
│   │   └── C_1M/
│   ├── CAT2/                        # Radiative configs
│   └── CAT3/                        # Full physics configs
├── results/
│   ├── CAT1/
│   │   ├── A_60k/
│   │   │   ├── b1p5/
│   │   │   │   ├── snapshot_*.csv
│   │   │   │   ├── energy.csv
│   │   │   │   └── checkpoints/
│   │   │   ├── b3/
│   │   │   └── b6/
│   │   ├── B_200k/
│   │   └── C_1M/
│   ├── CAT2/
│   └── CAT3/
├── animations/
│   └── (same structure as results)
├── scripts/
│   └── generate_category_configs.py  # Config generator
└── Makefile.categories               # Category system Makefile
```

## Initial Conditions

Each resolution requires its own relaxed Lane-Emden sphere. Generate them with:

```bash
make -f simulations/astrophysics/imbh_cloud/Makefile.relaxation relax_60k
make -f simulations/astrophysics/imbh_cloud/Makefile.relaxation relax_200k
make -f simulations/astrophysics/imbh_cloud/Makefile.relaxation relax_1M
```

## Configuration Generation

To regenerate all JSON configs:
```bash
python simulations/astrophysics/imbh_cloud/scripts/generate_category_configs.py
```

To generate specific configs:
```bash
python simulations/astrophysics/imbh_cloud/scripts/generate_category_configs.py --category CAT2 --resolution B
```

## Recommended Workflows

### 1. Quick Parameter Exploration (Level A)
```bash
# Test all impact parameters with adiabatic physics (~30 min total)
make -f simulations/astrophysics/imbh_cloud/Makefile.categories CAT1_A_b1p5
make -f simulations/astrophysics/imbh_cloud/Makefile.categories CAT1_A_b3
make -f simulations/astrophysics/imbh_cloud/Makefile.categories CAT1_A_b6
```

### 2. Physics Comparison Study
```bash
# Compare CAT1/CAT2/CAT3 at b=3pc with production resolution (~6 hours total)
make -f simulations/astrophysics/imbh_cloud/Makefile.categories compare_physics_b3
```

### 3. Resolution Convergence
```bash
# Test A/B/C for CAT2 at b=3pc
make -f simulations/astrophysics/imbh_cloud/Makefile.categories compare_resolution_CAT2_b3
```

### 4. Publication Run
```bash
# High-resolution full physics simulation (~36 hours)
make -f simulations/astrophysics/imbh_cloud/Makefile.categories CAT3_C_b3
```

## Physical Questions by Category

| Question | Best Category | Min Resolution |
|----------|---------------|----------------|
| Does tidal tail form? | CAT1 | A |
| What is the shock Mach number? | CAT1/CAT2 | A |
| How does cooling affect compression? | CAT2 vs CAT1 | B |
| Does the cloud fragment? | CAT3 | B+ |
| Can we reproduce Oka et al. P-V diagram? | CAT2/CAT3 | B |
| What is the mass-loss rate? | All | B (convergence) |

## Output Files

Each simulation produces:
- `snapshot_XXXX.csv` - Particle data at each output time
- `energy.csv` - Energy conservation diagnostics
- `checkpoints/` - Restart files

Each animation set produces:
- `tidal.gif` - Main tidal disruption animation
- `shock.gif` - Shock diagnostics
- `density.gif` - Density field
- `3d.gif` - 3D view with black hole

## Tips

1. **Start with Level A** to quickly explore parameter space
2. **Use Level B** for production runs and physics comparisons
3. **Reserve Level C** for final publication-quality results
4. **CAT1 is fastest** - use it to debug initial conditions
5. **CAT3 requires good resolution** - Level A may miss fragmentation
