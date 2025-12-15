# 3D Cloud Cooling Hydrostatic Test - Inoue & Inutsuka (2008)

Test configuration for validating the ISM cooling implementation on a 3D Lane-Emden
cloud, matching the setup used in IMBH-cloud encounter simulations.

## Physics Background

The Inoue & Inutsuka (2008) cooling function:
- **Net cooling**: `ρL = n(-Γ + nΛ)` [erg cm⁻³ s⁻¹]
- **Heating rate**: `Γ = 2×10⁻²⁶` erg s⁻¹ (constant)
- **Cooling coefficient**: `Λ/Γ = 10⁷ exp(-114800/(T+1000)) + 0.014√T exp(-92/T)`

### Equilibrium Phases

| Phase | Number Density | Temperature | Pressure |
|-------|---------------|-------------|----------|
| WNM   | n ~ 0.57 cm⁻³ | T ~ 6000 K  | P/k ~ 3500 K cm⁻³ |
| CNM   | n ~ 30 cm⁻³   | T ~ 100 K   | P/k ~ 3000 K cm⁻³ |

## Test: `cloud_cooling_61k.json`

**Purpose**: Test hydrostatic equilibrium of a 3D Lane-Emden cloud with ISM cooling.

- **Initial**: Lane-Emden polytrope (same as IMBH simulations)
- **Expected behavior**:
  - Dense core evolves toward CNM (T ~ 100 K)
  - Diffuse envelope evolves toward WNM (T ~ 6000 K)
  - Cloud maintains structural integrity

## Usage

```bash
# Run test
make -f simulations/astrophysics/imbh_cloud/Makefile.hydrostatic cooling_cloud

# Visualize only (after run)
make -f simulations/astrophysics/imbh_cloud/Makefile.hydrostatic cooling_cloud_viz

# Clean results
make -f simulations/astrophysics/imbh_cloud/Makefile.hydrostatic cooling_clean

# Show help
make -f simulations/astrophysics/imbh_cloud/Makefile.hydrostatic cooling_help
```

## Prerequisites

1. **DIM=3 build**:
   ```bash
   cd build && cmake -DDIM=3 .. && make -j
   ```

2. **Initial conditions** from relaxation:
   ```bash
   # Generate IC first
   make -f simulations/astrophysics/imbh_cloud/Makefile.relaxation relax_oneshot SIZE=61k
   ```

## Output

```
results/hydrostatic/cooling/cloud_61k/
├── snapshot_*.csv
├── energy.dat
└── run.log

animations/hydrostatic/cooling/cloud_61k/
└── (visualization outputs)
```

## Validation Criteria

- Core temperature decreases toward CNM (~100 K)
- Envelope temperature increases toward WNM (~6000 K)
- Cloud structure maintained (no collapse/explosion)
- Energy changes reflect cooling physics (thermal energy decreases)

## References

1. Inoue, T. & Inutsuka, S. (2008), ApJ
   "Two-Fluid MHD Simulations of Converging HI Flows in the ISM"
