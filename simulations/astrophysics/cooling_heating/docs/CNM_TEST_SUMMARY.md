# CNM Thermal Relaxation Test - Implementation Summary

## Overview
Successfully implemented CNM (Cold Neutral Medium) thermal relaxation test based on **Koyama & Inutsuka (2000) Model C10**.

## What Changed From Previous Version

### 1. Initial Conditions - From Density Gradient to Uniform CNM
**Old approach (incorrect):**
- Log-uniform density distribution: 0.1 - 1000 cm⁻³
- Uniform temperature: 8000 K (WNM-like)
- Variable mass particles (mass ∝ density)

**New approach (correct CNM test):**
- **Uniform density**: n = 10 cm⁻³ (CNM)
- **Uniform base temperature**: T = 107 K (±10% random perturbations)
- **Equal mass particles**: mass = ρ × dx (like Sod shock tube)

### 2. Kernel - From Wendland to Cubic Spline
**Problem:** Wendland kernel was returning w=0 for all distances (including r=0)
**Solution:** Switched to cubic_spline kernel
**Result:** Proper kernel values (w=133.333 at r=0, decreasing with distance)

### 3. Physical Units - Simplified Code Units
**Approach:** 
- density_to_n_H = 1.0 (set in config)
- ρ_code = n_H directly (number density)
- T_code = T_K directly (temperature in Kelvin)
- P = n_H × T (ideal gas, code units)

## Test Parameters (K&I Model C10)

| Parameter | Value | Description |
|-----------|-------|-------------|
| Initial density | 10 cm⁻³ | Uniform CNM |
| Initial temperature | 107 K | With ±10% random perturbations |
| Column density | 10²⁰ cm⁻² | For photoelectric heating |
| Domain | [0, 1] | 1D periodic box |
| Particles | 1000 | Equal mass initialization |
| Spacing | dx = 0.001 | Uniform |
| Mass | m = 0.01 | Equal for all particles |
| Gamma | 5/3 | Monoatomic gas |
| Kernel | cubic_spline | Fixed Wendland issue |
| Neighbors | 20 | Adequate for 1D |

## Physics

### Equal Mass Initialization (Sod Shock Tube Style)
```cpp
const real dx = domain_size / N;           // Uniform spacing: 0.001
const real mass = rho_code * dx;           // Equal mass: 10 × 0.001 = 0.01
const real rho_code = 10.0;                // Uniform density

for (int i = 0; i < N; ++i) {
    p[i].pos[0] = (i + 0.5) * dx;          // [0.0005, 0.0015, ..., 0.9995]
    p[i].mass = mass;                       // 0.01 (all particles equal)
    p[i].dens = rho_code;                   // 10 (will be smoothed)
    
    // Temperature perturbation: ±10%
    real T_pert = 107.0 * (0.9 + 0.2 * random);  // [96.3, 117.7] K
    p[i].pres = rho_code * T_pert;          // P = n × T
    p[i].ene = p[i].pres / ((gamma-1) * rho_code);
}
```

### Key Differences From Density Gradient Approach

| Aspect | Old (Gradient) | New (Uniform CNM) |
|--------|----------------|-------------------|
| Density | Log-uniform 0.1-1000 | Uniform 10 cm⁻³ |
| Mass | Variable (∝ density) | Equal (Sod-style) |
| Temperature | Uniform 8000 K | 107 K ± 10% |
| Purpose | Sample full phase space | Test thermal relaxation |
| K&I Model | None (custom) | Model C10 |

## What the Test Does

1. **Initialize**: 1000 particles with uniform density (10 cm⁻³), slightly perturbed temperatures (~107 K ±10%)

2. **Thermal Evolution**: Particles cool/heat via:
   - Photoelectric heating (dust grains)
   - Fine-structure line cooling (C II, O I, Si II, Fe II)
   - H₂ formation/dissociation
   - Lyman-α cooling
   - Gas-grain collisions

3. **Relaxation**: System evolves toward thermal equilibrium where heating = cooling

4. **Compare**: Final temperature distribution vs K&I Figure 1a equilibrium curve at n=10 cm⁻³

## Expected Results

For **n=10 cm⁻³** and **N_H = 10²⁰ cm⁻²**, K&I Figure 1a shows:
- **Equilibrium temperature**: T_eq ≈ 20-30 K (cold neutral medium)
- **Pressure**: P/k_B ≈ 2×10³ K cm⁻³

Particles should **cool** from T=107 K → T≈20-30 K over thermal timescale.

## Usage

```bash
# Show help
make -f simulations/astrophysics/cooling_heating/Makefile.cooling cnm_help

# Run simulation
make -f simulations/astrophysics/cooling_heating/Makefile.cooling cnm_run

# Generate plots
make -f simulations/astrophysics/cooling_heating/Makefile.cooling cnm_plot

# Create animation
make -f simulations/astrophysics/cooling_heating/Makefile.cooling cnm_animate

# One-shot: run + plot + animate
make -f simulations/astrophysics/cooling_heating/Makefile.cooling cnm_all

# Clean results
make -f simulations/astrophysics/cooling_heating/Makefile.cooling cnm_clean
```

## Files

```
simulations/astrophysics/cooling_heating/
├── Makefile.cooling                    # CNM test targets
├── cnm_relaxation.json                 # Configuration (Model C10)
├── scripts/
│   ├── plot_cnm_relaxation.py         # TODO: Create plotting script
│   └── animate_cnm_relaxation.py      # TODO: Create animation script
└── results/cnm_relaxation/
    ├── snapshot_*.csv                  # Simulation snapshots
    ├── energy.dat                      # Energy evolution
    └── plots/                          # Visualizations (TODO)
```

## Validation Checklist

✅ Equal mass particles initialized correctly
✅ Uniform density (dens=10 cm⁻³) maintained after initial_smoothing
✅ Temperature perturbations applied (96-117 K around 107 K)
✅ Cubic spline kernel working (w=133.333 at r=0)
✅ Periodic boundaries working
✅ Simulation running without NaN values
✅ Configuration matches K&I Model C10 parameters

⏳ TODO: Create visualization scripts
⏳ TODO: Compare final state with K&I Figure 1a
⏳ TODO: Measure thermal relaxation timescale

## References

- **Koyama & Inutsuka (2000)**: "Molecular Cloud Formation in Shock-Compressed Layers", ApJ 533, 793
- **Model C10**: Table 1, n_i=10 cm⁻³, T_i=107 K, N_H=10²⁰ cm⁻²
- **Figure 1a**: Equilibrium T(n) and P(n) curves for comparison

## Next Steps

1. Let simulation run to completion (endTime = 100 code units)
2. Create `plot_cnm_relaxation.py` to compare with K&I curves
3. Measure cooling timescale: τ_cool ~ t for T to reach T_eq
4. Verify equilibrium temperature matches K&I Figure 1a at n=10 cm⁻³
