# ISM Cooling 1D Test - Quick Start Guide

## One-Line Run

```bash
cd /Users/guo/Downloads/sphcode && make -f sample/cooling_heating/Makefile.cooling_heating cooling_all
```

## What It Does

Simulates 1D density gradient (0.1 → 1000 cm⁻³) evolving toward Koyama & Inutsuka (2000) thermal equilibrium using **GDISPH + ISM Cooling**.

## Quick Commands

```bash
# Show help
make -f sample/cooling_heating/Makefile.cooling_heating cooling_help

# Run simulation only (fast, ~100ms)
make -f sample/cooling_heating/Makefile.cooling_heating cooling_run

# Generate comparison plots
make -f sample/cooling_heating/Makefile.cooling_heating cooling_plot

# Run + plot + animate
make -f sample/cooling_heating/Makefile.cooling_heating cooling_all

# Clean
make -f sample/cooling_heating/Makefile.cooling_heating cooling_clean
```

## Expected Output

```
sample/cooling_heating/results/ism_cooling_1d/
├── snapshot_0000.csv          # Initial state (t=0)
├── snapshot_0001.csv          # Final state (t=10)
├── energy.dat                 # Energy evolution
└── plots/
    ├── snapshot_0000.png      # T(n), P(n) vs K&I curves
    └── snapshot_0001.png      # Final equilibrium state
```

## Physics

### Input
- **Density**: 0.1 - 1000 cm⁻³ (log-uniform)
- **Temperature**: 8000 K (uniform, WNM-like)
- **Pressure**: Varies (P = ρT)

### Output (at thermal equilibrium)
- **S-curve**: Temperature vs density follows K&I thermal bistability
- **Two phases**: 
  - WNM (n<1): T ~ 8000 K, P ~ 2000 K cm⁻³
  - CNM (n>30): T ~ 150 K, P ~ 4000 K cm⁻³
- **Unstable region**: 1 < n < 30 cm⁻³ → particles migrate away

## Comparison with Koyama & Inutsuka (2000)

The plot scripts overlay digitized curves from the original paper:
- ✅ Temperature T(n) - Panel (a), top
- ✅ Pressure P(n) - Panel (a), bottom  
- Future: Chemical fractions x_e, x_H2, x_CO - Panel (b)

## Test Configuration

**SPH Method**: GDISPH (Riemann solver + pressure-energy + cooling)

**Resolution**: 1000 particles

**Cooling**: Hardcoded K&I interpolation (no chemistry solving)

**Column Density**: N_H = 10^19 cm⁻² (low shielding)

**Thermal Relaxation**: 0.1 code units

## Makefile Style

Follows existing conventions:
- `cooling_run` → single simulation
- `cooling_plot` → visualization
- `cooling_animate` → movie generation
- `cooling_all` → complete workflow
- `cooling_clean` → remove outputs
- `cooling_help` → usage guide

Consistent with:
- `sedov_run`, `sedov_compare_all`, `sedov_kernel_all`
- `khi_run`, `khi_compare_all`
- `gresho_run`, `gresho_compare_all`

## Build Requirements

Already installed (thermal module compiled):
- ✅ C++ compiler with C++17
- ✅ OpenMP
- ✅ HDF5 (optional)
- ✅ Boost
- ✅ nlohmann/json

Python for visualization:
- Python 3
- numpy, pandas, matplotlib

## Validation

Check thermal equilibrium convergence:
```bash
# View energy evolution
cat sample/cooling_heating/results/ism_cooling_1d/energy.dat
# Columns: time  kinetic  thermal  potential  total

# Thermal energy should stabilize as particles relax
```

## Technical Details

See `ISM_COOLING_1D_INTEGRATION.md` for:
- Complete integration guide
- Code modifications
- Physical interpretation
- Debugging notes
