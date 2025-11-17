# SR-GSPH Riemann Solver Comparison Targets

## New Makefile Targets

### Quick Commands

```bash
# Run with ITERATIVE solver (includes animation)
make -f sample/sr_sod/Makefile.sr_sod sr_sod_iterative

# Run with EXACT solver (includes animation)
make -f sample/sr_sod/Makefile.sr_sod sr_sod_exact

# Compare both solvers (no animation)
make -f sample/sr_sod/Makefile.sr_sod sr_solver_compare

# Compare both solvers WITH animations
make -f sample/sr_sod/Makefile.sr_sod sr_solver_compare_viz
```

## New Preset Configurations

Created in `/Users/guo/Downloads/sphcode/sample/sr_sod/config/presets/`:

1. **`sr_sod_iterative.json`** - Uses `"riemannSolverSRGSPH": "ITERATIVE"`
   - Output: `sample/sr_sod/results/iterative/`
   - Uses van Leer relativistic iterative solver
   
2. **`sr_sod_exact.json`** - Uses `"riemannSolverSRGSPH": "EXACT"`
   - Output: `sample/sr_sod/results/exact/`
   - Uses Kitajima exact solver (default)

## Workflow Examples

### Compare Solvers Side-by-Side

```bash
# Run both and generate animations
make -f sample/sr_sod/Makefile.sr_sod sr_solver_compare_viz

# Results will be in:
# - sample/sr_sod/results/exact/
# - sample/sr_sod/results/iterative/
```

### Test Iterative Solver Standalone

```bash
# Quick test with animation
make -f sample/sr_sod/Makefile.sr_sod sr_sod_iterative

# Animation will be in:
# sample/sr_sod/results/iterative/sr_sod_evolution.gif
```

### Manual Comparison

```bash
# Run exact solver
make -f sample/sr_sod/Makefile.sr_sod sr_sod_exact

# Run iterative solver
make -f sample/sr_sod/Makefile.sr_sod sr_sod_iterative

# Compare animations:
open sample/sr_sod/results/exact/sr_sod_evolution.gif
open sample/sr_sod/results/iterative/sr_sod_evolution.gif
```

## What Gets Generated

Each target automatically:
1. ✅ Copies the appropriate preset configuration
2. ✅ Creates output directory
3. ✅ Runs the simulation
4. ✅ Generates visualization/animation (for individual targets)
5. ✅ Reports file locations

## Updated Help Menu

View all available targets:
```bash
make -f sample/sr_sod/Makefile.sr_sod sr_help
```

## Integration with Existing Workflow

The new targets work alongside existing ones:
- `sr_sod_all` - Still uses default configuration
- `sr_all_tests` - Runs all standard tests
- `sr_ultra_comparison` - Ultra-relativistic velocity comparison

All solver comparison targets are new additions.
