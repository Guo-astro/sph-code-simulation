# SR-GSPH Animation Scripts - Solver Detection Update

## Summary

Updated both animation scripts to automatically detect and display which Riemann solver (EXACT or ITERATIVE) was used in the simulation.

## Files Updated

### 1. `scripts/animate_sr_sod.py`
- Added `detect_solver_type()` function
- Automatically detects solver from directory name or config file
- Updates plot titles to show solver type
- Output filename includes solver type: `sr_sod_animation_iterative.gif` or `sr_sod_animation_exact.gif`

### 2. `scripts/visualize_sr_sod.py`
- Added `detect_solver_type()` function
- Updates all plot titles (animation, final state, evolution) to show solver type
- Consistent labeling across all visualizations

## Detection Logic

The scripts detect the solver type in this order:

1. **Directory name** - If path contains 'iterative' or 'exact'
   - `sample/sr_sod/results/iterative/` → ITERATIVE
   - `sample/sr_sod/results/exact/` → EXACT

2. **Config file** - Reads `riemannSolverSRGSPH` from JSON config
   - Falls back to checking parent directory for `*.json` files

3. **Default** - If no indication found, defaults to EXACT

## Usage

### Automatic (via Makefile)

```bash
# Iterative solver - auto-detects and labels
make -f sample/sr_sod/Makefile.sr_sod sr_sod_iterative

# Exact solver - auto-detects and labels
make -f sample/sr_sod/Makefile.sr_sod sr_sod_exact

# Compare both - creates labeled animations
make -f sample/sr_sod/Makefile.sr_sod sr_solver_compare_viz
```

### Manual

```bash
# Animate iterative solver results
python3 sample/sr_sod/scripts/animate_sr_sod.py sample/sr_sod/results/iterative

# Animate exact solver results
python3 sample/sr_sod/scripts/animate_sr_sod.py sample/sr_sod/results/exact

# Full visualization suite
python3 sample/sr_sod/scripts/visualize_sr_sod.py sample/sr_sod/results/iterative
```

## Plot Title Format

### Before
```
SR Sod Shock Tube (Special Relativistic Godunov SPH)
```

### After (ITERATIVE)
```
SR Sod Shock Tube (SR-GSPH)
ITERATIVE Riemann Solver (van Leer Relativistic)
```

### After (EXACT)
```
SR Sod Shock Tube (SR-GSPH)
EXACT Riemann Solver (Kitajima)
```

## Output Files

### animate_sr_sod.py
- **ITERATIVE**: `sr_sod_animation_iterative.gif`
- **EXACT**: `sr_sod_animation_exact.gif`

### visualize_sr_sod.py
- Same filenames for all solvers:
  - `sr_sod_evolution.gif` (animation)
  - `sr_sod_evolution.png` (multi-time plot)
  - `sr_sod_final.png` (final state)
- But titles clearly indicate solver type

## Testing

```bash
# Test detection function
python3 << 'EOF'
import sys
sys.path.insert(0, 'sample/sr_sod/scripts')
from visualize_sr_sod import detect_solver_type

print(detect_solver_type('sample/sr_sod/results/iterative'))
# Output: ('ITERATIVE', 'van Leer Relativistic')

print(detect_solver_type('sample/sr_sod/results/exact'))
# Output: ('EXACT', 'Kitajima')
EOF
```

## Benefits

1. ✅ **No manual labeling required** - Scripts auto-detect solver
2. ✅ **Clear identification** - Plots clearly show which solver was used
3. ✅ **Easy comparison** - Different output filenames prevent confusion
4. ✅ **Backward compatible** - Works with existing results directories
5. ✅ **Consistent labeling** - All plots use same format
