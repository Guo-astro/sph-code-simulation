# Quick Start: Testing GSPH Riemann Solvers

This guide shows how to easily switch between HLL and iterative Riemann solvers in GSPH.

## Quick Test

### Test HLL Solver (default, fast)
```bash
./build/sph sample/shock_tube/config/presets/shock_tube_1d_gsph_hll.json
```

### Test Iterative Solver (accurate, slower)
```bash
./build/sph sample/shock_tube/config/presets/shock_tube_1d_gsph_iterative.json
```

## Configuration

In your JSON config file, add:

```json
{
  "SPHType": "gsph",
  "riemannSolver": "hll"        // or "iterative"
}
```

**Options:**
- `"hll"` - Fast approximate solver (default if omitted)
- `"iterative"` - Accurate iterative solver from van Leer (1997)

## Verify Which Solver Is Running

Check the output log for:
```
SPH type: Godunov SPH (2nd order)
* Riemann solver: HLL
```
or
```
SPH type: Godunov SPH (2nd order)
* Riemann solver: Iterative (van Leer 1997)
```

## Performance Comparison

| Solver    | Speed | Accuracy | Best Use Case                |
|-----------|-------|----------|------------------------------|
| HLL       | Fast  | Good     | Production runs, general use |
| Iterative | Slower| Excellent| Validation, high-precision   |

## Example: Shock Tube Test Comparison

```bash
# Run with HLL
./build/sph sample/shock_tube/config/presets/shock_tube_1d_gsph_hll.json

# Run with Iterative
./build/sph sample/shock_tube/config/presets/shock_tube_1d_gsph_iterative.json

# Compare results
ls -lh sample/shock_tube/results/gsph_hll/
ls -lh sample/shock_tube/results/gsph_iterative/
```

## Creating Your Own Config

Copy an existing GSPH config and add the `riemannSolver` field:

```bash
# Copy base config
cp sample/shock_tube/config/presets/shock_tube_1d_gsph.json my_test.json

# Edit to add solver selection
nano my_test.json
```

Add:
```json
{
  ...
  "SPHType": "gsph",
  "use2ndOrderGSPH": true,
  "riemannSolver": "iterative",    // <-- Add this line
  ...
}
```

## Troubleshooting

### Solver not recognized
**Error:** Default HLL is used even with `"riemannSolver": "iterative"`

**Solution:** Check spelling - must be exactly `"iterative"` or `"hll"` (lowercase)

### Slow performance
**Issue:** Simulation is slower than expected

**Solution:** The iterative solver is 10-30% slower. Use HLL for production runs.

### Verification
```bash
# Check which solver is actually being used
./build/sph your_config.json 2>&1 | grep "Riemann solver"
```

## Implementation Details

See `docs/GSPH_RIEMANN_SOLVER.md` for:
- Full algorithm description
- Mathematical equations
- Performance analysis
- Code architecture

## References

- Cha & Whitworth (2003) - "Implementations and tests of Godunov-type particle hydrodynamics"
- van Leer (1997) - "Towards the Ultimate Conservative Difference Scheme"
