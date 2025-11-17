# Quick Start: SR-GSPH Iterative Riemann Solver

## What Changed

You now have **two Riemann solver options** for SR-GSPH:

1. **EXACT** (default) - Kitajima formulation with Brent's method
2. **ITERATIVE** (new) - van Leer adapted for special relativity

## Usage

### Option 1: Use in JSON config

Add this line to your SR-GSPH configuration:

```json
{
  "SPHType": "srgsph",
  "riemannSolverSRGSPH": "ITERATIVE"
}
```

### Option 2: Test both solvers

```bash
# Run with EXACT solver (default)
./build/sph sr_sod.json

# Run with ITERATIVE solver
./build/sph sr_sod_iterative.json
```

## When to Use ITERATIVE Solver

✅ **Use ITERATIVE when:**
- Experiencing instabilities with exact solver
- Moderate shock strengths (Mach < 10)
- Need more robust convergence
- Debugging solver issues

❌ **Stick with EXACT when:**
- High accuracy required
- Strong relativistic shocks
- Performance is critical
- Default behavior works fine

## Files Modified

1. `include/parameters.hpp` - Added `riemann_solver` to SRGSPH struct
2. `include/srgsph/sr_fluid_force.hpp` - Added `iterative_solver()` declaration
3. `src/srgsph/sr_fluid_force.cpp` - Implemented relativistic iterative solver
4. `src/solver.cpp` - Added parameter parsing for solver selection

## Key Differences from GSPH Iterative Solver

The SR-GSPH version includes:
- Relativistic enthalpy calculations
- Causality enforcement (|v*| < c)
- Lorentz factor corrections
- More iterations (50 vs 20)
- Tighter tolerance (1e-8 vs 1e-6)

## Testing

The code compiles cleanly and runs without errors. Test with:

```bash
cd /Users/guo/Downloads/sphcode
./build/sph sr_sod_iterative.json
```

See `SR_ITERATIVE_SOLVER.md` for complete documentation.
