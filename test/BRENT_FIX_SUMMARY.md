# Brent's Method Fix for Iterative Riemann Solver

## Problem Identified

The Newton-Raphson method in `sr_fluid_force.cpp` **catastrophically failed** for extreme pressure ratios:

- **SR Sod test**: 8.3 billion % error (converged to P*=1e-8 instead of 0.833)
- **Two Rarefactions**: 50% error  
- **Two Shocks**: 93% error
- **Left Shock + Right Rarefaction**: 42% error
- **With Tangential Velocity**: 34% error

### Root Cause

Newton-Raphson requires a **good initial guess**. For extreme pressure ratios (P_L/P_R > 1e6):
- Two-rarefaction formula gives: `P_guess = sqrt(13.33 * 1e-8) = 0.000365`
- Correct answer: `P* = 0.833`
- First Newton iteration: `P_new = -4.4` → clamped to 1e-9 → **converged to vacuum!**

Even "smart" adaptive initial guesses (geometric/arithmetic mean, safety bounds) **failed** because:
1. Newton can still go negative even from reasonable starting points
2. Clamping to stay positive creates artificial local minima
3. Under-relaxation helps but doesn't guarantee convergence

## Solution Implemented

### Replaced Newton-Raphson with **Brent's Method + Adaptive Bracketing**

Following the robust approach from `kitajima_solver.py`:

#### Step 1: Adaptive Bracketing
```cpp
// Start with arithmetic mean
real p_mid = 0.5 * (P_L + P_R);
real p_min = p_mid, p_max = p_mid;

// Expand bracket until dvel changes sign
for (int i = 0; i < 100; ++i) {
    p_min = 0.5 * max(p_min, 0.0);
    p_max = 2.0 * p_max;
    
    dvel_min = compute_dvel(p_min);
    dvel_max = compute_dvel(p_max);
    
    if (dvel_min * dvel_max <= 0.0) {
        bracket_found = true;
        break;
    }
}
```

#### Step 2: Brent's Root-Finding
```cpp
// Use Brent's method (combination of bisection + inverse quadratic interpolation)
// Guaranteed to converge within bracket [p_min, p_max]
// No derivatives needed, robust for poorly-behaved functions
```

### Why This Works

✅ **Bracketing guarantees root exists** (if dvel changes sign)  
✅ **Brent's method always converges** (combines bisection safety with interpolation speed)  
✅ **No initial guess needed** (only bracket required)  
✅ **No derivatives** (avoids finite-difference errors)  
✅ **Can't converge to vacuum** (bracket prevents it)  

## Files Modified

### `/Users/guo/Downloads/sphcode/src/srgsph/sr_fluid_force.cpp`

Lines ~390-550: Replaced entire Newton-Raphson section with:
1. Adaptive bracketing loop (lines 390-430)
2. Brent's method implementation (lines 455-550)
3. Final solution extraction (lines 555-585)

Key changes:
- **Removed**: `pguess`, two-rarefaction formula, Newton-Raphson iteration
- **Added**: `dvel_at_p` lambda, bracket expansion, Brent's algorithm
- **Preserved**: Same interface (inputs: P_L, P_R, v_L, v_R, etc.; outputs: P*, v*)

## Verification Status

### ✅ SPH Code (`sr_fluid_force.cpp`)
- Compiles successfully with Brent's method
- Simulations run without errors
- **Needs**: Full SR Sod run to t=0.35 to verify correctness

### ⚠️  Test Program (`test_iterative_riemann_solver.cpp`)
- **Still uses old Newton-Raphson code**
- Shows same 5/5 test failures
- **Needs**: Update to match new Brent's method OR separate verification

## Recommended Next Steps

1. **Update test program** to use Brent's method (for direct C++ vs Python comparison)
2. **Run full SR Sod simulation** and check final profiles match Python reference
3. **Benchmark performance**: Brent's method may be slightly slower than Newton (but much more robust)
4. **Test edge cases**: vacuum boundaries, strong shocks, ultra-relativistic flows

## Reference Implementation

The Python solver (`kitajima_solver.py`) uses this exact approach:
- Lines 705-721: Adaptive bracketing
- Lines 483-544: Brent's method
- This is the "gold standard" that our C++ should match

## Performance Notes

**Expected tradeoffs**:
- Brent's method: ~5-10 iterations typical, guaranteed convergence
- Newton-Raphson: ~3-5 iterations when it works, **diverges catastrophically** when it doesn't

For production SPH code, **robustness >> speed** at the particle-pair level.
