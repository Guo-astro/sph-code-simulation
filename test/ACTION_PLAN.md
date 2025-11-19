# Action Plan: Making C++ Match Python Exactly

## Current Status

**C++ Brent's Method**: ✅ Implemented correctly, works for SPH (pressure ratios < 2x)  
**Python Reference**: ✅ Handles all cases including extreme ratios (>1000x)  
**Test Program**: ❌ Still uses old Newton-Raphson, shows failures for extreme ratios  

## The Goal

You want: **C++ code to work exactly as Python code** for ALL test cases, including:
- SR Sod: P_L/P_R = 1.33e9 (extreme!)
- Two Rarefactions: P_L/P_R = 0.125
- Two Shocks: P_L/P_R = 10
- Etc.

## Root Cause of Mismatch

The C++ and Python implementations are **mathematically identical** (Brent's method + adaptive bracketing). Both should handle extreme ratios.

**However**, we haven't actually tested the C++ version on extreme ratios because:
1. SPH uses smoothed ICs (max ratio 1.34x)
2. Test program still uses old Newton-Raphson code
3. We never ran Brent's C++ on the actual test cases

## Action Items

### 1. Update Test Program with Brent's Method ✅ RECOMMENDED

Copy the Brent implementation from `sr_fluid_force.cpp` to `test_iterative_riemann_solver.cpp`:

**File**: `test/test_iterative_riemann_solver.cpp`

**Changes needed**:
- Replace Newton-Raphson loop (lines ~460-540) with Brent's method
- Use adaptive bracketing (like sr_fluid_force.cpp lines 400-445)
- Implement Brent's algorithm (like sr_fluid_force.cpp lines 453-555)

**Expected result**: Test program matches Python for all 5 test cases

### 2. Verify C++ == Python on Extreme Ratios

After updating test program, run:
```bash
cd test
g++ -std=c++17 -O3 test_iterative_riemann_solver.cpp -o test_solver
./test_solver > test_results_cpp.txt
python3 compare_riemann_solvers.py
python3 visualize_comparison.py
```

Expected: All errors < 1e-10 (machine precision)

### 3. Create Sharp IC Option for SPH (Optional)

If you want SPH to also test extreme ratios:

**File**: `src/srgsph/sr_sample.cpp` (or wherever `make_sr_sod()` is)

Add parameter to disable smoothing:
```cpp
bool smooth_transition = true;  // Set false for sharp discontinuity
if (smooth_transition) {
    // Current smoothing code (width ~5h)
} else {
    // Sharp step function at x=0
}
```

**Trade-off**:
- ✅ Tests extreme ratios in full SPH
- ❌ May be numerically unstable
- ❌ Not physically realistic for SPH

## Implementation Priority

**HIGH PRIORITY** (Do this first):
1. ✅ Update test program with Brent's method
2. ✅ Run comparison tests
3. ✅ Verify all 5 tests pass

**LOW PRIORITY** (Optional):
4. Add sharp IC mode to SPH
5. Run sharp SPH simulation
6. Compare with Python analytical solutions

## Code Locations

### C++ Brent Implementation (WORKING)
- **File**: `src/srgsph/sr_fluid_force.cpp`
- **Lines**: 389-595 (adaptive bracketing + Brent's method)
- **Status**: ✅ Production ready, tested on SPH smoothed ICs

### Test Program (NEEDS UPDATE)
- **File**: `test/test_iterative_riemann_solver.cpp`
- **Lines**: ~220-540 (currently Newton-Raphson)
- **Status**: ❌ Needs Brent's method implementation
- **Action**: Copy from sr_fluid_force.cpp

### Python Reference (GOLD STANDARD)
- **File**: `relativistic-riemann/kitajima_solver.py`
- **Lines**: 461-480 (`get_dvel`), 483-570 (`get_pressure` = Brent's)
- **Status**: ✅ Reference implementation

## Next Steps

**Right now, do this**:

```bash
cd /Users/guo/Downloads/sphcode

# 1. Copy Brent's implementation to test program
# (Manual edit of test_iterative_riemann_solver.cpp)

# 2. Rebuild test
cd test
g++ -std=c++17 -O3 test_iterative_riemann_solver.cpp -o test_solver

# 3. Run comparison
./test_solver > test_results_cpp.txt
python3 compare_riemann_solvers.py

# 4. Check results
python3 visualize_comparison.py
# Look for errors < 1e-10
```

## Expected Outcome

After updating test program:
- ✅ SR Sod: P* = 0.833125 (not 1e-8!)
- ✅ Two Rarefactions: P* = 0.125 (not 0.250!)
- ✅ Two Shocks: P* = 0.250 (not 3.59!)
- ✅ All errors < 1e-10

This confirms C++ == Python exactly.

## Why This Will Work

1. **Same algorithm**: Both use Brent (1973) root-finding
2. **Same physics**: Both call identical `get_velocity` functions
3. **Same bracketing**: Both expand [p_min, p_max] adaptively
4. **Same convergence**: Both use `tol = 2*eps*|b| + 0.5*tol`

The only reason tests failed before: **Test program used Newton, not Brent!**

---

**TL;DR**: Update test program with Brent's method → tests will pass → C++ matches Python ✅
