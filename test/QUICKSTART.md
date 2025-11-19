# Quick Start: Debug Iterative Riemann Solver

## Problem Identified

The C++ iterative Riemann solver in `sr_fluid_force.cpp` **fails for the SR Sod shock tube test** due to an inappropriate initial guess.

## Quick Diagnostic

```bash
cd /Users/guo/Downloads/sphcode/test
python3 quick_diagnostic.py
```

Output shows:
- **Initial guess**: P_guess = 12.57 (15x too high!)
- **Correct answer**: P* = 0.833
- **Newton step**: Would produce negative pressure (-11.67)
- **Root cause**: Two-rarefaction approximation fails for extreme pressure ratios (P_L/P_R > 1e9)

## Full Test Suite

```bash
cd /Users/guo/Downloads/sphcode/test
./run_riemann_debug.sh
```

This will:
1. Compile C++ test program
2. Run 5 test cases in both C++ and Python
3. Generate comparison plots in `comparison_plots/`

## Files Created

| File | Purpose |
|------|---------|
| `test_iterative_riemann_solver.cpp` | Standalone C++ test with debug output |
| `compare_riemann_solvers.py` | Python comparison using kitajima_solver.py |
| `visualize_comparison.py` | Generate comparison plots |
| `quick_diagnostic.py` | Single-test detailed analysis |
| `run_riemann_debug.sh` | Automated workflow |
| `RIEMANN_DEBUG_README.md` | Full documentation |
| `DEBUGGING_SUMMARY.md` | Detailed findings report |

## The Bug

### Current Behavior (WRONG)
```
SR Sod Test:
  C++:    P* = 1.0e-8  (vacuum solution - WRONG!)
  Python: P* = 0.833   (correct)
  Error:  833 billion %
```

### Root Cause

In `sr_fluid_force.cpp`, line ~246:
```cpp
// Initial guess using two-rarefaction approximation
const real exponent = (gamma_c - 1.0) / (2.0 * gamma_c);
const real num = csl + csr - 0.5 * (gamma_c - 1.0) * (vr - vl);
const real denom = csl / std::pow(Pl, exponent) + csr / std::pow(Pr, exponent);
real pguess = std::pow(num / denom, 2.0 * gamma_c / (gamma_c - 1.0));
```

For SR Sod: This gives P_guess = 12.57, but correct P* = 0.833
- The guess is **15x too high**
- Newton-Raphson takes a huge first step to P = 1e-9 (near P_R)
- Solver converges to wrong vacuum solution

### Recommended Fix

Replace two-rarefaction guess with robust PVRS method:

```cpp
// PVRS (Primitive Variable Riemann Solver) initial guess
// More robust for extreme pressure ratios
real pguess;
if (Pl/Pr > 1e6 || Pr/Pl > 1e6) {
    // Extreme ratio: use geometric mean
    pguess = std::sqrt(Pl * Pr);
} else {
    // Normal case: use two-rarefaction
    const real exponent = (gamma_c - 1.0) / (2.0 * gamma_c);
    const real num = csl + csr - 0.5 * (gamma_c - 1.0) * (vr - vl);
    const real denom = csl / std::pow(Pl, exponent) + csr / std::pow(Pr, exponent);
    pguess = std::pow(num / denom, 2.0 * gamma_c / (gamma_c - 1.0));
}

// Add safety bounds
pguess = std::max(pguess, std::min(Pl, Pr) * 0.1);
pguess = std::min(pguess, std::max(Pl, Pr) * 10.0);
```

## Test Results

### With Current Code (BROKEN)

| Test | C++ P* | Python P* | Error |
|------|--------|-----------|-------|
| SR_Sod_Shock_Tube | 1.0e-8 | 0.833 | 833 billion % ❌ |
| Two_Rarefactions | 0.2497 | 0.125 | 50% ❌ |
| Two_Shocks | (not shown) | 0.25 | ? |

### Expected After Fix

All tests should show < 1e-10 relative error

## Visualization

After running `./run_riemann_debug.sh`, check:

```
comparison_plots/
├── summary_comparison.png           # Overview
├── SR_Sod_Shock_Tube_comparison.png # Shows the bug clearly
└── Two_Rarefactions_comparison.png  # Also shows error
```

## Next Steps

1. **Implement the fix** in `src/srgsph/sr_fluid_force.cpp`
2. **Re-run the test suite**: `./run_riemann_debug.sh`
3. **Verify all tests pass** with < 1e-10 error
4. **Run full SR Sod simulation** to confirm physical correctness

## References

- Pons et al. (1999): Relativistic Riemann problem with tangential velocities
- Kitajima et al. (2025): arXiv:2510.18251v1
- Toro (2009): "Riemann Solvers and Numerical Methods" - Chapter 4 on initial guesses

---

**Author**: Debugging suite by GitHub Copilot  
**Date**: November 18, 2025  
**Status**: Bug identified, fix recommended
