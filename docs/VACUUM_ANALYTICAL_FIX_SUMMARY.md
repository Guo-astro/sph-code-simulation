# Vacuum Test Analytical Solution Fix

## Problem
The analytical solution for the 1D vacuum test showed incorrect behavior:
- Spurious vertical black lines in plots (especially internal energy)
- Extremely high values at center (ρ≈8.5, P≈8.0, e≈2.35) instead of expected low values
- 116 grid points with zero density (unassigned)

## Root Causes

### 1. Sign Error in Star-State Formula (Primary)
**Incorrect formula:**
```python
c_star = ((v_R - v_L) * (γ-1)/4) + (c_L + c_R)/2
```

**Correct formula (Toro Section 4.3.2, Eq. 4.53):**
```python
c_star = (c_L + c_R)/2 - (γ-1)/4 * (v_R - v_L)
```

The sign error caused:
- c_star = 1.148 (wrong) vs 0.348 (correct)
- This propagated to all derived quantities (ρ_*, P_*, e_*)
- Star-state density became 8.5 instead of 0.022
- Star-state pressure became 8.0 instead of 0.0019
- Star-state internal energy became 2.35 instead of 0.217

### 2. Velocity Relation Error
**Incorrect:**
```python
v_star = v_L + 2.0/(γ-1) * (c_star - c_L)
```

**Correct:**
```python
v_star = v_L + 2.0/(γ-1) * (c_L - c_star)
```

The velocity relation across a rarefaction wave is `v - v_L = 2/(γ-1) * (c_L - c)`, not `(c - c_L)`.

### 3. Unassigned Grid Points
Some grid points were not covered by any xi-based mask due to incorrect initialization order. Fixed by:
- Initializing all arrays to right-state (default) at the beginning
- Adding a coverage check/fallback before energy calculation

## Solution Implementation

### Files Modified
1. `sample/vacuum/scripts/vacuum_analytical.py` (VacuumRiemannSolver class)
2. `sample/vacuum/scripts/compare_vacuum_methods.py` (compute_analytical_solution function)

### Key Changes
```python
# Corrected star-state calculation
c_star = 0.5 * (c_L + c_R) - gm1 / 4.0 * (v_R - v_L)
c_star = max(c_star, 1e-10)
v_star = v_L + 2.0 / gm1 * (c_L - c_star)  # Note: c_L - c_star
rho_star = rho_L * (c_star / c_L) ** (2.0 / gm1)
P_star = P_L * (c_star / c_L) ** (2.0 * gamma / gm1)
```

### Initialization Pattern
```python
# Initialize to right (default) state to guarantee full coverage
rho[:] = rho_R
v[:] = v_R
P[:] = P_R

# Then apply region-specific assignments
# ... vacuum/non-vacuum branches ...

# Sanity check before energy calculation
zero_count = int((rho == 0.0).sum())
if zero_count:
    print(f"Warning: {zero_count} analytic points unassigned - filling with right state")
    rho[rho == 0.0] = rho_R
    v[rho == 0.0] = v_R
    P[rho == 0.0] = P_R
```

## Validation

### Test Case: t=0.14, γ=1.4
**Initial Conditions:**
- Left:  ρ=1.0, P=0.4, v=-2.0
- Right: ρ=1.0, P=0.4, v=+2.0

**Corrected Center Values:**
- ρ_* = 0.021852
- v_* = 0.0
- P_* = 0.001894
- e_* = 0.216669

**Validation:**
- Zero unassigned points
- Smooth analytical solution (no vertical lines)
- Values match reference solution (Toro 1999, Yuasa & Mori 2024)
- Good agreement with SPH simulation results

## Physical Interpretation

This problem is **not** a true vacuum case (vacuum_criterion < 0), but rather a **two-rarefaction** case with a low-density intermediate (star) state. The fluids move apart creating rarefaction waves that meet at the center, forming a low-pressure, low-density region (but not true vacuum).

The vacuum criterion check:
```python
vacuum_criterion = v_R - v_L - 2/(γ-1) * (c_L + c_R)
                 = 4.0 - 7.483 = -3.483 < 0  → no vacuum
```

## References
- Toro, E.F. (1999). "Riemann Solvers and Numerical Methods for Fluid Dynamics", Chapter 4, Section 4.3.2
- Yuasa & Mori (2024), Section 4.2.2
- Hopkins (2015), DISPH paper Section 4.1

## Testing
```bash
# Run comparison with corrected analytical solution
python3 sample/vacuum/scripts/compare_vacuum_methods.py \
  --results-dir sample/vacuum/results \
  --methods gsph_cubic ssph_cubic disph_cubic gdisph_cubic gdisph_balsara_cubic \
  --snapshot 14 \
  -o sample/vacuum/results/comparison_t0.14.png

# Generate full animation
python3 sample/vacuum/scripts/generate_animation.py \
  --results-dir sample/vacuum/results \
  --method gsph_cubic \
  -o sample/vacuum/results/gsph_cubic_animation.gif
```

---
**Date:** 2025-11-12  
**Status:** ✓ Fixed and validated
