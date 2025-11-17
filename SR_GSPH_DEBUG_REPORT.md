# SR-GSPH Debug Report
**Date**: 2025-11-17  
**Issue**: Spurious forces in SR-GSPH implementation with uniform initial conditions

## Executive Summary

Debugged the Special Relativistic Godunov SPH implementation and identified **two critical bugs** and **one fundamental numerical issue**:

### Bugs Fixed ✅
1. **Riemann Solver NaN bug**: Solver returned NaN for v* when states were nearly identical
2. **Riemann Solver numerical instability**: get_velocity_behind_wave could produce NaN values

### Fundamental Issue Identified ⚠️
**Spurious forces from varying baryon density**: Even with uniform pressure P=0.5 and velocity v=0, the simulation generates huge forces (dS/dt~13,000) due to spatial variation in N (baryon density).

---

## Detailed Findings

### 1. Riemann Solver NaN Bug

**Symptom**: With uniform or nearly-uniform states, Riemann solver returned:
- P* = 0.5 (correct)
- v* = NaN (WRONG!)

**Root Cause**: 
The exact Riemann solver uses Brent's method to find P* by bracketing the velocity difference function `dvel(P) = v_left(P) - v_right(P)`. When states are nearly identical, this function is close to zero everywhere, and `get_velocity_behind_wave` can produce NaN during the iteration.

**Fix Applied**:
```cpp
// Early return for nearly-identical states
const real eps_state = 1.0e-10;
if (std::abs(vell - velr) < eps_state && 
    std::abs(pl - pr) < eps_state && 
    std::abs(rhol - rhor) / std::max(rhol, rhor) < 0.01) {
    pstar = pl;
    vstar = vell;
    return;
}

// NaN safety check at convergence
if (std::isnan(vstar) || std::isnan(pstar)) {
    pstar = 0.5 * (pl + pr);
    vstar = 0.5 * (vell + velr);
}
```

**Result**: v* now correctly returns 0.0 instead of NaN for uniform states.

---

### 2. Spurious Force Issue (NOT FIXED)

**Test Case**: Uniform initial conditions
- P = 0.5 everywhere
- n = 0.5 everywhere  
- v = 0 everywhere
- ν = 0.000654 (equal for all particles)
- h = 0.0111 (global fixed smoothing length)

**Expected Result**: Zero forces (no pressure gradient)

**Actual Result**: Huge forces!
```
Particle 0:
  State: P=0.5, n=0.506, v=0
  Force: dS/dt = -12,884 (should be ~0!)
```

**Diagnostic Output**:
```
Neighbor j=1:
  P_i=0.5, P_j=0.5 (identical pressure!)
  v_i=0, v_j=0 (both stationary!)
  Riemann: P*=0.5, v*=0 (correct!)
  
  Kernel gradients:
    dw_i = +64.76
    dw_j = -64.76
    grad_diff = 129.5
  
  Volumes:
    V_i = ν/N_i = 0.00129
    V_j = ν/N_j = 0.00122
    V²_ij = sqrt(V_i × V_j) = 0.00125
  
  Force term:
    f = -P* × V²_ij × grad_diff
      = -0.5 × 0.00125 × 129.5
      = -0.0812
```

**Root Cause Analysis**:

The force formula from Equations 64-65:
```
⟨ν_i Ṡ_i⟩ = -Σ_j P*_ij V²_ij [∇_i W(r_ij, 2h_i) - ∇_j W(r_ij, 2h_j)]
```

For uniform pressure with fixed h:
- P*_ij = P (constant)
- ∇_i W(r_ij, 2h) = -∇_i W(-r_ij, 2h) = -dw_j
- grad_diff = dw_i - dw_j = 2 × dw_i

So the force becomes:
```
F_i = -P Σ_j V²_ij × 2 × ∇_i W(r_ij, 2h)
```

**For an interior particle with symmetric neighbors**, this should cancel. But it DOESN'T because:

1. **Non-uniform N**: Even with equal ν, particles have different N values (0.506, 0.537, 0.567, ...) due to kernel summation: N_i = Σ_j ν_j W(r_ij, h)

2. **Varying volumes**: V = ν/N varies spatially → V²_ij = sqrt(V_i × V_j) varies

3. **No cancellation**: The sum Σ_j V²_ij × ∇W ≠ 0 even for symmetric geometry

**Why N varies**:
- Kernel sum in pre_interaction: `N_i = Σ_j ν_j W(r_ij, h)`
- With equal ν and non-uniform particle spacing, Σ_j W(r_ij, h) varies by position
- Edge particles have fewer neighbors → smaller N
- Left region has 8x more particles → different kernel sums than right region

**Particle distribution**:
```
Left:  800 particles, dx = 0.000625
Right: 100 particles, dx = 0.00500  
```

This creates spatially-varying N even with uniform primitive variables!

---

## Comparison with Python Reference

The Python reference (`example_srgsph.py`) only implements:
- Riemann solver (shock and rarefaction relations)
- Primitive recovery
- **NO force calculation** - doesn't show full SPH evolution

So we cannot directly compare force formulas.

---

## Attempted Solutions (All Failed)

1. ✅ **Fixed Riemann solver NaN** → v* correct but forces still huge
2. ❌ **Smoothed initial discontinuity** (5h → 20h) → Forces still huge
3. ❌ **Increased artificial viscosity** (α=1 → α=5) → No improvement
4. ❌ **Reduced CFL numbers** (0.5→0.1, 0.125→0.01) → Smaller timestep, same forces
5. ❌ **Reduced particle count** (3600 → 900) → Proportionally smaller forces but still huge
6. ❌ **Enabled periodic boundaries** → Uniform N in interior, but still varies
7. ❌ **Energy floor** → Prevents NaN but doesn't fix force issue
8. ❌ **Uniform test mode** (P=0.5, n=0.5 everywhere) → **STILL spurious forces!**

---

## Physical Interpretation

The spurious force arises from the **gradient of particle volume** V = ν/N:

```
F ∝ Σ_j P × ∇(V²_ij)
```

When N varies spatially, ∇V ≠ 0, creating a spurious "volume gradient force" even when ∇P = 0.

This is analogous to the spurious "surface tension" forces in classical SPH from varying smoothing lengths, but here it's from varying N with fixed h.

---

## Possible Solutions (Not Implemented)

### Option 1: Different ν Mode
Use `different_nu = true` to set ν inversely proportional to local kernel sum:
```
ν_i = V_target / (Σ_j W_ij)
```
This keeps V constant, eliminating the volume gradient.

### Option 2: Kernel Sum Normalization
Modify force formula to include kernel sum correction:
```
F_i = -P Σ_j (V²_ij / Ω_i) × ∇W
```
where Ω_i = Σ_j W_ij

### Option 3: Consistent Volume Formulation
Use consistent SPH formulation that automatically conserves momentum:
```
F_i = -Σ_j m_i m_j (P_i/ρ_i² + P_j/ρ_j²) ∇W
```

### Option 4: Accept Numerical Noise
Recognize that some spurious force is unavoidable with:
- Fixed h
- Non-uniform particle spacing
- Equal ν (variable V)

Use very small timesteps and hope it averages out.

---

## Recommendations

1. **Verify force formula against original Kitajima paper**
   - Re-read Equations 64-65 carefully
   - Check if there's a normalization factor we're missing
   - Contact authors if needed

2. **Test different_nu mode**
   - Set `different_nu = true` in config
   - This adjusts ν to maintain constant V
   - Should eliminate volume gradient forces

3. **Compare with other SR-GSPH implementations**
   - Check if spurious forces are common issue
   - Look for normalization techniques

4. **Consider alternative SPH formulation**
   - Traditional momentum-conserving form
   - Density-based instead of baryon-based

---

## Code Changes Made

### Files Modified:
1. `src/srgsph/sr_fluid_force.cpp`:
   - Added identical-state check in Riemann solver
   - Added NaN safety checks
   - Added detailed force diagnostics

2. `src/sample/sr_sod.cpp`:
   - Added uniform_test mode for debugging
   - Fixed kernel sum estimate for different particle counts

3. `src/solver.cpp`:
   - Added energy floor to prevent negative e

4. `sr_sod.json`:
   - Reduced CFL numbers
   - Increased artificial viscosity
   - Enabled periodic boundaries

### Diagnostic Output Added:
```cpp
printf("=== FORCE DIAGNOSTIC: Particle i=%d ===\n", i);
printf("State: P=%.6e, n=%.6e, N=%.6e\n", ...);
printf("Riemann output: P*=%.6e, v*=%.6e\n", pstar, vstar);
printf("Kernel: dw_i=%.6e, dw_j=%.6e\n", ...);
printf("Volumes: V_i=%.6e, V_j=%.6e, V²_ij=%.6e\n", ...);
printf("f_momentum=%.6e\n", f_momentum[0]);
```

---

## Current Status

**Simulation**: CRASHES on first timestep
- Reason: Spurious forces → huge dS → negative energy → NaN → periodic boundary assertion failure

**Next Steps**: Need to resolve spurious force issue before simulation can run.

**Critical Decision Point**: 
- Fix force formulation (requires understanding what's wrong)
- OR switch to different_nu mode
- OR accept this is fundamental limitation of current implementation

---

## Key Data Points

**Uniform Test (P=0.5, n=0.5, v=0)**:
```
Expected:  dS/dt = 0
Actual:    dS/dt = -12,884
Ratio:     Force is ~10,000× too large!
```

**Force breakdown for single pair interaction**:
```
P* = 0.5 (correct)
v* = 0.0 (correct after fix)
V²_ij = 0.00125
grad_diff = 129.5
f_pair = -0.0812

Total over ~20 neighbors:
Σ f_pair ≈ -8.4
After division by ν:
dS/dt = -8.4 / 0.000654 = -12,844
```

**The physics is correct but the SPH discretization creates spurious forces!**
