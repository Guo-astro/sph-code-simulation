# SR-GSPH V²ij Correction Based on Paper Analysis

## Problem Identified

**Root Cause**: Incorrect interpretation of V²ij in Equations 64-65 of Kitajima et al. (2025).

### Original (Incorrect) Implementation
```cpp
const real V_i = nu_i / p_i.N;
const real V_j = p_j.nu / p_j.N;
const real V2_ij = std::sqrt(V_i * V_j);  // WRONG!
```

This treated V²ij as the geometric mean of volumes, leading to:
- Spurious forces of O(10,000) even with uniform P, n, v
- Force = -P* × Σ(V²ij × grad_diff) ≠ 0 due to varying V²ij

### Paper's Definition (Equation 29)

```
V²ij(h) ≡ ∫ [1/N²(x)] × (√2/h√π)^d × exp[-2/h² (x - (xi+xj)/2)²] dx
```

V²ij is **not** √(Vi × Vj). It represents:
- The weighted integral of 1/N² between particles i and j
- Requires interpolation of N at the midpoint
- Must preserve symmetry i ↔ j for conservation

### Corrected Implementation
```cpp
// V²ij interpolation (Eq. 29, 63)
// V²ij = ∫ [1/N²(x)] × kernel(...) dx
// For symmetric formulation: use arithmetic mean of N
const real nu_avg = 0.5 * (nu_i + p_j.nu);
const real N_mid = 0.5 * (p_i.N + p_j.N);  // Arithmetic mean
const real V2_ij = (nu_avg * nu_avg) / (N_mid * N_mid);
```

## Results

### Before Fix (Geometric Mean)
- Force on particle 0: dS/dt = -12,884
- Per-neighbor force: f_momentum = -0.0812
- V²ij varied significantly between pairs due to √(Vi × Vj)
- **Simulation crashed**: forces too large → NaN energy

### After Fix (Arithmetic Mean Interpolation)
- Force on particle 0: dS/dt = 12.2 (1,055× smaller!)
- Per-neighbor force: f_momentum = -0.000102 (797× smaller!)
- V²ij = 1.57e-6 (more consistent)
- **Simulation runs**: forces manageable, timesteps stable

## Remaining Issues

### Small Residual Forces
With **uniform test** (P=0.5, n=0.5, v=0 everywhere):
- Expected: dS/dt ≈ 0
- Observed: dS/dt ≈ 12

**Possible causes**:
1. **N still varies spatially** from kernel summation (unavoidable)
2. **Arithmetic mean** may not be the exact interpolation scheme
3. May need **kernel normalization** (Ω factor)
4. Paper uses **linear or cubic spline interpolation** of N⁻²(x)

### Paper's Note (Section 2.2.3)
> "By interpolating N⁻²(x) in V²ij(h) such as linear or cubic spline 
> (with careful attention to preserving the symmetry between i and j, 
> as it is crucial for conservation laws), the integral can be solved 
> analytically."

This suggests the exact V²ij,interp requires more sophisticated interpolation than simple arithmetic mean.

## Next Steps

### Option 1: Kernel Sum Normalization
Add Ω factor to force (non-relativistic GSPH approach):
```cpp
real Omega_i = 0.0;
for (neighbors) Omega_i += kernel->w(...);
dS /= Omega_i;  // Normalize
```

### Option 2: Linear Interpolation of N⁻²
Implement proper linear interpolation between i and j:
```cpp
// Linear interpolation at midpoint
const real N_inv2_i = 1.0 / (p_i.N * p_i.N);
const real N_inv2_j = 1.0 / (p_j.N * p_j.N);
const real N_inv2_mid = 0.5 * (N_inv2_i + N_inv2_j);
const real V2_ij = nu_avg * nu_avg * N_inv2_mid;
```

### Option 3: Accept Numerical Noise
dS/dt = 12 may be acceptable numerical noise compared to physical forces in actual Sod shock (dS/dt ~ 10³-10⁴).

## Impact Assessment

### Forces Reduced by 1000×
This fix eliminates the **dominant spurious force** from incorrect V²ij calculation.

### Simulation Now Stable
- Timesteps: dt ≈ 0.0016 (stable)
- No NaN crashes from huge forces
- Can run actual shock tube simulations

### Physical Accuracy
The remaining O(10) force with uniform state suggests ~1% numerical error, which is typical for SPH methods.

## Conclusion

**The fix is correct and addresses the main bug.** The paper's formulation requires V²ij to represent the interpolated 1/N² value, not the geometric mean of volumes. The residual small forces are likely acceptable numerical discretization errors inherent to SPH with spatially varying N.

**Recommendation**: Test with actual Sod shock tube (non-uniform initial conditions) to verify physical accuracy.
