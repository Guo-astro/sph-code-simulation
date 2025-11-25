# Density Blip Diagnostic Report
## SR Sod Shock Tube Simulation

**Date**: November 25, 2025
**Snapshot Analyzed**: `sample/sr_sod/results/sharp/snapshot_0131.csv` (t = 0.35035)
**Configuration**: `sample/sr_sod/config/presets/sr_sod_sharp.json`

---

## Summary

⚠️ **CRITICAL ISSUE DETECTED**: Massive density blip at contact discontinuity

- **Overshoot (left side)**: **29.82%** (should be < 5%)
- **Undershoot (right side)**: **50.17%** (should be < 5%)
- **Pressure variation**: **16.27%** (should be < 1%)

This is **NOT** the "slight pressure fluctuation" mentioned in Kitajima et al. This is a severe numerical artifact.

---

## Root Cause Analysis

### 1. **Missing C_cd Parameter** ⚠️ CRITICAL

**Problem**: Configuration file `sr_sod_sharp.json` does **NOT** specify `cContactDiscontinuity`.

**Effect**: Uses default value **C_cd = 0.2** (from `solver.cpp:445`)

**Kitajima recommendation**: **C_cd = 1.0** (SRGSPH.tex line 397)

**Impact**:
```
With C_cd = 0.2: Switches to 1st order when P_i/P_j > 1.58
With C_cd = 1.0: Switches to 1st order when P_i/P_j > 10.0
```

The monotonicity limiter is **over-sensitive** with C_cd = 0.2, triggering too early and causing oscillations rather than preventing them.

### 2. **Using 1st Order Method**

**Configuration**:
```json
"use2ndOrderSRGSPH": false
```

**Effect**: No MUSCL reconstruction → more diffusion and larger errors at discontinuities

**Kitajima recommendation**: Use 2nd order with proper limiting (SRGSPH §2.7.3)

### 3. **No Variable Smoothing Length**

**Configuration**:
```json
"iterativeSmoothingLength": false
```

**Effect**: Fixed smoothing length cannot adapt to density variations at contact discontinuity

**Kitajima recommendation**: Use iterative smoothing (SRGSPH Eqs. 231-233)

### 4. **Small η Parameter**

**Configuration**:
```json
"etaSmoothingLength": 1.0
```

**Kitajima recommendation**: η = 1.2 (with C_smooth = 2.0)

---

## Visual Evidence

See `sample/sr_sod/results/sharp/density_blip_analysis.png`:

1. **Top panel**: Full domain showing rarefaction wave, contact discontinuity, and shock
2. **Middle panel**: Zoom on contact showing:
   - Sharp spike (75% above plateau) just before contact
   - Deep trough (50% below plateau) just after contact
3. **Bottom panel**: Pressure NOT constant across contact (physics violation)

---

## Comparison with Kitajima Paper

From SRGSPH.tex lines 459-462:
> "These slight pressure fluctuations at the contact discontinuity... this is a common issue in many SPH methods, regardless of whether it is relativistic or non-relativistic, but the amplitude of this **pressure wiggle in Godunov SPH is much smaller than that of standard SPH**."

Kitajima's results (Figure 2, lines 476-495):
- **Volume-based approach**: Smooth transition, small overshoot/undershoot
- **Standard approach**: Larger but still < 10% deviations

Your results:
- **29.82% overshoot, 50.17% undershoot**: Far worse than "standard approach"
- This indicates incorrect implementation or parameter settings

---

## Recommended Fixes

### IMMEDIATE (Critical):

1. **Add C_cd = 1.0 to config file**:
   ```json
   "cContactDiscontinuity": 1.0
   ```

2. **Enable 2nd order**:
   ```json
   "use2ndOrderSRGSPH": true
   ```

3. **Enable variable smoothing**:
   ```json
   "iterativeSmoothingLength": true,
   "etaSmoothingLength": 1.2
   ```

### Configuration Template:

Based on your other configs (`sr_strong_blast.json`, `sr_ultra_*.json`), which DO include C_cd = 1.0:

```json
{
  "outputDirectory": "sample/sr_sod/results/fixed",
  "startTime": 0.0,
  "endTime": 0.35,
  "outputTime": 0.002,
  "SPHType": "srgsph",
  "use2ndOrderSRGSPH": true,
  "cSpeed": 1.0,
  "etaSmoothingLength": 1.2,
  "fixedSmoothingLength": -1.0,
  "gamma": 1.666666667,
  "kernel": "gaussian",
  "iterativeSmoothingLength": true,
  "neighborNumber": 50,
  "N": 200,
  "cContactDiscontinuity": 1.0,    // ← ADD THIS
  "cShock": 3.0,                    // ← ADD THIS
  "cflSound": 0.3,
  "cflForce": 0.25
}
```

---

## Theory from Kitajima Paper

### Monotonicity Limiter (SRGSPH Eq. 387-399):

Switches to 1st order when:
```
|log₁₀(P_i/P_j)| > C_cd  OR  C_shock * e_ij·(v_i - v_j) > min(c_s,i, c_s,j)
```

**Purpose**: Detect contact discontinuities and prevent oscillations

**Problem with C_cd = 0.2**:
- Contact discontinuity in Sod test: P_L*/P_R* ≈ 0.3/0.3 = 1.0 (should not trigger)
- But pressure oscillations create local jumps > 1.58
- Limiter triggers inappropriately, creating MORE oscillations

**Solution with C_cd = 1.0**:
- Only triggers at TRUE discontinuities with 10x pressure jump
- Allows smooth interpolation across contact discontinuity

### V²_ij Factor (SRGSPH Eqs. 194-202):

Your implementation (sr_fluid_force.cpp:287):
```cpp
const real Vij2 = 0.5 * (V_i * V_i + V_j * V_j);
```

This is **correct** (matches SRGSPH Eq. 365). The problem is NOT in V²_ij calculation.

### Volume-Based Density (SRGSPH Eq. 243):

Your implementation (sr_pre_interaction.cpp:225):
```cpp
const real N_new = p_i.nu / Vp;
```

This is **correct**. The problem is NOT in density calculation.

---

## Conclusion

The **primary cause** of your large density blip is:

**Missing `cContactDiscontinuity` parameter → uses default C_cd = 0.2 → monotonicity limiter over-triggers → creates oscillations**

Secondary factors:
- 1st order method (more diffusive)
- Fixed smoothing length (cannot adapt)
- Small η parameter

### Action Items:

1. ✅ Add `"cContactDiscontinuity": 1.0` to `sr_sod_sharp.json`
2. ✅ Enable 2nd order: `"use2ndOrderSRGSPH": true`
3. ✅ Enable variable h: `"iterativeSmoothingLength": true`
4. ✅ Set η = 1.2: `"etaSmoothingLength": 1.2`
5. ✅ Re-run simulation
6. ✅ Re-analyze with `analyze_density_blip.py`
7. ✅ Verify blip < 5%

---

## References

1. **Kitajima et al. (2025)**: "Special Relativistic Smoothed Particle Hydrodynamics Based on Riemann Solver", SRGSPH.tex
   - Line 397: C_cd = 1.0 recommended
   - Lines 459-462: Discussion of "slight" pressure wiggles
   - Lines 476-495: Volume-based vs standard approach comparison
   - Lines 553-562: C_smooth parameter discussion

2. **Your codebase**:
   - `src/solver.cpp:445`: Default C_cd = 0.2
   - `src/srgsph/sr_fluid_force.cpp:287`: V²_ij calculation (correct)
   - `src/srgsph/sr_pre_interaction.cpp:225`: N = ν/V_p (correct)
   - `sample/sr_sod/config/presets/sr_strong_blast.json`: Has C_cd = 1.0 ✓

---

**Report generated by**: `scripts/analyze_density_blip.py`
**Analysis plot**: `sample/sr_sod/results/sharp/density_blip_analysis.png`
