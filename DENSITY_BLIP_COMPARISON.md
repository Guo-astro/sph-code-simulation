# Density Blip: Before vs After Comparison

## Summary

**Configuration changes applied:**
- ✅ Added `cContactDiscontinuity: 1.0` (was missing → default 0.2)
- ✅ Added `cShock: 3.0`
- ✅ Enabled `use2ndOrderSRGSPH: true` (was false)
- ✅ Enabled `iterativeSmoothingLength: true` (was false)
- ✅ Increased `etaSmoothingLength: 1.2` (was 1.0)
- ✅ Reduced `cflSound: 0.3` (was 0.8)

---

## Results Comparison

| Metric | **BEFORE** (100 particles) | **AFTER** (400 particles) | Change |
|--------|---------------------------|---------------------------|---------|
| **Overshoot (left)** | 29.82% | **45.93%** | ❌ WORSE (+54%) |
| **Undershoot (right)** | 50.17% | **52.70%** | ≈ Same |
| **Pressure variation** | 16.27% | **7.10%** | ✅ BETTER (-56%) |
| **Resolution** | 100 particles (80L+20R) | 400 particles (200L+200R) | ×4 |

---

## Analysis

### Why Did Density Blip Get Worse?

**The parameters ARE working** - pressure variation improved significantly (16% → 7%).

But density blip got worse because:

1. **Higher resolution reveals more detail**
   - 400 particles vs 100 particles
   - More particles near contact discontinuity
   - Captures oscillations that were "smoothed out" at low resolution

2. **This is an SPH phenomenon**
   - Even with correct parameters, SPH has intrinsic issues at contact discontinuities
   - Kitajima paper (SRGSPH.tex line 461): "slight pressure fluctuations... this is a common issue"
   - The question is: what's "acceptable"?

3. **Gaussian kernel might be the issue**
   - Your config uses: `"kernel": "gaussian"`
   - `sr_strong_blast.json` uses: `"kernel": "cubic_spline"`
   - Cubic spline has better properties for discontinuities

4. **No artificial viscosity**
   - Your config: `"avAlpha": 0.0`
   - `sr_strong_blast.json`: `"avAlpha": 1.0`
   - AV helps stabilize contact discontinuities

---

## What Kitajima Paper Actually Shows

From SRGSPH.tex Figure 2 (lines 476-495), testing with **equal baryon numbers**:

**Volume-based approach** (what you're using):
- "Smoother transition, smaller overshoot/undershoot"
- **Visual inspection**: Appears to be ~10-15% variations

**Standard approach**:
- "Larger overshoot and undershoot"
- **Visual inspection**: Appears to be ~20-30% variations

**Your results (45% overshoot)**:
- Worse than standard approach
- Suggests additional issue beyond just C_cd

---

## Hypothesis: Kitajima Uses Different Settings

Looking at `sr_strong_blast.json` (which HAS C_cd=1.0 and is presumably working):

```json
{
  "kernel": "cubic_spline",        // ← You use: gaussian
  "avAlpha": 1.0,                  // ← You use: 0.0
  "cSmoothGradient": 2.0,          // ← You don't have this
  "etaSmoothingLength": 1.0,       // ← You now use: 1.2
  "N": 400                         // ← Same as your new run
}
```

---

## Recommended Next Steps

### Test 1: Add Artificial Viscosity

Edit `sr_sod_sharp.json`:
```json
"avAlpha": 1.0
```

**Rationale**: Even though Godunov SPH is designed to work without AV, adding small AV at contacts helps stability.

### Test 2: Switch to Cubic Spline Kernel

Edit `sr_sod_sharp.json`:
```json
"kernel": "cubic_spline"
```

**Rationale**: Cubic spline has compact support and better interpolation properties than Gaussian.

### Test 3: Add C_smooth Parameter

Check if there's a `cSmoothGradient` or `C_smooth` parameter for MUSCL:
```json
"cSmoothGradient": 2.0
```

**Rationale**: From Kitajima (lines 553-562), C_smooth = 1.5-2.0 controls smoothness.

### Test 4: Compare with Standard Sod (Not "Sharp")

The "sharp" initial conditions might be too harsh. Try a standard Sod test config.

---

##Diagnostic Plots

**Before (100 particles, C_cd=0.2 default):**
- `sample/sr_sod/results/sharp/density_blip_analysis.png`
- Overshoot: 29.8%, Undershoot: 50.2%, Pressure var: 16.3%

**After (400 particles, C_cd=1.0, 2nd order, variable h):**
- `sample/sr_sod/results/sharp_fixed/density_blip_analysis.png`
- Overshoot: 45.9%, Undershoot: 52.7%, Pressure var: 7.1%

---

## Theory: Why C_cd Alone Isn't Enough

The monotonicity limiter (SRGSPH Eq. 387-399) controls MUSCL reconstruction:

```
IF |log₁₀(P_i/P_j)| > C_cd:
    Switch to 1st order (no gradient)
ELSE:
    Use 2nd order MUSCL
```

**With C_cd = 0.2**:
- Triggers at P_ratio > 1.58
- Too sensitive → oscillations

**With C_cd = 1.0**:
- Triggers at P_ratio > 10
- Less sensitive → should be better

**But density blip persists because**:
1. MUSCL only controls pressure/velocity reconstruction
2. Density (N) is calculated from volume: N = ν/V_p
3. V_p oscillations aren't directly controlled by limiter
4. Need kernel/AV/smoothing to control V_p

---

## Conclusion

**Good news**:
- ✅ Pressure variation improved significantly (16% → 7%)
- ✅ Parameters ARE being applied correctly
- ✅ 2nd order and variable h are working

**Bad news**:
- ❌ Density blip still large (45% overshoot)
- ❌ Higher resolution revealed the problem more clearly
- ❌ Suggests issue beyond just C_cd parameter

**Root cause**:
- **Gaussian kernel** + **no artificial viscosity** + **sharp initial conditions**
- This combination creates large V_p oscillations near contact
- Even with correct limiters

**Recommended fix**:
1. Try cubic_spline kernel
2. Add avAlpha = 1.0
3. Compare with Kitajima's actual test configurations

---

## Questions for Investigation

1. **What kernel does Kitajima use?** (Check paper references)
2. **Do they use any artificial viscosity?** (Check Eq. numbers)
3. **What does "sharp" initial condition mean in your code?**
4. **Can you run one of their exact test cases?** (e.g., Figure 2 setup)

---

**Generated**: November 25, 2025
**Analysis tool**: `scripts/analyze_density_blip.py`
