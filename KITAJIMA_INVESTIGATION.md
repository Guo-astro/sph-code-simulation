# Investigation: Why Large Density Blip Persists

**Date**: November 25, 2025

---

## Summary

Even with Kitajima's exact particle count (1800+1800 = 3600), density blip remains **54%**.

This is MUCH larger than expected for a working implementation.

---

## Results at Different Resolutions

| Config | Particles | C_cd | Overshoot | Undershoot | Pressure var |
|--------|-----------|------|-----------|------------|--------------|
| Original | 100 (80+20) | 0.2 | 29.8% | 50.2% | 16.3% |
| Fixed | 400 (200+200) | 1.0 | 45.9% | 52.7% | 7.1% |
| **Kitajima** | **3600 (1800+1800)** | **1.0** | **54.9%** | **54.0%** | **5.0%** |

**Key observation**: Density blip INCREASES with resolution, but pressure variation DECREASES.

This suggests:
- ✅ Pressure physics is correct (improves with resolution + C_cd)
- ❌ Density calculation has fundamental issue

---

## What Kitajima Paper Says

### Line 461 (SRGSPH.tex)
> "This **slight** pressure fluctuations at the contact discontinuity..."

- We achieve 5% pressure variation → ✅ **SLIGHT**
- Paper calls it "slight" explicitly

### Line 478 (SRGSPH.tex)
> "the volume-based method exhibits **smaller** overshoot and undershoot than the standard approach"

- Says "smaller", not "small" or "negligible"
- Implies overshoot/undershoot exists in BOTH methods
- Volume-based is just BETTER than standard

**Question**: What magnitude of overshoot/undershoot is normal?

---

## Critical Missing Information

1. **What does Kitajima's Figure 2 actually show?**
   - Need to see the actual plots
   - How large is the density blip in their results?

2. **Is 50% density blip NORMAL for SPH at contacts?**
   - Paper doesn't give quantitative values
   - Only says "smaller" than standard approach

3. **Are we using the correct volume-based formula?**
   - Implementation: N = ν/V_p where V_p = 1/Σ W
   - Is this exactly what Kitajima uses?

4. **What other parameters might matter?**
   - Time integrator?
   - Kernel width?
   - Neighbor count?
   - C_smooth value?

---

## Code Implementation Check

### Your sr_sod.cpp (lines 74-75)
```cpp
const real nu_left = n_left * dx_left;
const real nu_right = n_right * dx_right;
```

This creates **different baryon numbers** per particle (ν_L ≠ ν_R), which matches Kitajima Figure 2.

### Your sr_pre_interaction.cpp (line 225)
```cpp
const real N_new = p_i.nu / Vp;
```

This is the **volume-based approach** (SRGSPH Eq. 243).

**Both implementations look correct!**

---

## Hypothesis

**Possibility 1**: Kitajima's results ALSO show ~50% density blip
- The figures in the paper might not be visible in .tex source
- "Smaller overshoot" could mean 50% vs 70%
- This might be acceptable for SPH

**Possibility 2**: Missing implementation detail
- Some smoothing or limiting we're not applying
- Parameter value not documented in paper
- Difference between paper formulation and actual code

**Possibility 3**: Boundary treatment
- How particles near x=0 are handled
- Discontinuity initialization method
- Ghost particles or special treatment

---

## Next Steps

### 1. Check if there's a working reference implementation

```bash
find . -name "*.py" | xargs grep -l "sod\|Sod"
```

### 2. Try alternative approaches mentioned in paper

From line 478, there's a "standard approach" vs "volume-based approach".

Maybe we need to compare both to see if volume-based is actually better?

### 3. Contact the authors?

The paper is recent (2025). They might have:
- Reference code
- Actual figure values
- Additional implementation notes

### 4. Try exact Figure 1 setup

Figure 1 uses 3200+400 particles (unequal spacing).
Maybe that works better?

---

## Conclusion

**The density blip problem persists regardless of:**
- ✅ C_cd = 1.0 (correct)
- ✅ 2nd order MUSCL (enabled)
- ✅ Variable smoothing length (enabled)
- ✅ High resolution (3600 particles)
- ✅ Volume-based approach (implemented)

**Pressure variation is good** (5%), but **density blip is bad** (54%).

This suggests either:
1. SPH inherently has large density blips at contacts (expected behavior)
2. We're missing a critical implementation detail
3. The paper's figures show better results than we can reproduce

**Need**: Visual comparison with Kitajima's actual Figure 2 results, or reference code.

---

**Investigation**: In progress
**Status**: Unexplained discrepancy between implementation and expected results
