# Koyama & Inutsuka (2000) Figure Analysis - CNM Identification

## User Question
"Which figure is for CNM?"

## Answer

**Figure 1 shows the thermal equilibrium curves that include ALL phases:**
- **WNM (Warm Neutral Medium)**: n ~ 0.1-1 cm⁻³, T ~ 5000-8000 K (stable)
- **Unstable branch**: n ~ 1-30 cm⁻³, T ~ 500-2000 K (thermally unstable)
- **CNM (Cold Neutral Medium)**: n > 30-100 cm⁻³, T ~ 20-100 K (stable)

### Critical Finding

At **n_H = 10 cm⁻³**, Figure 1 shows:
- **T_eq ≈ 1000-2000 K** (the thermally UNSTABLE branch)
- NOT the CNM equilibrium temperature (~25K)

### Model C10 from Table 1

**Initial conditions:**
- n_H = 10 cm⁻³
- T_init = 107 K  
- N_H = 10²⁰ cm⁻²

**This is a SHOCK PROPAGATION simulation, NOT a thermal relaxation test!**

The gas starts at n=10, T=107K which is:
- Above the equilibrium curve (not in thermal equilibrium)
- In the CNM density range but NOT at equilibrium temperature
- Will be compressed and heated by the shock (Mach 10)

## Why the Confusion?

The current cooling implementation assumes:
1. Single-valued T_eq(n) relationship
2. Gas relaxes from current T to T_eq at given density

But K&I (2000) equilibrium curve is:
1. **Multi-valued** (S-shaped curve with 3 branches)
2. At n=10 cm⁻³, there are potentially 3 equilibrium temperatures:
   - WNM: T ~ 6000K (if extended)
   - Unstable: T ~ 1200K (shown in Figure 1)
   - CNM: T ~ 25K (exists but not shown at this density in the data)

## What the Hardcoded Data Shows

The PostScript digitization from Figure 1a gives:
- At n=10 cm⁻³: **T_eq ≈ 1200 K**
- This is the **thermally unstable branch**
- Marked with **dashed lines (LT2)** in PostScript (unstable)
- **Solid lines (LT1)** show stable WNM and CNM branches

## Correct Interpretation for CNM Tests

For a **true CNM thermal relaxation test**, you need:

### Option 1: Use higher density where CNM is the only equilibrium
- n_H > 100 cm⁻³
- At n ~ 100-1000 cm⁻³, T_eq ~ 10-50 K (stable CNM)
- No multi-valued ambiguity

### Option 2: Implement multi-phase cooling
- Track which branch the gas is on
- CNM branch: Use T_eq ~ (25K) * (n/10)^(-0.1) for n ~ 10-100 cm⁻³
- This requires modifying the cooling module to handle S-curve properly

### Option 3: Use simplified cooling law
- Prescribe T_eq = 25K at n=10 cm⁻³ directly
- Ignore the equilibrium curve multi-valuedness
- Cooling rate: du/dt = (u_eq - u) / τ_cool with τ_cool ~ 0.1 code time

## Figure 1 PostScript Structure

```
LT1 (solid lines)  = Thermally STABLE equilibrium
                     - WNM branch at low n
                     - CNM branch at high n

LT2 (dashed lines) = Thermally UNSTABLE equilibrium  
                     - Intermediate densities
                     - Where Figure 1 shows T~1200K at n=10
```

## Recommendation

The current test setup (n=10, T_init=107K, expecting cooling to 25K) is **fundamentally incompatible** with the K&I equilibrium data because:

1. The equilibrium curve at n=10 shows T_eq ~ 1200K (unstable branch)
2. The CNM equilibrium at ~25K exists at **much higher densities** (n > 100 cm⁻³)
3. Model C10 is a **shock simulation**, not thermal relaxation

**To make a CNM thermal relaxation test work:**
- Change to n_H = 100-1000 cm⁻³ where CNM is the only equilibrium
- Or implement multi-phase aware cooling
- Or use prescribed T_eq=25K (bypassing the K&I data)

## Files Examined

1. `/Users/guo/Downloads/sphcode/docs/papers/cooling-heating/ms.tex`
   - Table 1: Model definitions
   - Section 3.2: Model C10 description (shock into CNM)
   - Figure 1 caption: Equilibrium curves explanation

2. `/Users/guo/Downloads/sphcode/docs/papers/cooling-heating/f1a.ps`
   - Temperature and Pressure equilibrium curves
   - LT1 (solid) = stable, LT2 (dashed) = unstable

3. Current hardcoded data: `include/thermal/koyama_inutsuka_data.hpp`
   - Digitized from Figure 1
   - Shows T_eq ~ 1200K at n=10 (unstable branch)
   - Correctly represents what Figure 1 shows!

## Conclusion

**The hardcoded data is CORRECT for what K&I Figure 1 shows.**

**The test design is WRONG - it expects CNM behavior at n=10 where only the unstable branch exists in the equilibrium curve.**

To fix: Either change density to n>100 cm⁻³, or implement CNM-aware cooling that doesn't rely on the single-valued equilibrium curve.
