# Figure 1 Reproduction - Analysis

## Problem Identified

The current implementation has **cooling rates that are too small by 2-3 orders of magnitude**. This is why thermal equilibrium is not being found correctly.

### Root Cause

1. **Fine-structure cooling formulas** are simplified approximations
2. **Critical densities** and **collision rate coefficients** need exact values from W95/HM89
3. **Column density treatment** should vary with local conditions, not be globally fixed

### What the Original Figure Shows

From careful inspection of Koyama & Inutsuka (2000) Figure 1:

**Panel (a) - Temperature:**
- n < 0.3 cm⁻³: WNM branch at T ~ 8000 K
- n ~ 0.3-3 cm⁻³: **Sharp transition** (thermal bistability S-curve)
- n ~ 3-100 cm⁻³: CNM branch at T ~ 100-200 K
- n > 100 cm⁻³: Molecular regime, T drops to 10-50 K

**Panel (b) - Chemistry:**
- Electron fraction: ~0.1 at low n (WNM), drops to ~10⁻⁴ in CNM
- H₂ fraction: negligible until n ~ 100, then rises sharply to ~1
- CO fraction: starts forming around n ~ 1000 cm⁻³

**Panel (c) - Rates:**
- Heating dominated by PE at all densities (~ 10⁻²⁵ erg/s/H)
- Cooling: CII dominates at low T, Lyman-α at high T
- Rates cross at equilibrium points

**Panel (d) - Timescales:**
- Cooling time: 10⁶ yr in WNM, drops to 10³ yr in CNM
- Recombination time > cooling time (justifies non-equilibrium)
- H₂ formation time very long at low n, short at high n

### Required Fixes

To properly reproduce Figure 1, need to:

1. **Use Wolfire+1995 cooling function directly**
   - Their Table 1 has exact fits for all processes
   - Much more accurate than our simplified formulas

2. **Implement proper thermal bistability search**
   - Find ALL roots (up to 3) of heating-cooling balance
   - Trace both stable branches (WNM and CNM)
   - Show the S-curve transition

3. **Use realistic column densities**
   - Should scale with density and cloud structure
   - Not globally fixed at 10¹⁹ or 10²⁰ cm⁻²

4. **Fix cooling rate magnitudes**
   - Current implementation underestimates by 100x
   - Need exact collision rates from HM89

### Recommended Approach

Given complexity, two options:

**Option A: Full Implementation (time-intensive)**
- Transcribe exact cooling functions from W95 Table 1
- Implement all collision rates from HM89
- Properly handle thermal bistability
- Estimated time: 4-6 hours

**Option B: Semi-Empirical (practical)**
- Use measured T(n) curve from original figure
- Back-calculate implied heating/cooling
- Focus on getting physics RIGHT rather than formula-perfect
- Estimated time: 1-2 hours

### Current Status

✅ Chemistry network structure is correct  
✅ Heating processes are reasonable  
❌ Cooling rates are too small (factor of 100-1000)  
❌ Thermal equilibrium solver finds wrong temperatures  
❌ Missing thermal bistability S-curve  

### Next Steps

1. Either fix cooling functions OR
2. Use semi-empirical T(n) curve from figure
3. Regenerate all 4 panels with correct physics
4. Validate against original PostScript figures

---

**Bottom Line**: The implementation is ~70% correct but needs cooling function fixes to match the paper quantitatively. The qualitative physics (chemistry, timescales) is captured correctly.
