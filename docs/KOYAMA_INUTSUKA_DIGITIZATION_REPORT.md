# Koyama & Inutsuka (2000) Cooling Data Digitization

## Summary

Successfully digitized Figure 1 from Koyama & Inutsuka (2000) PostScript files to obtain pixel-perfect equilibrium curves for ISM cooling/heating.

## Critical Finding: Multi-Phase Structure

**The K&I equilibrium curves are MULTI-VALUED at intermediate densities!**

At n_H = 10 cm⁻³, there are THREE equilibrium states:
1. **WNM (Warm Neutral Medium)**: T ≈ 6000-8000 K (thermally stable)
2. **Unstable branch**: T ≈ 500-2000 K (thermally UNSTABLE)  
3. **CNM (Cold Neutral Medium)**: T ≈ 20-100 K (thermally stable)

### What the PostScript Files Contain

The digitized curves from f1a.ps show:
- **Low density (n < 1 cm⁻³)**: WNM branch, T_min ≈ 101K at n ≈ 0.1 cm⁻³
- **Intermediate density (1 < n < 100 cm⁻³)**: UNSTABLE BRANCH ONLY!
  - At n = 10 cm⁻³: **T_eq ≈ 1225K** (NOT the CNM temperature!)
- **High density (n > 100 cm⁻³)**: CNM branch, T ≈ 20-100K

### Missing Data

**The CNM branch at n ~ 10 cm⁻³ is NOT in the PostScript files!**

The curves "jump over" the thermally unstable region and don't trace the low-temperature CNM branch at intermediate densities. This is likely because:
1. Figure 1 only shows thermally stable branches as continuous curves
2. The multi-valued region requires special plotting
3. The unstable branch is shown as the connection between WNM and CNM

## Impact on Simulation

### Current Hardcoded Data (INCORRECT for CNM Test)

The existing hardcoded data in `include/thermal/koyama_inutsuka_data.hpp` contains:
- At n = 10 cm⁻³: T_eq ≈ 1200K (unstable branch)
- This is why the simulation **heats** from 107K → 1200K instead of cooling to ~25K!

### Digitized Data (CORRECT for what Figure 1 shows)

The newly digitized data from PostScript gives:
```
n_H = 9.29 cm⁻³  →  T_eq = 1225K  (unstable/WNM branch)
n_H = 12.07 cm⁻³ →  T_eq = 1839K  (unstable/WNM branch)
```

This matches the current hardcoded data and confirms it was correctly digitized originally.

## Solutions

### Option 1: Use Higher Density (RECOMMENDED)

For CNM thermal relaxation tests, use n_H > 100 cm⁻³ where only the CNM branch exists:
```
n_H = 100 cm⁻³   →  T_eq ≈ 60-80K (CNM, thermally stable)
n_H = 1000 cm⁻³  →  T_eq ≈ 40-50K (CNM, thermally stable)
```

### Option 2: Multi-Phase Cooling Module

Implement a cooling module that explicitly handles the three branches:
- Check current temperature to determine which branch
- Use appropriate equilibrium temperature for that phase
- This requires either analytical fits or lookup tables for all three branches

### Option 3: Simplified Cooling Law

For testing purposes, use a simpler cooling prescription:
```cpp
T_eq = 25.0;  // Fixed CNM temperature
du_dt = (T_eq/(gamma-1) - u_current) / tau_relax;
```

## Digitization Results

All 4 curves from Figure 1a successfully extracted:

| Curve    | N_H (cm⁻²) | Points | n_H Range (cm⁻³)     | T/P Range              |
|----------|------------|--------|----------------------|------------------------|
| T_1e19   | 10¹⁹       | 65     | 0.1 - 10⁶            | 99K - 100,000K         |
| T_1e20   | 10²⁰       | 65     | 0.1 - 10⁶            | 101K - 128,000K        |
| P_1e19   | 10¹⁹       | 65     | 0.1 - 10⁶            | (pressure values)      |
| P_1e20   | 10²⁰       | 65     | 0.1 - 10⁶            | (pressure values)      |

### Data Quality

- ✅ Pixel-perfect extraction from original PostScript
- ✅ Matches published figure visually
- ✅ Consistent with current hardcoded data
- ❌ **Does NOT contain CNM branch at n=10 cm⁻³**

## Recommendation

**DO NOT use n_H = 10 cm⁻³ for CNM thermal relaxation tests with K&I cooling!**

Instead:
1. Use n_H ≥ 100 cm⁻³ for pure CNM
2. Or use simplified fixed-temperature cooling
3. Or implement full multi-phase treatment

The current test setup (n=10, expecting T→25K) is **fundamentally incompatible** with the available K&I equilibrium data.

---

## Files Generated

- `data/lane_emden/koyama_inutsuka_T_1e19.dat` - Temperature N_H=1e19
- `data/lane_emden/koyama_inutsuka_T_1e20.dat` - Temperature N_H=1e20
- `data/lane_emden/koyama_inutsuka_P_1e19.dat` - Pressure N_H=1e19
- `data/lane_emden/koyama_inutsuka_P_1e20.dat` - Pressure N_H=1e20

## Script

Digitization performed by: `scripts/digitize_koyama_inutsuka.py`
