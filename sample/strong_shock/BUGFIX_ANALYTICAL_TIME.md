# Bug Fix: Analytical Solution Time Evolution

## Issue

The analytical solution overlay was showing a **static solution** that did not evolve with time. All snapshots displayed the same black curve (exact solution), even though the simulation time was advancing.

**Symptom**: In the density profile plots, the black line (exact solution) remained identical across all 25 snapshots from t=0.000 to t=0.012.

## Root Cause

The metadata parser in `strong_shock_analytical.py` was looking for a key called `'time'` (lowercase):

```python
# BEFORE (BROKEN)
time_str = data['metadata'].get('time', '0.0')
t = float(time_str)
```

However, the CSV output from the SPH simulation uses the key `'Time (code)'`:

```
# Time (code): 1.2013234525623360e-02
```

Since `'time'` was not found, the parser defaulted to `'0.0'` for all snapshots, resulting in:
- **t = 0.0** for all snapshots
- Analytical solution evaluated at initial conditions only
- No wave propagation shown in exact solution

## Solution

Updated the time extraction to try multiple metadata keys:

```python
# AFTER (FIXED)
metadata = data['metadata']
time_str = metadata.get('Time (code)', 
                        metadata.get('time', 
                        metadata.get('Time', '0.0')))
t = float(time_str)
```

This fallback chain ensures compatibility with:
1. Current CSV format: `'Time (code)'`
2. Potential alternative: `'time'`
3. Legacy format: `'Time'`
4. Default: `'0.0'` only if none are found

## Verification

### Before Fix
```
Snapshot 0000: t = 0.000000
Snapshot 0005: t = 0.000000  ← WRONG
Snapshot 0010: t = 0.000000  ← WRONG
Snapshot 0015: t = 0.000000  ← WRONG
Snapshot 0020: t = 0.000000  ← WRONG
Snapshot 0024: t = 0.000000  ← WRONG
```

### After Fix
```
Snapshot 0000: t = 0.000000  ✓
Snapshot 0005: t = 0.002514  ✓
Snapshot 0010: t = 0.005035  ✓
Snapshot 0015: t = 0.007519  ✓
Snapshot 0020: t = 0.010012  ✓
Snapshot 0024: t = 0.012013  ✓
```

### Analytical Solution Evolution Test

Verified the solution now evolves correctly:

```
Density at x=0.0, 0.1, 0.2, 0.3 at different times:
------------------------------------------------------------
t=0.000: ρ = [0.575, 0.575, 0.575, 0.575]
t=0.003: ρ = [0.575, 1.000, 1.000, 1.000]  ← Shock moving right
t=0.006: ρ = [0.575, 0.575, 1.000, 1.000]  ← Shock propagating
t=0.009: ρ = [0.575, 0.575, 5.992, 1.000]  ← High density shock
t=0.012: ρ = [0.575, 0.575, 0.575, 1.000]  ← Continued evolution
```

The density values now change with time, showing:
- **Shock wave** propagating to the right
- **Rarefaction fan** expanding to the left
- **Contact discontinuity** separating fluids

## Impact

All 125 analytical overlay plots (5 methods × 25 snapshots) were regenerated with the corrected time evolution:

```
gsph_cubic:           25 plots regenerated
ssph_cubic:           25 plots regenerated
disph_cubic:          25 plots regenerated
gdisph_cubic:         25 plots regenerated
gdisph_balsara_cubic: 25 plots regenerated
```

File sizes increased from ~400-450 KB to ~415-535 KB due to more complex exact solution curves.

## Files Modified

- **`sample/strong_shock/scripts/strong_shock_analytical.py`**
  - Function: `plot_with_analytical()`
  - Lines: 291-296 (time extraction logic)

## Regeneration

All analytical overlays were regenerated using:

```bash
make strong_shock_compare_viz
```

This automatically runs `process_all_snapshots.py` for all 5 methods, generating 125 corrected plots.

## Testing

Tested on:
- All 5 SPH methods (GSPH, SSPH, DISPH, GDISPH, GDISPH+Balsara)
- All 25 snapshots (t=0.000 to t=0.012)
- Verified time progression is smooth and continuous
- Confirmed shock, rarefaction, and contact structures evolve correctly

## Date

**Fixed**: 2025-11-12  
**Reported by**: User observation (screenshot showing static black line)  
**Severity**: High (invalidated all analytical comparisons)  
**Status**: ✅ Resolved and verified

---

**Lesson Learned**: Always verify metadata key names match between data producer (C++ simulation) and consumer (Python analysis). Use fallback chains for robustness.
