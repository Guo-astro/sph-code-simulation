# SR-GSPH Solver Simplification (2025-11 Update)

The SR-GSPH implementation now follows the Python reference verbatim and always
uses the exact relativistic Riemann solver (with an automatic HLLC fallback).
The previous `riemannSolverSRGSPH` toggle and the iterative solver variant have
been removed to eliminate divergent code paths and hard-to-debug behavior.

## Key Points

- ✅ **Always EXACT:** Every SR-GSPH run now uses the exact solver, matching the
  `srg_sph.py` reference implementation.
- 🗑️ **No config knob:** Remove any `"riemannSolverSRGSPH"` entries from JSON
  presets—they are ignored as of this change.
- 🧠 **Consistent limiter values:** Defaults follow the Python reference
  (`C_shock = 3.0`, `C_cd = 0.2`). Override them only if you know you need to.
- 🧪 **Validation strategy:** Compare against the Python driver or the published
  SR Sod snapshots rather than toggling solvers.

## Migrating Existing Runs

1. Delete the obsolete field from configs:

   ```json
   {
     "SPHType": "srgsph",
     "use2ndOrderSRGSPH": true
   }
   ```

2. Rebuild and rerun using `./build/sph <preset>.json`.
3. (Optional) Regenerate animations via
   `python3 sample/sr_sod/scripts/animate_sr_sod.py sample/sr_sod/results/sharp`.

## Historical Note

The iterative solver prototype served as a debugging aid while porting the
Python formulation. With the C++ force pipeline now mirroring the reference, a
single solver path keeps the results reproducible and the codebase easier to
maintain. This document is retained to explain why the toggle disappeared and
how to update old presets.
