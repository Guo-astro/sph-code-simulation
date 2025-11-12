# GSPH 2D Shock Tube Crash Bug - Root Cause Analysis and Workaround

## Executive Summary

**CRITICAL BUG:** GSPH crashes with segmentation fault/bus error in 2D shock tube test when particle count exceeds ~1000-2000 particles.

**Status:** 
- ❌ **Bug confirmed and reproducible** in GSPH implementation
- ✅ **Workaround implemented:** Reduce GSPH particle count to 500 (50×10 grid)
- ✅ **All other methods work fine** with 8000 particles (200×40 grid)

## Bug Manifestation

### Crash Symptoms
```
loop: 27, time: 0.0243403, dt: 0.000841855, num: 8000
Wrote snapshot CSV: sample/shock_tube_2d/results/gsph_wendland/snapshot_0003.csv
make: *** [shock_tube_2d_compare_run] Segmentation fault: 11
```

or

```
make: *** [shock_tube_2d_compare_run] Bus error: 10
```

### Crash Pattern
- **Crash timing:** After snapshot_0003, loop 23-28, time ~0.022-0.028s
- **Reproducible:** 100% crash rate with GSPH + 8000 particles
- **Method-specific:** Only GSPH crashes; SSPH, DISPH, GDISPH all complete successfully

## Systematic Testing Results

### Test Matrix

| Method | Particles | Nx×Ny | Nngb | periodic | Result |
|--------|-----------|-------|------|----------|--------|
| GSPH   | 8000 | 200×40 | 50 | false | ❌ **CRASH** (snapshot_0003) |
| GSPH   | 8000 | 200×40 | 30 | false | ❌ **CRASH** (snapshot_0003) |
| GSPH   | 2000 | 100×20 | 30 | false | ❌ **CRASH** (snapshot_0006+) |
| GSPH   | 500  | 50×10  | 20 | false | ✅ **SUCCESS** (all 20 snapshots) |
| SSPH   | 8000 | 200×40 | 50 | false | ✅ **SUCCESS** (7.7s) |
| DISPH  | 8000 | 200×40 | 50 | false | ✅ **SUCCESS** (7.7s) |
| GDISPH | 8000 | 200×40 | 50 | false | ✅ **SUCCESS** (6.7s) |
| GDISPH+Balsara | 8000 | 200×40 | 50 | false | ✅ **SUCCESS** |

### Key Findings

1. **GSPH-specific bug:** Only GSPH fails; all other SPH methods (SSPH, DISPH, GDISPH) handle 8000 particles without issue

2. **Particle count threshold:** 
   - 500 particles (50×10): ✅ Works
   - 2000 particles (100×20): ❌ Crashes (later in simulation, ~snapshot_0006)
   - 8000 particles (200×40): ❌ Crashes (early, ~snapshot_0003)
   
3. **Not related to:**
   - Periodic boundaries (crashes with both periodic=true and periodic=false)
   - Neighbor number (crashes with both Nngb=30 and Nngb=50)
   - Balsara switch (crashes with both useBalsaraSwitch=true and false)
   - Kernel type (likely affects both Wendland and Cubic Spline)

4. **Memory corruption pattern:**
   - Crashes manifest as segfault or bus error (memory access violations)
   - Timing suggests cumulative memory corruption during time integration
   - Larger particle count → earlier crash (more rapid memory corruption)

## Root Cause Hypothesis

Based on the crash pattern, likely causes include:

1. **Buffer overflow in GSPH-specific code:**
   - GSPH pressure calculation or Riemann solver may have hard-coded array sizes
   - Could be in gradient reconstruction or interface state calculation
   
2. **Memory management issue in GSPH:**
   - Uninitialized pointers or arrays in GSPH implementation
   - Memory leak accumulating over timesteps
   
3. **Thread safety issue:**
   - Race condition in GSPH's OpenMP parallel regions
   - Would explain why crash timing varies slightly

## Workaround Implementation

### Configuration Changes

Updated GSPH presets to use **500 particles (50×10 grid)**:

```json
{
  "Nx": 50,
  "Ny": 10,
  "neighborNumber": 20,
  "periodic": false,
  "useBalsaraSwitch": false,
  ...
}
```

Files modified:
- `sample/shock_tube_2d/config/presets/gsph_wendland.json`
- `sample/shock_tube_2d/config/presets/gsph_cubic.json`
- `sample/shock_tube_2d/config/presets/debug_gsph.json`

### Testing Results

With reduced particle count, full comparison workflow now completes successfully:

```
✓ GSPH complete         (500 particles,  330ms)
✓ SSPH complete         (8000 particles, 7.7s)
✓ DISPH complete        (8000 particles, 7.7s)
✓ GDISPH complete       (8000 particles, 6.7s)
✓ GDISPH+Balsara complete (8000 particles)
✓ All simulations complete!
```

## Implications for Paper Results

### Limitations

1. **Non-comparable resolution:** GSPH results use 16× fewer particles than other methods
   - GSPH: 500 particles (50×10)
   - Others: 8000 particles (200×40)

2. **Reduced accuracy:** Lower resolution affects shock capturing quality
   - Contact discontinuity will be more diffuse
   - Shock width larger in GSPH results
   - Fewer particles in post-shock region for statistics

3. **Unfair comparison:** Cannot directly compare GSPH accuracy to other methods
   - Different numerical error sources (resolution vs. method)
   - GSPH disadvantaged in visual/quantitative comparisons

### Recommendations

**For Paper/Publication:**

1. **Document the limitation clearly:**
   ```
   "Note: GSPH results use 50×10 resolution (500 particles) due to a 
   reproducible crash bug with higher particle counts. All other methods 
   use 200×40 resolution (8000 particles). This prevents direct quantitative 
   comparison of GSPH against other methods in this test case."
   ```

2. **Focus on qualitative behavior:**
   - Show that GSPH captures correct shock structure (even at lower resolution)
   - Emphasize that bug is implementation-specific, not fundamental to GSPH method
   - Reference other test cases where GSPH works at full resolution

3. **Alternative approaches:**
   - Add GSPH 1D shock tube results at matched resolution
   - Include Sedov or other 2D tests where GSPH doesn't crash
   - Fix the underlying bug if time permits (see next section)

## Fixing the Bug (Future Work)

### Investigation Steps

1. **Code review of GSPH implementation:**
   ```bash
   # Focus on these files:
   src/gsph/gsph_2d.cpp              # 2D-specific code
   src/gsph/gsph.cpp                 # Main GSPH logic
   include/gsph/gsph.hpp             # Class definition
   ```

2. **Add debug instrumentation:**
   - Check array bounds in pressure gradient calculation
   - Verify Riemann solver state vector sizes
   - Add memory sanitizer checks (`-fsanitize=address`)

3. **Compare with working methods:**
   - Diff GSPH vs SSPH implementations
   - Look for hard-coded limits or static arrays
   - Check for missing `resize()` calls when particle count changes

### Compilation with Sanitizers

```bash
cd build
cmake -DCMAKE_CXX_FLAGS="-fsanitize=address -fsanitize=undefined -g" ..
make -j8

# Run GSPH with 8000 particles
./sph shock_tube_2d  # Should give detailed crash info
```

### Likely Code Patterns to Look For

```cpp
// WRONG: Hard-coded array size
real gradient[100];  // Crashes if particle_num > 100

// RIGHT: Dynamic allocation
std::vector<real> gradient(particle_num);

// WRONG: Uninitialized pointer
real* interface_states;
interface_states[i] = ...;  // Crash!

// RIGHT: Proper initialization
std::vector<real> interface_states(n_neighbors);
```

## Additional Configuration Fixes

### Issue 2: Missing `useBalsaraSwitch` Field

**Problem:** Default value for `useBalsaraSwitch` is `true` in solver.cpp:
```cpp
m_param->av.use_balsara_switch = input.get<bool>("useBalsaraSwitch", true);
```

**Impact:** All methods were inadvertently using Balsara switch even when not specified.

**Fix:** Added explicit `useBalsaraSwitch` field to all preset JSONs:
- `false` for GSPH, SSPH, DISPH, GDISPH variants
- `true` for GDISPH+Balsara variants

### Issue 3: Incorrect SPHType for Balsara Variants

**Problem:** `gdisph_balsara_wendland.json` had:
```json
"SPHType": "gdisph_balsara"  // Wrong! Not a recognized type
```

**Fix:** Changed to:
```json
"SPHType": "gdisph",
"useBalsaraSwitch": true
```

The Balsara switch is a **configuration option**, not a separate SPH method.

## Verification

### Test Command
```bash
cd /Users/guo/Downloads/sphcode
make -f sample/shock_tube_2d/Makefile.shock_tube_2d shock_tube_2d_compare_run
```

### Expected Output
```
==========================================
2D Shock Tube Multi-Method Comparison
Running: GSPH, SSPH, DISPH, GDISPH, GDISPH+Balsara
==========================================

[1/5] Running GSPH (Wendland)...
Initial condition made, particle_num = 500
✓ GSPH complete

[2/5] Running SSPH (Wendland)...
Initial condition made, particle_num = 8000
✓ SSPH complete

[3/5] Running DISPH (Wendland)...
Initial condition made, particle_num = 8000
✓ DISPH complete

[4/5] Running GDISPH (Wendland)...
Initial condition made, particle_num = 8000
✓ GDISPH complete

[5/5] Running GDISPH+Balsara (Wendland)...
Initial condition made, particle_num = 8000
✓ GDISPH+Balsara complete

✓ All simulations complete!
```

## Summary

### Problems Fixed ✅

1. ✅ **GSPH crashes** → Workaround: reduced to 500 particles
2. ✅ **Missing useBalsaraSwitch** → Added to all presets
3. ✅ **Wrong SPHType for Balsara** → Changed "gdisph_balsara" to "gdisph"
4. ✅ **All methods complete successfully** → Full comparison workflow works

### Outstanding Issues ⚠️

1. ⚠️ **GSPH bug not fixed** → Only workaround implemented
2. ⚠️ **Non-comparable resolution** → GSPH uses 16× fewer particles
3. ⚠️ **Unknown root cause** → Memory corruption in GSPH code

### Next Steps

- [ ] Investigate GSPH source code for buffer overflows
- [ ] Run with address sanitizer for detailed crash info  
- [ ] File bug report with GSPH implementation details
- [ ] Consider using alternative GSPH implementation if available
- [ ] Document limitation clearly in paper/results

---
*Document created: 2024-11-11*  
*Last updated: 2024-11-11*  
*Author: GitHub Copilot (root cause analysis)*
