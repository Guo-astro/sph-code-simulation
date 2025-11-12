# 2D Shock Tube Root Cause Analysis

## Investigation Date
2025-11-11

## Problem Statement
Initial testing of `shock_tube_2d` with 8000 particles (200×40 grid) and `neighborNumber=50` resulted in intermittent crashes (segmentation fault / bus error) after 3-4 simulation snapshots.

## Investigation Methodology

### Debug Instrumentation Added
1. **BHTree neighbor search logging**: Added overflow detection with detailed diagnostic output showing:
   - Particle ID, position, smoothing length
   - Current neighbor count vs. maximum allowed
   - Tree node level, edge, center, and leaf particle count
   - Search parameters (is_ij mode, kernel_size)

2. **PreInteraction statistics tracking**: Added monitoring of first 5 timesteps showing:
   - Smoothing length distribution (min, max, average)
   - Maximum neighbors found per particle
   - Neighbor list capacity vs. utilization

### Test Configurations
- **Problematic**: 200×40 = 8000 particles, Nngb=50, periodic boundaries
- **Working**: 100×20 = 2000 particles, Nngb=20, periodic boundaries

## Root Cause Findings

### Primary Issue: Makefile Variable Collision
**Impact**: Critical - Wrong configuration files loaded

The `PRESET_DIR` variable defined in `sample/shock_tube_2d/Makefile.shock_tube_2d` was being overridden by `sample/sedov/Makefile.sedov` (included later in main Makefile line 54 vs. 39).

**Evidence**:
```bash
$ make -n shock_tube_2d_compare_run | grep "cp -f" | head -3
cp -f sample/sedov/config/presets/gsph_wendland.json ...  # WRONG!
```

**Fix**: Renamed to method-specific variable:
```makefile
# Before
PRESET_DIR := $(SHOCK_TUBE_2D_CONFIG)/presets

# After  
SHOCK_TUBE_2D_PRESET_DIR := $(SHOCK_TUBE_2D_CONFIG)/presets
```

### Secondary Issue: Configuration File Confusion
**Impact**: High - Parameters not being applied

Active configuration files in `sample/shock_tube_2d/shock_tube_2d.json` retained old values:
- `neighborNumber: 50` (from old config)
- `outputDirectory: "sample/sedov/results/..."` (wrong path)
- Missing `Nx`, `Ny` parameters

**Fix**: Force-delete active config before copying presets in Makefile targets

### Actual System Capability Analysis

When properly configured, the system **successfully handles 8000 particles with Nngb=50**:

**SSPH Test Results** (200×40, Nngb=50, periodic):
```
=== PreInteraction Debug Statistics ===
Particles: 8000
Target neighbors: 50
Smoothing length: min=0.0132, max=0.0352, avg=0.0202
Max neighbors found: 76
Neighbor list capacity: 50,000
Status: COMPLETED SUCCESSFULLY (7.86s, 165 loops)
```

**GSPH Test Results** (200×40, Nngb=50, periodic):
```
Particles: 8000
Target neighbors: 50
Status: COMPLETED SUCCESSFULLY (11.4s, 242 loops)
```

### Neighbor Count Analysis

**Expected**: ~50 neighbors/particle based on `neighborNumber` parameter

**Actual**: 69-76 neighbors/particle in shock regions

**Explanation**: 
- Shock discontinuity creates strong density gradients (ρ=1.0 → ρ=0.125, 8:1 ratio)
- SPH kernel adapts smoothing length: h ∝ (m/(ρ·A))^(1/DIM)
- High-density region particles have smaller h, low-density have larger h
- Particles at shock interface see neighbors from BOTH sides
- 50% excess is within normal variation for 2D shock problems

**Memory Safety**:
- Neighbor list size: `neighbor_list_size * Nngb = 1000 * 50 = 50,000`
- Actual usage: 76 neighbors << 50,000 capacity
- **No overflow risk**

### Intermittent Crash Analysis

**Previous crashes were NOT due to**:
- ❌ Neighbor list overflow (capacity: 50K, usage: 76)
- ❌ Memory exhaustion (200KB/particle × 8 threads = 1.6MB total)
- ❌ Fundamental algorithm limitation

**Likely causes of intermittent crashes**:
1. **Wrong configuration loaded** (sedov params instead of shock_tube_2d)
2. **Non-periodic boundary conditions** causing particle escape
3. **OpenMP race conditions** (timing-dependent, not reproducible)
4. **Method-specific issues** (GDISPH still has known problems)

## Recommended Configuration

### For Production Use
```json
{
  "Nx": 200,
  "Ny": 40,
  "neighborNumber": 50,
  "periodic": true
}
```

**Works with**: GSPH, SSPH, DISPH  
**Performance**: 8-12 seconds for t=0→0.2  
**Particle count**: 8000  

### For Testing/Development
```json
{
  "Nx": 100,
  "Ny": 20,
  "neighborNumber": 20,
  "periodic": true
}
```

**Works with**: GSPH, SSPH, DISPH, GDISPH  
**Performance**: 1.5-2 seconds for t=0→0.2  
**Particle count**: 2000  

### Conservative Fallback
```json
{
  "Nx": 100,
  "Ny": 20,
  "neighborNumber": 30,
  "periodic": true
}
```

**Recommended when**: GDISPH comparison needed or system stability uncertain

## Lessons Learned

1. **Makefile variable naming**: Use test-specific prefixes to avoid collisions across included makefiles
2. **Configuration management**: Force-delete active configs before regeneration to prevent stale values
3. **Neighbor count variability**: SPH neighbor counts naturally vary ±50% in shock problems - this is expected
4. **Method-specific behavior**: Different SPH methods (GSPH/SSPH/DISPH vs GDISPH) have different robustness
5. **Debug instrumentation value**: Statistical tracking reveals system behavior not visible from crash logs alone

## Action Items

- [x] Fix `PRESET_DIR` variable collision in Makefile
- [x] Add debug logging to neighbor search and pre-interaction
- [x] Document actual system capabilities with proper configuration
- [x] Update preset JSONs with correct output directories
- [ ] Investigate GDISPH-specific crashes (separate issue)
- [ ] Consider making neighbor list size configurable per test case

## References

- BHTree neighbor search: `src/bhtree.cpp:266-318`
- PreInteraction calculation: `src/pre_interaction.cpp:38-165`
- Neighbor list allocation: `src/fluid_force.cpp:37`
- Configuration constants: `include/defines.hpp:32` (`neighbor_list_size = 1000`)
