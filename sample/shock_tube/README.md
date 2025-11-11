# Shock Tube Test

1D Sod shock tube test - a classic benchmark for compressible SPH solvers.

## Important: Requires 1D Build

⚠️ **This simulation requires `DIM=1` compilation.** The shock tube test is a 1D problem.

Before running, you must:

1. **Reconfigure build for 1D**:
   ```bash
   # Edit include/defines.hpp and change DIM to 1
   sed -i '' 's/#define DIM 3/#define DIM 1/' include/defines.hpp
   
   # Rebuild
   cd build && make -j8 && cd ..
   ```

2. **Run shock tube**:
   ```bash
   make shock_tube_run
   ```

3. **Visualize results**:
   ```bash
   make shock_tube_animate
   ```
   
   This generates:
   - `shock_tube_sod_final.png` - Final state (4 subplots)
   - `shock_tube_sod_comparison.png` - Time evolution comparison
   - `shock_tube_sod_animation.gif` - Animated evolution

4. **Resume from snapshot**:
   ```bash
   make shock_tube_resume
   make shock_tube_resume CHECKPOINT=sample/shock_tube/results/snapshot_0010.csv
   ```

5. **Restore 3D build** (for lane_emden and other simulations):
   ```bash
   # Edit include/defines.hpp and change DIM back to 3
   sed -i '' 's/#define DIM 1/#define DIM 3/' include/defines.hpp
   
   # Rebuild
   cd build && make -j8 && cd ..
   ```

## Physics

**Classic Riemann Problem:**
- Left state: ρ=1.0, P=1.0
- Right state: ρ=0.125, P=0.1
- Discontinuity at x=0.5
- γ=1.4 (diatomic gas)

**Expected solution:**
- Shock wave propagating right
- Contact discontinuity
- Rarefaction wave propagating left

## Configuration

Presets are located in `config/presets/`:
- `shock_tube_1d_sod.json` - Initial run
- `shock_tube_1d_sod_resume.json` - Resume template

## Workflow (SSOT System)

All physics parameters are stored in snapshot metadata (Single Source of Truth).

**Initial run:**
```bash
make shock_tube_run
```

**Resume from latest:**
```bash
make shock_tube_resume
```

**Resume from specific snapshot:**
```bash
make shock_tube_resume CHECKPOINT=sample/shock_tube/results/snapshot_0010.csv
```

**Visualize:**
```bash
make shock_tube_animate
```

Visualization generates 3 files showing all physics variables:
- **Density** (ρ) - Should show shock, contact discontinuity, rarefaction
- **Velocity** (u) - Should show characteristic step pattern
- **Pressure** (P) - Should show pressure jump at shock
- **Internal Energy** (e) - Should track thermal energy distribution

Resume config only needs:
```json
{
  "checkpoint": {"resumeFile": "path/to/snapshot.csv"},
  "outputDirectory": "sample/shock_tube/results",
  "endTime": 0.4,
  "outputTime": 0.01
}
```

All physics parameters (γ, neighbor count, kernel, etc.) are loaded from snapshot metadata.
