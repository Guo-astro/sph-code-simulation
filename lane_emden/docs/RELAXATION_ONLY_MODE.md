# Lane-Emden Relaxation-Only Mode

## Overview

The relaxation-only mode allows you to run ONLY the relaxation step without running the subsequent simulation. This lets you examine the relaxed particle configuration and all physical variables before deciding whether to proceed with a full simulation.

## Usage

### Quick Start

Run relaxation-only mode with a single command:

```bash
make lane_emden_relax_only
```

This will:
1. Automatically configure the JSON to enable relaxation-only mode
2. Run 100 relaxation steps
3. Save the relaxed particle configuration to `sample/lane_emden/results/00000.dat`
4. Generate comprehensive analysis visualizations

### Output Files

After running `make lane_emden_relax_only`, you'll get:

- **`00000.dat`** - Relaxed particle data (positions, velocities, accelerations, densities, etc.)
- **`relaxed_analysis.png`** - 12-panel comprehensive visualization showing:
  - Spatial distributions (XY projection, 3D view, density map)
  - Physical variables vs radius (density, pressure, energy, smoothing length)
  - Velocities and accelerations (magnitudes and distributions)
- **`relaxed_profiles.png`** - Radially averaged profiles of key quantities
- **`energy.dat`** - Energy evolution during relaxation

### What to Look For

**Good Relaxation:**
- Velocities should be very small (near zero)
- Accelerations should be converging towards hydrostatic balance
- Density profile should be smooth
- Particles should maintain radial distribution

**Issues to Check:**
- Large velocities indicate instability
- Non-converging accelerations suggest relaxation needs more steps
- Irregular density profiles may indicate numerical problems

### Available Make Targets

```bash
make lane_emden              # Run WITHOUT relaxation + visualization
make lane_emden_relax        # Run WITH relaxation (100 steps) + simulation + visualization
make lane_emden_relax_only   # Run ONLY relaxation (no simulation) + analysis
```

## Manual Configuration

To manually enable/disable relaxation-only mode, edit `sample/lane_emden/lane_emden.json`:

```json
{
  "useRelaxation": true,        // Enable relaxation
  "relaxationSteps": 100,       // Number of relaxation steps
  "relaxationOnly": true        // Exit after relaxation (no simulation)
}
```

Then run:
```bash
./build/sph lane_emden
```

## Understanding the Relaxed State

The relaxation process:
1. Calculates SPH forces (pressure, gravity)
2. Subtracts analytical equilibrium forces
3. Zeros all velocities
4. Repeats for N steps

The goal is to reduce accelerations and achieve better force balance before starting the simulation.

### Typical Output Statistics

For n=1.5 Lane-Emden sphere with 5400 particles:

```
Radius range: [~0.01, ~0.99]
Velocity magnitude: ~O(1e-5) to O(1e2)  (should be small after relaxation)
Acceleration magnitude: ~O(1e-3) to O(1e2)  (monitors convergence)
Density range: ~0.008 to ~0.8 (ratio ~100:1 for equal-mass particles)
Total mass: Should be ~1.0 (check for conservation)
```

## Visualization Details

### Comprehensive Analysis (`relaxed_analysis.png`)

**Row 1: Spatial Distributions**
- XY projection colored by radius
- 3D particle distribution
- Density map (log scale)
- Radial histogram

**Row 2: Physical Variables**
- Density profile (log scale)
- Pressure profile (log scale)
- Internal energy profile
- Smoothing length profile

**Row 3: Dynamics**
- Velocity magnitude (should be ~0)
- Acceleration magnitude
- Velocity distribution (log scale)
- Acceleration distribution

### Radial Profiles (`relaxed_profiles.png`)

Shows radially averaged quantities:
- Density profile
- Pressure profile
- Internal energy profile
- Acceleration profile

These should match the analytical Lane-Emden solution.

## Technical Details

### Implementation

The relaxation-only mode is controlled by three parameters:

1. **`useRelaxation`** (bool): Enable relaxation module
2. **`relaxationSteps`** (int): Number of relaxation iterations
3. **`relaxationOnly`** (bool): Exit after relaxation without simulation

### Code Flow

1. `Solver::run()` → `Solver::initialize()`
2. Initialize particles with Lane-Emden profile
3. If `useRelaxation == true`:
   - Run relaxation loop
   - If `relaxationOnly == true`:
     - Output relaxed state
     - Return early from `initialize()`
     - Return early from `run()`
4. Otherwise continue to main simulation

### Files Modified

- `include/solver.hpp` - Added `m_relaxation_only` member
- `src/solver.cpp` - Added relaxation-only logic in `run()` and `initialize()`
- `sample/lane_emden/lane_emden.json` - Added `relaxationOnly` parameter
- `Makefile` - Added `lane_emden_relax_only` target
- `scripts/analyze_relaxed_state.py` - New comprehensive analysis script

## Troubleshooting

**Q: Relaxation doesn't seem to improve the configuration**
- This is expected for 3D Lane-Emden spheres with extreme density gradients
- Relaxation works better for 2D disk geometries
- Consider it as a diagnostic tool rather than a fix

**Q: Velocities are not zero after relaxation**
- Small velocities (O(1e-5)) are acceptable
- Large velocities indicate fundamental instability
- Try increasing `relaxationSteps` or check initial conditions

**Q: Accelerations don't converge**
- For Lane-Emden n=1.5, accelerations may stay at ~193
- This is a known limitation for 3D spheres
- The discrete SPH representation cannot perfectly match the analytical solution

## Next Steps

After examining the relaxed state:

1. If satisfied, set `relaxationOnly: false` and run full simulation
2. If velocities/accelerations are too large, increase `relaxationSteps`
3. If results look unstable, review initial particle placement
4. Use the analysis plots to identify problematic regions

---

Created: 2025-11-10
Last Updated: 2025-11-10
