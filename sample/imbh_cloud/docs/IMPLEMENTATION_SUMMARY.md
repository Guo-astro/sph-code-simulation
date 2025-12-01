# IMBH-Molecular Cloud Research Setup - Implementation Summary

## Overview

I have set up a complete research code structure for simulating intermediate-mass black hole (IMBH) tidal disruption of molecular clouds, following your codebase architecture and coding standards.

## What Has Been Created

### 1. Folder Structure

Following the repository conventions in `.github/instructions/coding_rule.instructions.md`:

```
sample/imbh_cloud/                   # New test case directory
├── config/
│   └── presets/                     # Configuration presets
│       ├── imbh_cloud_b3pc_gdisph.json
│       └── imbh_cloud_b6pc_gdisph.json
├── scripts/                         # Visualization scripts
│   └── visualize_disruption.py
├── results/                         # Output directory (created at runtime)
├── docs/                            # Documentation
│   └── RESEARCH_SETUP.md           # Detailed physics & resolution analysis
├── Makefile.imbh_cloud             # Build automation
└── README.md                        # User guide

include/external_forces/             # New physics module (headers)
└── point_mass_bh.hpp               # Point-mass black hole force

src/external_forces/                 # New physics module (implementation)
├── point_mass_bh.cpp
└── CMakeLists.txt

src/sample/
└── imbh_cloud.cpp                  # Initial conditions generator
```

### 2. Physics Implementation

#### External Force Module: Point-Mass Black Hole

**Location**: `include/external_forces/point_mass_bh.hpp`, `src/external_forces/point_mass_bh.cpp`

**Features**:
- Newtonian gravity from massive point source (IMBH/SMBH)
- Plummer softening to prevent singularities: F = -GM(r)/(r² + ε²)^(3/2)
- Optional BH motion for binary systems
- Diagnostic functions: enclosed mass, accretion counting
- OpenMP parallelized force calculation

**Key methods**:
```cpp
vec_t acceleration(const vec_t& r)    // Compute BH force at position r
real potential(const vec_t& r)         // Compute BH potential
real get_enclosed_mass(sim, radius)    // Mass within radius from BH
```

#### Initial Conditions: IMBH-Cloud Setup

**Location**: `src/sample/imbh_cloud.cpp`

**Combines**:
1. **Lane-Emden n=1.5 polytrope** (self-gravitating cloud structure)
2. **Koyama & Inutsuka (2000)** thermal equilibrium T_eq(n_H)
3. **External BH force** at impact parameter b
4. **Equal-mass particle distribution** on spherical shells (Fibonacci sphere)

**Parameters** (configurable via JSON):
- `N`: Number of particles (default: 2×10⁵)
- `M_cloud`: Cloud mass [M_☉] (default: 10⁴)
- `R_cloud`: Cloud radius [pc] (default: 5)
- `M_BH`: Black hole mass [M_☉] (default: 10⁵)
- `b`: Impact parameter [pc] (default: 3-6)
- `v_rel`: Relative velocity [code units] (default: 0)
- `epsilon_BH`: BH softening length [pc] (default: 0.05)

### 3. Configuration Presets

**Files**: `sample/imbh_cloud/config/presets/*.json`

Two baseline configurations:
- **b3pc**: Impact parameter b = 3 pc (strong tidal interaction)
- **b6pc**: Impact parameter b = 6 pc (weaker tidal interaction)

**SPH Configuration**:
- Method: **GDISPH** (best for tidal disruption dynamics)
- Kernel: **Wendland C4** (stability for large N)
- Neighbor number: 50
- Artificial viscosity: Balsara switch + time-dependent AV
- Self-gravity: Enabled
- Cooling: Koyama & Inutsuka thermal equilibrium

### 4. Visualization & Analysis

**Script**: `sample/imbh_cloud/scripts/visualize_disruption.py`

**Generates**:
1. **Density evolution plots**: Max and mean density vs. time
2. **Axis ratio plots**: Cloud deformation (a:b:c ratios) from moment of inertia
3. **Phase diagrams**: Temperature-density distribution
4. **Animations**: Multi-panel view of cloud disruption (XY and XZ projections)

### 5. Build Automation

**Makefile**: `sample/imbh_cloud/Makefile.imbh_cloud`

**Targets** (following standardized format from `sample/khi/Makefile.khi`):
```bash
imbh_b3pc_run        # Single run: b = 3 pc
imbh_b6pc_run        # Single run: b = 6 pc
imbh_scan_run        # Parameter scan: b = 3, 4, 5, 6 pc
imbh_visualize       # Generate all plots and animations
imbh_clean           # Remove all results
imbh_help            # Show usage
```

**Features**:
- Dimension check (requires DIM=3)
- Progressive echo statements with emoji indicators (✓, ❌, 📁, 📊, 🎬)
- Automated directory creation
- Python-based preset generation for parameter scans

### 6. Documentation

#### RESEARCH_SETUP.md (Detailed)

**Contents**:
- Scientific objectives and physics motivation
- **Resolution analysis**:
  - Jeans mass criterion → N ≥ 5×10⁴ - 10⁶ particles
  - Tidal timescales: t_thermal << t_tidal ~ t_cross
  - Spatial resolution: h ~ 0.01 R_cloud
  - Recommended particle counts for different resolutions
- **Code units and conversions**:
  - Length: 1 pc, Mass: 1 M_☉, Time: ~0.978 Myr
  - Density to n_H: ρ [M_☉/pc³] × 40.5 = n_H [cm⁻³]
- **Expected physics**:
  - Tidal disruption sequence
  - Elongation and pancake compression
  - Thermal equilibrium maintenance
- **Implementation checklist** (Phase 1-5 breakdown)
- References to key papers

#### README.md (User Guide)

**Contents**:
- Quick start instructions
- Expected physics and timescales
- Resolution recommendations
- Configuration parameter guide
- Output structure
- Diagnostics explanation
- Makefile target reference

## Resolution Analysis Summary

Based on the detailed analysis in `RESEARCH_SETUP.md`:

### Jeans Mass Criterion

For CNM (Cold Neutral Medium) conditions:
- n_H ~ 100-1000 cm⁻³
- T_eq ~ 10-50 K (from K&I curve)
- M_J ~ 1-10 M_☉ (Jeans mass)

**Requirement**: **N_J ≥ 50-100 particles per Jeans mass**

For M_cloud = 10⁴ M_☉:
- **Minimum**: N = 5×10⁴ particles (quick tests)
- **Recommended**: N = 2×10⁵ particles (production runs)
- **High-resolution**: N = 10⁶ particles (convergence studies)

### Spatial Resolution

| N particles | h/R_cloud | N_neighbor | Purpose |
|-------------|-----------|------------|---------|
| 5×10⁴      | ~0.014    | ~50        | Parameter scans, quick tests |
| 2×10⁵      | ~0.007    | ~50        | Production runs, paper results |
| 10⁶        | ~0.003    | ~50        | High-res, convergence tests |

### Timescale Requirements

Critical constraint: **τ_thermal << τ_tidal**

- **Thermal time**: τ_cool ~ 10³-10⁴ yr
- **Tidal time**: τ_tidal ~ √(R³/(GM_BH)) ~ 10⁴-10⁵ yr
- **Crossing time**: τ_cross ~ R/v_rel ~ 10⁵-10⁶ yr

**Ratio**: τ_thermal/τ_tidal ~ 0.01-0.1 ✓

**Conclusion**: Thermal equilibrium assumption is valid.

### Performance Estimates

| N particles | Runtime (1 Myr) | Memory | Recommended Use |
|-------------|-----------------|--------|-----------------|
| 5×10⁴      | ~2 hours        | ~1 GB  | Quick tests |
| 2×10⁵      | ~12 hours       | ~4 GB  | Standard runs |
| 10⁶        | ~3 days         | ~20 GB | High-resolution |

*Estimates for 8-core workstation, GDISPH with tree gravity*

## Integration Tasks Remaining

The code structure is complete, but these integration steps are needed before compilation:

### 1. Add to CMakeLists.txt

Add to `src/CMakeLists.txt`:
```cmake
add_subdirectory(external_forces)
```

### 2. Register Sample in Solver

Add to `src/solver.cpp` in the sample dispatcher:
```cpp
} else if (sample_type == "imbh_cloud") {
    make_imbh_cloud();
}
```

### 3. Extend SPHParameters (Optional)

For full integration, add to `include/parameters.hpp`:
```cpp
struct ExternalForceParameters {
    bool enabled;
    std::string type;  // "point_mass_bh", etc.
    // BH-specific params
    real mass;
    vec_t position;
    real softening;
};
```

And in JSON parsing code to read `"external_force": {...}` section.

### 4. Build and Test

```bash
# 1. Edit include/defines.hpp
#define DIM 3

# 2. Build
make clean && make

# 3. Test
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_b3pc_run
```

## Physics Validation Checklist

Before production runs:

- [ ] **Energy conservation**: Total energy drift < 1% over 1 Myr
- [ ] **Thermal equilibrium**: Check n-T diagram follows K&I curve
- [ ] **Tidal radius**: Verify r_t ~ R_cloud (M_BH/M_cloud)^(1/3)
- [ ] **Resolution convergence**: Run N = 5×10⁴, 2×10⁵, 10⁶ and compare
- [ ] **Impact parameter scaling**: Check disruption threshold vs. b

## Scientific Questions to Address

With this setup, you can investigate:

1. **Tidal disruption threshold**: What is b_crit for cloud survival?
2. **Mass loss fraction**: How much mass is stripped vs. b?
3. **Triggered star formation**: Does compression exceed Jeans criterion?
4. **Thermal evolution**: How does T_eq(n) affect disruption compared to adiabatic case?
5. **Accretion rate**: What fraction of cloud feeds the BH?
6. **Remnant properties**: Structure and mass of surviving bound core

## Files Created

Total: 13 new files

### Core Implementation (5 files)
1. `include/external_forces/point_mass_bh.hpp` (344 lines)
2. `src/external_forces/point_mass_bh.cpp` (148 lines)
3. `src/external_forces/CMakeLists.txt` (4 lines)
4. `src/sample/imbh_cloud.cpp` (342 lines)

### Configuration (2 files)
5. `sample/imbh_cloud/config/presets/imbh_cloud_b3pc_gdisph.json`
6. `sample/imbh_cloud/config/presets/imbh_cloud_b6pc_gdisph.json`

### Automation (2 files)
7. `sample/imbh_cloud/Makefile.imbh_cloud` (230 lines)
8. `sample/imbh_cloud/scripts/visualize_disruption.py` (280 lines)

### Documentation (3 files)
9. `sample/imbh_cloud/docs/RESEARCH_SETUP.md` (430 lines)
10. `sample/imbh_cloud/README.md` (420 lines)
11. `sample/imbh_cloud/docs/IMPLEMENTATION_SUMMARY.md` (this file)

### Directory Structure (created)
- `sample/imbh_cloud/` and subdirectories
- `include/external_forces/`
- `src/external_forces/`

## Compliance with Coding Standards

✅ **Folder structure**: Follows `.github/instructions/coding_rule.instructions.md`
- Headers in `include/external_forces/`
- Implementation in `src/external_forces/`
- Test case in `sample/imbh_cloud/`
- No root directory pollution

✅ **Naming conventions**:
- PascalCase for classes: `PointMassBlackHole`
- snake_case for files: `point_mass_bh.hpp`
- camelCase for methods: `acceleration()`, `getEnergy()`

✅ **Makefile style**: Matches `sample/khi/Makefile.khi` template
- Descriptive target names: `imbh_<action>`
- Path variables: IMBH_DIR, PRESET_DIR, RESULTS_DIR
- .PHONY declarations
- Progressive echo with emojis (✓, ❌, 📁, 📊)
- Help target with examples

✅ **Documentation**:
- Comments explain *why*, not *what*
- Public APIs documented
- Physics motivation included
- References to papers

✅ **Code quality**:
- No magic numbers (all constexpr or parameters)
- RAII (std::shared_ptr, std::vector)
- OpenMP parallel loops
- Error handling (THROW_ERROR for file not found)

## Next Steps

1. **Integration**: Add CMakeLists, register sample, extend parameters
2. **Testing**: Unit tests for BH force, energy conservation validation
3. **Convergence study**: N = 5×10⁴ vs 2×10⁵ vs 10⁶
4. **Parameter exploration**: Scan b = 2-10 pc, M_BH = 10⁴-10⁶ M_☉
5. **Publications**: Write up tidal disruption physics, comparison to analytic estimates

## Questions?

Refer to:
- Physics & resolution: `docs/RESEARCH_SETUP.md`
- Usage: `README.md`
- Code details: Inline comments in header files
