# Lane-Emden Multi-Method Benchmark Structure

## Overview
This document describes the folder structure and workflow for running 5-method Lane-Emden benchmarks.

## Folder Structure

```
lane_emden/
├── 2d/
│   ├── config/
│   │   ├── active_config.json              # Active config (dynamically generated)
│   │   └── presets/
│   │       ├── benchmark_ssph.json         # SSPH preset
│   │       ├── benchmark_gsph.json         # GSPH preset
│   │       ├── benchmark_disph.json        # DISPH preset
│   │       ├── benchmark_gdisph.json       # GDISPH preset
│   │       └── benchmark_gdisph_balsara.json  # GDISPH+Balsara preset
│   ├── results/
│   │   └── benchmark_n1_5/                 # Benchmark results root
│   │       ├── ssph/                       # SSPH results
│   │       │   ├── snapshot_*.csv
│   │       │   └── energy.dat
│   │       ├── gsph/                       # GSPH results
│   │       ├── disph/                      # DISPH results
│   │       ├── gdisph/                     # GDISPH results
│   │       ├── gdisph_balsara/             # GDISPH+Balsara results
│   │       └── comparison/                 # Comparison plots
│   │           ├── density_comparison.png
│   │           ├── pressure_comparison.png
│   │           ├── animation_comparison.gif
│   │           └── error_analysis.png
│   ├── Makefile.2d                         # 2D workflow automation
│   └── scripts/
│       ├── compare_methods.py              # Multi-method comparison plots
│       └── benchmark_analysis.py           # Error analysis vs analytical
│
└── 3d/
    ├── config/                             # Same structure as 2d/
    ├── results/
    │   └── benchmark_n1_5/
    │       ├── ssph/
    │       ├── gsph/
    │       ├── disph/
    │       ├── gdisph/
    │       ├── gdisph_balsara/
    │       └── comparison/
    ├── Makefile.3d
    └── scripts/
```

## Method Configurations

### 1. SSPH (Standard SPH)
- **Type**: `ssph`
- **Balsara**: `false`
- **Description**: Baseline density-energy formulation

### 2. GSPH (Godunov SPH)
- **Type**: `gsph`
- **2nd Order**: `true`
- **Description**: Riemann solver based method

### 3. DISPH (Density Independent SPH)
- **Type**: `disph`
- **Description**: Pressure-energy formulation

### 4. GDISPH (Godunov + DISPH)
- **Type**: `gdisph`
- **2nd Order**: `true`
- **Balsara**: `false`
- **Description**: Combined Riemann+pressure-energy

### 5. GDISPH+Balsara
- **Type**: `gdisph`
- **2nd Order**: `true`
- **Balsara**: `true`
- **Description**: GDISPH with Balsara viscosity switch

## Makefile Targets

### Benchmark Workflow (2D)
```bash
make -f Makefile.2d benchmark_run        # Run all 5 methods
make -f Makefile.2d benchmark_viz        # Generate comparison plots
make -f Makefile.2d benchmark_analysis   # Error analysis vs analytical
make -f Makefile.2d benchmark_all        # Complete benchmark workflow
make -f Makefile.2d benchmark_clean      # Clean benchmark results
```

### Individual Method Testing (2D)
```bash
make -f Makefile.2d ssph_run             # Run SSPH only
make -f Makefile.2d gsph_run             # Run GSPH only
make -f Makefile.2d disph_run            # Run DISPH only
make -f Makefile.2d gdisph_run           # Run GDISPH only
make -f Makefile.2d gdisph_balsara_run   # Run GDISPH+Balsara only
```

### 3D Benchmarks
```bash
make -f Makefile.3d benchmark_all        # Complete 3D benchmark
```

## Performance Improvements

### Binary Search Optimization
**File**: `src/relaxation/lane_emden_data.cpp`

**Change**: Replaced O(n) linear search with O(log n) binary search using `std::lower_bound`

**Impact**:
- Old: 18 million O(n) searches for 900 particles × 10,000 steps
- New: 18 million O(log n) searches
- **Expected speedup**: ~100-1000× for relaxation phase

**Code**:
```cpp
// Binary search for interpolation interval (O(log n) vs O(n))
auto it = std::lower_bound(m_xi_array.begin(), m_xi_array.end(), xi);
```

## Comparison Metrics

### Quantitative Metrics
1. **Density Error**: L2 norm vs analytical Lane-Emden solution
2. **Pressure Error**: RMS deviation from hydrostatic equilibrium
3. **Energy Conservation**: Total energy drift over simulation time
4. **Computational Cost**: Wall-clock time per method

### Qualitative Metrics
1. **Density Profile**: ρ(r) vs analytical θ^n
2. **Pressure Profile**: P(r) vs analytical equilibrium
3. **Particle Distribution**: Symmetry and clustering
4. **Shock Handling**: Stability at boundaries

## Workflow Steps

### Step 1: Run Benchmark
```bash
cd /Users/guo/Downloads/sphcode/lane_emden/2d
make -f Makefile.2d benchmark_run
```
- Runs 5 methods sequentially
- Progress: [1/5] SSPH, [2/5] GSPH, etc.
- Output: Separate directories per method

### Step 2: Generate Comparisons
```bash
make -f Makefile.2d benchmark_viz
```
- Overlays density, pressure profiles
- Side-by-side snapshots
- Animated comparison

### Step 3: Error Analysis
```bash
make -f Makefile.2d benchmark_analysis
```
- Compares each method to Lane-Emden analytical solution
- Generates error plots and tables

## Output Organization

### Results Directory
```
lane_emden/2d/results/benchmark_n1_5/
├── ssph/snapshot_0000.csv               # t=0 initial condition
├── ssph/snapshot_0050.csv               # t=5 final state
├── ssph/energy.dat
├── comparison/density_t0.png            # 5-method density at t=0
├── comparison/density_t5.png            # 5-method density at t=5
└── comparison/error_summary.txt         # Quantitative comparison
```

## Naming Convention

- **Folders**: `method_name` (ssph, gsph, disph, gdisph, gdisph_balsara)
- **Configs**: `benchmark_<method>.json`
- **Consistent underscore usage** (no dots in folder names)

## Notes

- 2D uses analytical ICs (relaxation disabled due to segfault)
- 3D can use relaxation (100-10000 steps)
- CSV-only output (HDF5 disabled for simplicity)
- All methods use same physics: γ=5/3, n=1.5 polytrope
