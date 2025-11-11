# 1D Shock Tube Multi-Method Comparison

This directory contains a complete setup for comparing four SPH formulations on the classic 1D Sod shock tube test:
- **GSPH** (Godunov SPH) - Riemann solver based
- **SSPH** (Standard SPH) - Density-energy formulation
- **DISPH** (Density Independent SPH) - Pressure-energy formulation
- **GDISPH** (Godunov DISPH) - Hybrid combining DISPH pressure-energy with GSPH Riemann solver

## Quick Start

### Prerequisites

⚠️ **This simulation requires `DIM=1` compilation.** The shock tube test is a 1D problem.

```bash
# 1. Configure build for 1D
sed -i '' 's/#define DIM 3/#define DIM 1/' include/defines.hpp

# 2. Rebuild
cd build && make -j8 && cd ..
```

### Run Complete Comparison Study

```bash
# Run all four SPH methods + generate all visualizations
make shock_tube_compare_all
```

This single command will:
1. Run GSPH simulation → `results/gsph/`
2. Run SSPH simulation → `results/ssph/`
3. Run DISPH simulation → `results/disph/`
4. Run GDISPH simulation → `results/gdisph/`
5. Generate comparison plots → `results/comparison/*.png`
6. Generate comparison animation → `results/comparison/comparison_animation.gif`

**Estimated time:** 2-3 minutes total

### View Results

After completion, check:
```bash
ls -lh sample/shock_tube/results/comparison/
```

You'll find:
- `comparison_final.png` - Final state comparison
- `comparison_t*.png` - Snapshots at different times
- `comparison_animation.gif` - Animated evolution

## Manual Step-by-Step Workflow

If you prefer to run each step separately:

### 1. Run All Simulations

```bash
make shock_tube_compare_run
```

This runs all four methods sequentially. Output directories:
- `sample/shock_tube/results/gsph/` - GSPH results
- `sample/shock_tube/results/ssph/` - SSPH results  
- `sample/shock_tube/results/disph/` - DISPH results
- `sample/shock_tube/results/gdisph/` - GDISPH results

### 2. Generate Comparison Plots

```bash
make shock_tube_compare_viz
```

Creates static comparison plots showing all four methods side-by-side.

### 3. Generate Comparison Animation

```bash
make shock_tube_compare_animate
```

Creates an animated GIF showing time evolution of all three methods.

### 4. Clean Up

```bash
make shock_tube_compare_clean
```

Removes all comparison results (gsph/, ssph/, disph/, gdisph/, comparison/)

## Understanding the Results

### Expected Physics

The 1D Sod shock tube is a Riemann problem with:
- **Left state:** ρ=1.0, P=1.0, u=0
- **Right state:** ρ=0.125, P=0.1, u=0
- **Discontinuity:** x=0.5
- **Gas:** γ=1.4 (diatomic)

**Expected solution features:**
1. **Shock wave** - Propagates right
2. **Contact discontinuity** - Separates fluids
3. **Rarefaction wave** - Propagates left

### Comparing SPH Methods

**What to look for in the plots:**

1. **Density Profile:**
   - Sharp shock front capture
   - Contact discontinuity preservation
   - Rarefaction smoothness

2. **Velocity Profile:**
   - Characteristic step pattern
   - Plateau between shock and contact

3. **Pressure Profile:**
   - Pressure jump at shock
   - Decay through rarefaction

4. **Energy Profile:**
   - Energy conservation
   - Distribution pattern

**Method Characteristics:**

- **GSPH** - Sharp features due to Riemann solver, uses artificial viscosity
- **SSPH** - Classic density-energy approach, may show some diffusion
- **DISPH** - Pressure-energy formulation, better contact discontinuities
- **GDISPH** - Hybrid combining DISPH pressure-energy with GSPH Riemann solver (no artificial viscosity in Case 1)

## Advanced Usage

### Run Individual Methods

```bash
# GSPH only
cp sample/shock_tube/config/presets/shock_tube_1d_gsph.json sample/shock_tube/shock_tube.json
./build/sph shock_tube

# SSPH only
cp sample/shock_tube/config/presets/shock_tube_1d_ssph.json sample/shock_tube/shock_tube.json
./build/sph shock_tube

# DISPH only
cp sample/shock_tube/config/presets/shock_tube_1d_disph.json sample/shock_tube/shock_tube.json
./build/sph shock_tube

# GDISPH only
cp sample/shock_tube/config/presets/shock_tube_1d_gdisph.json sample/shock_tube/shock_tube.json
./build/sph shock_tube

# DISPH only
cp sample/shock_tube/config/presets/shock_tube_1d_disph.json sample/shock_tube/shock_tube.json
./build/sph shock_tube
```

### Customize Parameters

Edit config files in `config/presets/`:

```json
{
  "outputDirectory": "sample/shock_tube/results/gsph",
  "endTime": 0.2,           // Simulation end time
  "outputTime": 0.01,       // Snapshot interval
  "neighborNumber": 4,      // Neighbors (1D: use 4)
  "gamma": 1.4,             // Adiabatic index
  "kernel": "cubic_spline", // Kernel function
  "N": 50,                  // Number of particles
  "avAlpha": 1.0,           // Artificial viscosity
  "SPHType": "gsph"         // Method: gsph/ssph/disph
}
```

### Visualize Single Method

```bash
# Visualize GSPH results
python3 sample/shock_tube/scripts/visualize_shock_tube.py \
    sample/shock_tube/results/gsph \
    gsph_output

# Visualize SSPH results
python3 sample/shock_tube/scripts/visualize_shock_tube.py \
    sample/shock_tube/results/ssph \
    ssph_output
```

## Folder Structure

```
sample/shock_tube/
├── Makefile.shock_tube                 # Build targets
├── README_COMPARISON.md                # This file
├── README.md                           # Original single-method docs
├── shock_tube.json                     # Active config (auto-generated)
│
├── config/
│   └── presets/
│       ├── shock_tube_1d_gsph.json    # GSPH configuration
│       ├── shock_tube_1d_ssph.json    # SSPH configuration
│       ├── shock_tube_1d_disph.json   # DISPH configuration
│       ├── shock_tube_1d_gdisph.json  # GDISPH configuration
│       ├── shock_tube_1d_sod.json     # Legacy (same as SSPH)
│       └── shock_tube_1d_sod_resume.json
│
├── scripts/
│   ├── visualize_shock_tube.py        # Single-method visualization
│   ├── compare_shock_tube.py          # Multi-method comparison plots
│   └── animate_comparison.py          # Multi-method animation
│
└── results/
    ├── gsph/                          # GSPH outputs
    │   ├── snapshot_*.csv
    │   └── energy.dat
    ├── ssph/                          # SSPH outputs
    │   ├── snapshot_*.csv
    │   └── energy.dat
    ├── disph/                         # DISPH outputs
    │   ├── snapshot_*.csv
    │   └── energy.dat
    ├── gdisph/                        # GDISPH outputs
    │   ├── snapshot_*.csv
    │   └── energy.dat
    └── comparison/                    # Comparison outputs
        ├── comparison_final.png
        ├── comparison_t*.png
        └── comparison_animation.gif
```

## Makefile Targets Reference

### Single-Method Targets
```bash
make shock_tube_run        # Run SSPH (default)
make shock_tube_resume     # Resume from checkpoint
make shock_tube_animate    # Visualize single method
make shock_tube_clean      # Clean single method results
```

### Multi-Method Targets
```bash
make shock_tube_compare_run      # Run GSPH + SSPH + DISPH + GDISPH
make shock_tube_compare_viz      # Generate comparison plots
make shock_tube_compare_animate  # Generate comparison animation
make shock_tube_compare_all      # Run all + visualizations
make shock_tube_compare_clean    # Clean all comparison results
```

## Troubleshooting

### "ERROR: Shock tube requires DIM=1 build"

You need to rebuild for 1D:
```bash
sed -i '' 's/#define DIM 3/#define DIM 1/' include/defines.hpp
cd build && make -j8 && cd ..
```

### "No snapshot files found"

Run the simulations first:
```bash
make shock_tube_compare_run
```

### "No common snapshots found across all methods"

This means the three methods produced different numbers of snapshots. Check that:
- All configs have the same `endTime` and `outputTime`
- All simulations completed successfully

### Missing Python Dependencies

```bash
pip install numpy matplotlib pillow tqdm
```

## Restoring 3D Build

After shock tube testing, restore the 3D build for other simulations:

```bash
# Restore DIM=3
sed -i '' 's/#define DIM 1/#define DIM 3/' include/defines.hpp

# Rebuild
cd build && make -j8 && cd ..
```

## Scientific Background

### The Sod Shock Tube Test

The Sod shock tube (Sod, 1978) is a standard benchmark for compressible hydrodynamics codes. It tests:
- Shock capturing ability
- Contact discontinuity handling
- Rarefaction wave resolution
- Overall numerical stability

### SPH Method Comparison

**GSPH (Godunov SPH):**
- Uses Riemann solver at particle interfaces
- Generally sharpest shock capture
- More computationally expensive

**SSPH (Standard SPH):**
- Traditional density-energy formulation
- Widely used, well-understood
- May show some diffusion at discontinuities

**DISPH (Density Independent SPH):**
- Pressure-energy formulation
- Better handling of pressure gradients
- Reduced tensile instability

## References

- Sod, G. A. (1978). "A survey of several finite difference methods for systems of nonlinear hyperbolic conservation laws"
- Monaghan, J. J. (2005). "Smoothed particle hydrodynamics"
- Price, D. J. (2012). "Smoothed particle hydrodynamics and magnetohydrodynamics"

## See Also

- `README.md` - Original single-method shock tube documentation
- `../../docs/` - General SPH code documentation
- `config/presets/` - Configuration file templates
