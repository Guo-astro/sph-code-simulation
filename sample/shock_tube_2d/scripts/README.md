# 2D Shock Tube Visualization Scripts

## Overview

This directory contains Python scripts for visualizing and analyzing 2D shock tube (Sod problem) simulation results from the SPH code.

## Scripts

### 1. `shock_tube_2d_analytical.py`
Main script for creating comparison plots between SPH results and the exact Riemann solution.

**Features:**
- Computes exact analytical solution using Riemann solver
- Bins 2D particle data in x-direction (averages over y)
- Creates 4-panel comparison plots: density, velocity, pressure, internal energy
- Overlays SPH results on analytical solution

**Usage:**
```bash
python3 shock_tube_2d_analytical.py <snapshot.csv> -o <output.png>
python3 shock_tube_2d_analytical.py snapshot_0010.csv -o comparison.png --no-show
```

### 2. `process_all_snapshots.py`
Batch processing script that creates analytical comparison plots for all snapshots in a results directory.

**Usage:**
```bash
python3 process_all_snapshots.py <results_dir>
python3 process_all_snapshots.py sample/shock_tube_2d/results/gsph_wendland
```

**Output:**
- Creates `plots/` subdirectory in results directory
- Generates `snapshot_XXXX_comparison.png` for each snapshot

### 3. `plot_2d_spatial.py`
Creates 2D spatial visualizations showing particle distributions in the x-y plane.

**Features:**
- 4-panel layout: density, pressure, x-velocity, internal energy
- Includes velocity vector field overlay
- Shows initial discontinuity location (x=0.5)
- Proper aspect ratio for shock tube geometry (1.0 × 0.2)

**Usage:**
```bash
python3 plot_2d_spatial.py <snapshot.csv> -o <output.png>
python3 plot_2d_spatial.py snapshot_0010.csv -o spatial_2d.png --no-show
```

### 4. `generate_all_animations.py`
Creates animated GIFs showing time evolution of SPH simulations with analytical overlay.

**Features:**
- Generates animations for all SPH methods found in results directory
- Each frame shows 4-panel comparison with analytical solution
- Configurable frame rate

**Usage:**
```bash
python3 generate_all_animations.py --results-dir sample/shock_tube_2d/results --fps 5
```

**Output:**
- Creates `animations/` subdirectory
- Generates `<method_name>_animation.gif` for each method

## Makefile Integration

The scripts are integrated into the Makefile workflow:

```bash
# Single method visualization
make shock_tube_2d_plot_all         # Process all snapshots for default method

# Multi-method comparison
make shock_tube_2d_compare_viz      # 1D analytical plots for all methods
make shock_tube_2d_compare_2d       # 2D spatial plots for all methods
make shock_tube_2d_compare_animate  # Animations for all methods
make shock_tube_2d_compare_all      # Complete workflow: run + viz + 2D + animate
```

## Analytical Solution

The exact Riemann solver computes the analytical solution for the Sod shock tube problem:

**Initial Conditions:**
- Left state (x < 0.5): ρ=1.0, p=1.0, u=0.0
- Right state (x ≥ 0.5): ρ=0.125, p=0.1, u=0.0
- Adiabatic index: γ=1.4

**Solution Features:**
- Left rarefaction wave
- Contact discontinuity
- Right shock wave
- Star region between waves

The solver uses iterative Newton-Raphson to find the pressure in the star region, then computes the full wave structure.

## Dependencies

Required Python packages:
```bash
pip3 install numpy matplotlib pandas
```

For animations:
```bash
pip3 install pillow  # For GIF generation
```

## Output Structure

```
sample/shock_tube_2d/results/
├── gsph_wendland/
│   ├── snapshot_*.csv              # Simulation data
│   ├── plots/                       # 1D analytical comparison plots
│   │   └── snapshot_*_comparison.png
│   └── plots_2d/                    # 2D spatial plots
│       └── snapshot_*_2d.png
├── ssph_wendland/
│   └── ...
├── animations/                      # Animated GIFs
│   ├── gsph_wendland_animation.gif
│   └── ...
└── comparison/                      # Multi-method comparisons
```

## Notes

- 2D simulations average particles over y-direction for 1D comparison
- Scripts automatically extract time and gamma from snapshot metadata
- All scripts support `--no-show` flag for batch processing
- Plots use consistent color schemes across all methods
