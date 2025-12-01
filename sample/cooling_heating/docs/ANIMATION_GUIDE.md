# CNM Thermal Relaxation - Animation and Plotting Guide

## Animation Scripts Location

The animation and plotting scripts are in:
```
sample/cooling_heating/scripts/
├── animate_cnm_relaxation.py   # Create GIF animation
└── plot_cnm_relaxation.py      # Create static plots
```

## Quick Commands

### 1. Generate Plots
```bash
cd /Users/guo/Downloads/sphcode
make -f sample/cooling_heating/Makefile.cooling cnm_plot
```
or directly:
```bash
python3 sample/cooling_heating/scripts/plot_cnm_relaxation.py sample/cooling_heating/results/cnm_relaxation
```

### 2. Create Animation
```bash
cd /Users/guo/Downloads/sphcode
make -f sample/cooling_heating/Makefile.cooling cnm_animate
```
or directly:
```bash
python3 sample/cooling_heating/scripts/animate_cnm_relaxation.py sample/cooling_heating/results/cnm_relaxation
```

### 3. One-Shot (Run + Plot + Animate)
```bash
cd /Users/guo/Downloads/sphcode
make -f sample/cooling_heating/Makefile.cooling cnm_all
```

## Output Locations

After running, you'll find:

### Simulation Data
```
sample/cooling_heating/results/cnm_relaxation/
├── snapshot_0000.csv
├── snapshot_0001.csv
├── ...
└── energy.dat
```

### Plots and Animation
```
sample/cooling_heating/results/cnm_relaxation/plots/
├── temperature_evolution.png   # T(t) showing cooling
├── phase_diagram.png            # T vs n_H comparison
└── thermal_evolution.gif        # Animated evolution
```

## What the Animation Shows

The animation (`thermal_evolution.gif`) displays 4 panels:

1. **Temperature (top-left)**: Shows thermal relaxation from 107 K → ~25 K
   - Blue line: SPH simulation
   - Gray dashed: Initial T = 107 K
   - Red dashed: K&I equilibrium ~25 K (from Figure 1a)

2. **Pressure (top-right)**: Pressure evolution as system cools

3. **Density (bottom-left)**: Should remain ~10 cm⁻³ (uniform CNM)
   - Green line: SPH simulation
   - Gray dashed: Target density

4. **Internal Energy (bottom-right)**: Energy evolution from cooling

## Requirements

Make sure you have Python packages installed:
```bash
pip install numpy matplotlib pillow
```

## Checking Progress

To see if simulation is still running:
```bash
ps aux | grep sph
```

To see how many snapshots exist:
```bash
ls sample/cooling_heating/results/cnm_relaxation/snapshot_*.csv | wc -l
```

To check latest snapshot time:
```bash
tail -1 sample/cooling_heating/results/cnm_relaxation/energy.dat
```

## Expected Results

- **Initial**: T = 107 K, n = 10 cm⁻³ (uniform CNM)
- **Final**: T ≈ 20-30 K (K&I equilibrium at n=10 cm⁻³)
- **Cooling timescale**: Few to tens of code units
- **Density**: Should stay approximately constant (hydrostatic)

## View Animation

After creating the animation:
```bash
open sample/cooling_heating/results/cnm_relaxation/plots/thermal_evolution.gif
```

Or on Linux:
```bash
xdg-open sample/cooling_heating/results/cnm_relaxation/plots/thermal_evolution.gif
```
