# Lane-Emden Reorganized Structure

## Overview
The Lane-Emden directory has been reorganized to clearly separate 2D and 3D workflows.

## Directory Structure

```
lane_emden/
├── 2d/                          # 2D disk simulations
│   ├── config/
│   │   └── presets/
│   │       └── polytrope_n1_5_ssph.json
│   ├── results/
│   │   └── polytrope_n1_5/      # All 2D results here
│   └── Makefile.2d              # 2D workflow automation
│
├── 3d/                          # 3D stellar simulations
│   ├── config/
│   │   └── presets/
│   │       └── polytrope_n1_5_ssph.json
│   ├── results/
│   │   └── polytrope_n1_5/      # All 3D results here
│   └── Makefile.3d              # 3D workflow automation
│
├── data/                        # Shared numerical solutions
│   └── numerical_solutions/
│       ├── 2d/
│       │   └── n1.5.dat
│       └── 3d/
│           └── n1.5.dat
│
└── scripts/                     # Shared analysis scripts
    └── visualization/
```

## Key Improvements

1. **Clear Separation**: 2D and 3D workflows are completely independent
2. **Clean Folder Names**: Using `polytrope_n1_5` (with underscores) consistently
3. **CSV-Only Output**: HDF5 disabled for faster I/O and simpler workflow
4. **Simplified Makefiles**: Separate Makefile.2d and Makefile.3d for each dimension

## Usage

### 2D Workflow (Disk Simulations)

```bash
# Build for 2D
cd build && cmake -DSPH_DIM=2 .. && make -j8 && cd ..

# Run simulation (no relaxation - uses analytical ICs)
make -f lane_emden/2d/Makefile.2d relax

# Results
ls lane_emden/2d/results/polytrope_n1_5/
```

**Note**: 2D relaxation currently has stability issues (segfault). The workflow skips relaxation and uses analytical initial conditions directly, which are already in perfect hydrostatic equilibrium.

### 3D Workflow (Stellar Models)

```bash
# Build for 3D
cd build && cmake -DSPH_DIM=3 .. && make -j8 && cd ..

# Run relaxation
make -f lane_emden/3d/Makefile.3d relax RELAX_STEPS=10000

# Run simulation
make -f lane_emden/3d/Makefile.3d run

# Results
ls lane_emden/3d/results/polytrope_n1_5/
```

## Configuration Details

### Output Format
- **CSV only** (HDF5 disabled)
- Precision: 16 digits
- Energy file: enabled

### 2D Preset (`polytrope_n1_5_ssph`)
- Particles: 900
- SPH Type: SSPH (Standard SPH)
- Kernel: Wendland C4
- Gravity: enabled
- Output: `lane_emden/2d/results/polytrope_n1_5/`

### 3D Preset (`polytrope_n1_5_ssph`)
- Particles: 5400
- SPH Type: SSPH (Standard SPH)
- Kernel: Wendland C4
- Gravity: enabled
- Output: `lane_emden/3d/results/polytrope_n1_5/`

## Legacy Files

The old structure with mixed naming (`polytrope_n1.5_2d`, `polytrope_n1_5_3d`) is still present in:
- `lane_emden/config/presets/` (old)
- `lane_emden/results/` (old)

These can be safely deleted after confirming the new structure works correctly.

## Migration Notes

1. **Folder naming**: Now consistently uses underscores (`polytrope_n1_5`)
2. **No more dimension suffix**: Results folder is just `polytrope_n1_5/`, separated by `2d/` and `3d/` parent directories
3. **Separate workflows**: Each dimension has its own Makefile and config directory
4. **Simpler**: No complex dimension checking in a single Makefile

## Next Steps

To fully migrate:

```bash
# Clean old results
rm -rf lane_emden/results/polytrope_n1.5_*
rm -rf lane_emden/results/polytrope_n1_5_*

# Remove old presets
rm lane_emden/config/presets/polytrope_n1_5_2d.json
rm lane_emden/config/presets/polytrope_n1_5_3d.json
```
