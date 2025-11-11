# Sample Test Cases

This directory contains benchmark test cases for SPH simulations, organized by dimensionality and physics.

## Available Test Cases

### 1D Hydrodynamic Tests
- **shock_tube** - Sod shock tube test (1D Riemann problem)
  - Single method: `make shock_tube_run`
  - Multi-method comparison: `make shock_tube_compare_all`
  - See `sample/shock_tube/Makefile.shock_tube` for all targets

### 2D Hydrodynamic Tests
- **khi** - Kelvin-Helmholtz instability
  - Single method: `make khi_run`
  - Multi-method comparison: `make khi_compare_all`
  - Kernel comparison: `make khi_kernel_all`
  - See `sample/khi/Makefile.khi` for all targets

- **gresho_chan_vortex** - Gresho-Chan vortex test
  - Single method: `make gresho_run`
  - Multi-method comparison: `make gresho_compare_all`
  - Kernel comparison: `make gresho_kernel_all`
  - See `sample/gresho_chan_vortex/Makefile.gresho` for all targets

- **hydrostatic** - Hydrostatic equilibrium test
  - Single method: `make hydrostatic_run`
  - Multi-method comparison: `make hydrostatic_compare_all`
  - Kernel comparison: `make hydrostatic_kernel_all`
  - See `sample/hydrostatic/Makefile.hydrostatic` for all targets

- **pairing_instability** - Pairing instability test
  - Single method: `make pairing_run`
  - Multi-method comparison: `make pairing_compare_all`
  - Kernel comparison: `make pairing_kernel_all`
  - See `sample/pairing_instability/Makefile.pairing_instability` for all targets

### 3D Tests
- **sedov** - Sedov blast wave (3D spherical explosion)
  - Single method: `make sedov_run`
  - Multi-method comparison: `make sedov_compare_all`
  - Kernel comparison: `make sedov_kernel_all`
  - See `sample/sedov/Makefile.sedov` for all targets

- **evrard** - Evrard collapse test (equilibrium after shock)
  - Coming soon

- **lane_emden** - Lane-Emden hydrostatic sphere (located in `lane_emden/`)
  - Single method: `make lane_emden`
  - Animation: `make lane_emden_animate`
  - Resume: `make lane_emden_resume`
  - See `lane_emden/Makefile.lane_emden` for all targets

## Makefile System

All test cases follow a consistent Makefile structure with standardized targets:

### Standard Targets

**Single Method Execution:**
- `make <test>_run` - Run simulation with default method/preset
- `make <test>_clean` - Clean results directory
- `make <test>_help` - Show available targets

**Multi-Method Comparison (5 SPH methods with default kernel):**
- `make <test>_compare_run` - Run GSPH, SSPH, DISPH, GDISPH, GDISPH+Balsara
- `make <test>_compare_viz` - Generate comparison plots
- `make <test>_compare_animate` - Generate comparison animations
- `make <test>_compare_all` - Run all comparisons + visualizations
- `make <test>_compare_clean` - Clean comparison results

**Kernel Comparison (5 methods × 2 kernels = 10 runs):**
- `make <test>_kernel_run` - Run all methods with Wendland and Cubic Spline kernels
- `make <test>_kernel_viz` - Generate kernel comparison plots
- `make <test>_kernel_animate` - Generate kernel comparison animations
- `make <test>_kernel_all` - Run all + visualizations
- `make <test>_kernel_clean` - Clean kernel comparison results

### Examples

```bash
# Run single simulation with default settings
make khi_run
make sedov_run

# Run complete multi-method comparison workflow
make khi_compare_all          # Runs 5 methods + generates plots + animations
make sedov_compare_all

# Run complete kernel comparison (10 simulations)
make khi_kernel_all           # 5 methods × 2 kernels + plots + animations
make sedov_kernel_all

# Individual steps
make khi_compare_run          # Just run the 5 simulations
make khi_compare_viz          # Just generate comparison plots
make khi_compare_animate      # Just generate animations

# Clean up results
make khi_clean                # Remove all results
make khi_compare_clean        # Remove only comparison results
make khi_kernel_clean         # Remove only kernel comparison results
```

## Test Case Structure

Each test case follows this directory structure:

```
sample/<test>/
├── Makefile.<test>           # Makefile with all targets
├── README.md                 # Test-specific documentation
├── <test>.json               # Active configuration (auto-generated)
├── config/
│   └── presets/              # Preset configurations
│       ├── <test>_gsph.json
│       ├── <test>_ssph.json
│       ├── <test>_disph.json
│       ├── <test>_gdisph.json
│       ├── <test>_gdisph_balsara.json
│       ├── <test>_gsph_wendland.json (if kernel comparison supported)
│       └── ...
├── scripts/                  # Visualization and analysis scripts
│   ├── compare_<test>.py     # Multi-method comparison plots
│   ├── animate_<test>.py     # Comparison animations
│   ├── compare_kernels.py    # Kernel comparison plots (if supported)
│   └── animate_kernel_comparison.py
└── results/                  # Output directory (created on run)
    ├── gsph/
    ├── ssph/
    ├── disph/
    ├── gdisph/
    ├── gdisph_balsara/
    ├── comparison/           # Comparison plots and animations
    └── kernel_comparison/    # Kernel comparison results (if applicable)
```

## SPH Methods

All test cases support comparison between these SPH methods:

1. **GSPH** - Godunov SPH (Inutsuka 2002)
2. **SSPH** - Standard SPH (Monaghan 1992)
3. **DISPH** - Density-Independent SPH (Saitoh & Makino 2013)
4. **GDISPH** - Generalized Density-Independent SPH (Hopkins 2013)
5. **GDISPH+Balsara** - GDISPH with Balsara switch for artificial viscosity

## Kernel Functions

Most 2D/3D tests support comparison between:

1. **Cubic Spline** (M4 kernel) - Traditional SPH kernel
2. **Wendland C4** - Compact support kernel with better stability

## Dimension Requirements

- **1D tests** require `DIM=1` build (shock_tube)
- **2D tests** require `DIM=2` build (khi, gresho, hydrostatic, pairing_instability)
- **3D tests** require `DIM=3` build (sedov, evrard, lane_emden)

To change build dimension:
```bash
# Change to 2D
sed -i '' 's/#define DIM [0-9]/#define DIM 2/' include/defines.hpp
cd build && make -j8 && cd ..

# Change to 3D
sed -i '' 's/#define DIM [0-9]/#define DIM 3/' include/defines.hpp
cd build && make -j8 && cd ..
```

## Quick Start

```bash
# 1. Ensure correct build dimension
grep "^#define DIM" include/defines.hpp

# 2. Run a single test
make khi_run

# 3. Run complete comparison workflow
make khi_compare_all

# 4. Run kernel comparison
make khi_kernel_all

# 5. Get help for specific test
make khi_help
make sedov_help
```

## Results and Visualization

Results are saved to `sample/<test>/results/` with the following structure:

- **Snapshots**: `snapshot_NNNN.csv` - Particle data at each output time
- **Plots**: Generated in `results/comparison/` or `results/kernel_comparison/`
- **Animations**: `.gif` files showing time evolution

All visualization scripts are in `sample/<test>/scripts/` and are automatically called by the Makefile targets.

## Notes

- All Makefiles follow consistent naming and structure conventions
- Each test includes comprehensive documentation in its own README
- Preset configurations are stored in `config/presets/` for reproducibility
- Output directories are created automatically by Makefile targets
- Use `make <test>_help` to see all available targets for any test case
