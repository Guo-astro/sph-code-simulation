# Hydrostatic & KHI Infrastructure - Setup Complete ✓

## Overview
Successfully created comprehensive comparison infrastructure for both Hydrostatic and Kelvin-Helmholtz Instability benchmarks, matching the pairing instability setup.

## Created Files Summary

### Hydrostatic Test (sample/hydrostatic/)
**Configuration Presets** (10 files):
- Cubic Spline: `hydrostatic_gsph.json`, `hydrostatic_ssph.json`, `hydrostatic_disph.json`, `hydrostatic_gdisph.json`, `hydrostatic_gdisph_balsara.json`
- Wendland: `hydrostatic_gsph_wendland.json`, `hydrostatic_ssph_wendland.json`, `hydrostatic_disph_wendland.json`, `hydrostatic_gdisph_wendland.json`, `hydrostatic_gdisph_balsara_wendland.json`

**Build System**:
- `Makefile.hydrostatic` - Complete build targets for all workflows

**Visualization Scripts** (4 files):
- `scripts/compare_hydrostatic.py` - Method comparison plots
- `scripts/animate_hydrostatic.py` - Method comparison animation
- `scripts/compare_kernels.py` - Kernel comparison plots
- `scripts/animate_kernel_comparison.py` - Kernel comparison animation

**Documentation**:
- `README.md` - Complete usage guide

### KHI Test (sample/khi/)
**Configuration Presets** (10 files):
- Cubic Spline: `khi_gsph.json`, `khi_ssph.json`, `khi_disph.json`, `khi_gdisph.json`, `khi_gdisph_balsara.json`
- Wendland: `khi_gsph_wendland.json`, `khi_ssph_wendland.json`, `khi_disph_wendland.json`, `khi_gdisph_wendland.json`, `khi_gdisph_balsara_wendland.json`

**Build System**:
- `Makefile.khi` - Complete build targets for all workflows

**Visualization Scripts** (4 files):
- `scripts/compare_khi.py` - Method comparison plots
- `scripts/animate_khi.py` - Method comparison animation
- `scripts/compare_kernels.py` - Kernel comparison plots
- `scripts/animate_kernel_comparison.py` - Kernel comparison animation

**Documentation**:
- `README.md` - Complete usage guide

### Main Makefile Updates
- Added `hydrostatic_help` and `khi_help` to main help menu
- Included `sample/hydrostatic/Makefile.hydrostatic`
- Included `sample/khi/Makefile.khi`

## Test Configurations

### Hydrostatic Test
- **Resolution**: 64×64 = 4,096 particles
- **Domain**: [-0.5, 0.5]² with periodic boundaries
- **Duration**: t = 0 → 8.0
- **Output**: Every 0.2 time units (41 snapshots)
- **Physics**: Tests pressure gradient accuracy and equilibrium maintenance

### KHI Test
- **Resolution**: 256×256 = 65,536 particles
- **Domain**: [0, 1]² with periodic boundaries
- **Duration**: t = 0 → 3.0
- **Output**: Every 0.1 time units (31 snapshots)
- **Physics**: Tests shear instability, vortex evolution, mixing
- **Features**: Time-dependent artificial viscosity enabled

## Available Workflows

### Single Method Run
```bash
# Hydrostatic
make hydrostatic_run
make hydrostatic_run HYDROSTATIC_PRESET=hydrostatic_gsph

# KHI
make khi_run
make khi_run KHI_PRESET=khi_gsph
```

### Multi-Method Comparison (5 runs)
```bash
# Hydrostatic
make hydrostatic_compare_all      # Complete workflow
make hydrostatic_compare_run      # Just simulations
make hydrostatic_compare_viz      # Just plots
make hydrostatic_compare_animate  # Just animation

# KHI
make khi_compare_all              # Complete workflow
make khi_compare_run              # Just simulations
make khi_compare_viz              # Just plots
make khi_compare_animate          # Just animation
```

### Kernel Comparison (10 runs)
```bash
# Hydrostatic
make hydrostatic_kernel_all       # Complete workflow
make hydrostatic_kernel_run       # 5 methods × 2 kernels
make hydrostatic_kernel_viz       # Side-by-side plots
make hydrostatic_kernel_animate   # Comparison animation

# KHI
make khi_kernel_all               # Complete workflow
make khi_kernel_run               # 5 methods × 2 kernels
make khi_kernel_viz               # Side-by-side plots
make khi_kernel_animate           # Comparison animation
```

### Cleanup
```bash
# Hydrostatic
make hydrostatic_clean            # Clean all
make hydrostatic_compare_clean    # Clean comparison only
make hydrostatic_kernel_clean     # Clean kernel comparison only

# KHI
make khi_clean                    # Clean all
make khi_compare_clean            # Clean comparison only
make khi_kernel_clean             # Clean kernel comparison only
```

## Makefile Target Summary

### Hydrostatic Targets
- `hydrostatic_run` - Single method run
- `hydrostatic_compare_run` - Run all 5 methods (cubic spline)
- `hydrostatic_compare_viz` - Generate comparison plots
- `hydrostatic_compare_animate` - Generate comparison animation
- `hydrostatic_compare_all` - Complete comparison workflow
- `hydrostatic_kernel_run` - Run 10 simulations (5 methods × 2 kernels)
- `hydrostatic_kernel_viz` - Generate kernel comparison plots
- `hydrostatic_kernel_animate` - Generate kernel comparison animation
- `hydrostatic_kernel_all` - Complete kernel comparison workflow
- `hydrostatic_clean` - Clean all results
- `hydrostatic_compare_clean` - Clean comparison results
- `hydrostatic_kernel_clean` - Clean kernel results
- `hydrostatic_help` - Show help

### KHI Targets
- `khi_run` - Single method run
- `khi_compare_run` - Run all 5 methods (cubic spline)
- `khi_compare_viz` - Generate comparison plots
- `khi_compare_animate` - Generate comparison animation
- `khi_compare_all` - Complete comparison workflow
- `khi_kernel_run` - Run 10 simulations (5 methods × 2 kernels)
- `khi_kernel_viz` - Generate kernel comparison plots
- `khi_kernel_animate` - Generate kernel comparison animation
- `khi_kernel_all` - Complete kernel comparison workflow
- `khi_clean` - Clean all results
- `khi_compare_clean` - Clean comparison results
- `khi_kernel_clean` - Clean kernel results
- `khi_help` - Show help

## Visualization Features
- Hot colormap with log₁₀(density) for better dynamic range
- Black background for enhanced contrast
- 5-panel layout for method comparison
- 5×2 grid for kernel comparison (Wendland left, Cubic Spline right)
- Automatic snapshot selection at key time points
- Progress bars (with tqdm installed)
- Publication-quality output (150 DPI for plots)

## Performance Estimates

### Hydrostatic (4,096 particles)
- Single simulation: ~1-3 seconds
- 5 method comparison: ~5-15 seconds
- 10 kernel runs: ~10-30 seconds
- Plot generation: ~10-30 seconds
- Animation: ~1-2 minutes

### KHI (65,536 particles)
- Single simulation: ~30-90 seconds
- 5 method comparison: ~2.5-7.5 minutes
- 10 kernel runs: ~5-15 minutes
- Plot generation: ~30-60 seconds
- Animation: ~2-5 minutes

## Scientific Applications

### Hydrostatic Test
- Validates pressure gradient accuracy
- Tests equilibrium maintenance
- Measures dissipation characteristics
- Compares long-term stability
- Assesses kernel smoothness effects on stability

### KHI Test
- Tests shear flow instability evolution
- Compares vortex resolution quality
- Evaluates mixing and turbulence development
- Measures numerical diffusion
- Assesses small-scale structure preservation
- Compares kernel impact on instability growth

## Directory Structure
```
sample/
├── hydrostatic/
│   ├── config/presets/          # 10 JSON presets
│   ├── scripts/                 # 4 Python scripts
│   ├── results/                 # Output directories (created on run)
│   ├── Makefile.hydrostatic
│   ├── README.md
│   └── hydrostatic.json         # Active config (updated by presets)
│
└── khi/
    ├── config/presets/          # 10 JSON presets
    ├── scripts/                 # 4 Python scripts
    ├── results/                 # Output directories (created on run)
    ├── Makefile.khi
    ├── README.md
    └── khi.json                 # Active config (updated by presets)
```

## Integration with Main Build System
Both benchmarks are now fully integrated into the main Makefile:
```bash
make                          # Shows all available workflows
make hydrostatic_help         # Hydrostatic help
make khi_help                 # KHI help
```

## Requirements
- DIM=2 build configuration (required for both tests)
- Python 3 with numpy, matplotlib
- Optional: tqdm (progress bars), pillow (animations)

## Next Steps
1. Run initial tests: `make hydrostatic_compare_all` and `make khi_compare_all`
2. Analyze method differences in equilibrium vs instability scenarios
3. Compare kernel effects across different flow regimes
4. Generate publication figures
5. Consider adding more benchmarks (e.g., Gresho vortex, shock tube variations)

## Status: Ready for Use ✓
All infrastructure is complete and ready for scientific simulations. Both benchmarks now have the same comprehensive tooling as the pairing instability test.
