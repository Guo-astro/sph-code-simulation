# Lane-Emden Multi-Method Comparison Guide

This guide explains how to run relaxation followed by multi-method hydrodynamic simulations (SSPH, GSPH, DISPH, GDISPH) for Lane-Emden polytropes.

## 📋 Prerequisites

1. **Build for 3D**: Lane-Emden requires DIM=3
   ```bash
   cd build
   cmake -DSPH_DIM=3 ..
   make -j8
   cd ..
   ```

2. **Verify available presets**:
   ```bash
   make -f lane_emden/Makefile.lane_emden lane_emden_list
   ```

## 🚀 Quick Start: Complete Workflow

### Option 1: Automated Complete Workflow (Recommended)
Run relaxation + all 4 SPH methods + visualizations in one command:

```bash
make -f lane_emden/Makefile.lane_emden lane_emden_workflow RELAX_STEPS=10000
```

This will:
1. ✅ Run relaxation (10,000 steps) to prepare equilibrium initial conditions
2. ✅ Run SSPH simulation from relaxed state
3. ✅ Run GSPH simulation from relaxed state
4. ✅ Run DISPH simulation from relaxed state
5. ✅ Run GDISPH simulation from relaxed state
6. ✅ Generate comparison plots
7. ✅ Generate comparison animations

### Option 2: Step-by-Step Workflow
For more control, run each step separately:

```bash
# Step 1: Run relaxation to prepare initial conditions
make -f lane_emden/Makefile.lane_emden lane_emden_relax RELAX_STEPS=10000

# Step 2: Run all 4 SPH methods + visualizations
make -f lane_emden/Makefile.lane_emden lane_emden_compare_all
```

### Option 3: Manual Step-by-Step
For maximum control:

```bash
# Step 1: Relaxation
make -f lane_emden/Makefile.lane_emden lane_emden_relax RELAX_STEPS=10000

# Step 2: Run all SPH methods
make -f lane_emden/Makefile.lane_emden lane_emden_compare_run

# Step 3: Generate plots
make -f lane_emden/Makefile.lane_emden lane_emden_compare_viz

# Step 4: Generate animations
make -f lane_emden/Makefile.lane_emden lane_emden_compare_animate
```

## 📊 Understanding the Workflow

### Phase 1: Relaxation
**Purpose**: Create stable, equilibrium initial conditions

```bash
make -f lane_emden/Makefile.lane_emden lane_emden_relax RELAX_STEPS=10000
```

**What happens**:
- Particles are initialized from Lane-Emden analytical solution
- System is evolved with damping to remove transients
- Final relaxed state saved to `snapshot_final.csv`
- No gravity, no hydrodynamics - just settling into equilibrium

**Output**: `lane_emden/results/polytrope_n1.5_3d/snapshot_final.csv`

### Phase 2: Multi-Method Hydrodynamic Simulations
**Purpose**: Compare different SPH formulations on the same initial conditions

```bash
make -f lane_emden/Makefile.lane_emden lane_emden_compare_run
```

**What happens**:
- Reads `snapshot_final.csv` as initial conditions for ALL methods
- Runs 4 independent simulations:
  1. **SSPH** (Standard SPH)
  2. **GSPH** (Godunov SPH)
  3. **DISPH** (Density-Independent SPH)
  4. **GDISPH** (Godunov DISPH hybrid)
- Each method evolves the same relaxed star with gravity enabled
- Results saved to separate directories for comparison

**Output**:
```
lane_emden/results/polytrope_n1.5_3d/comparison/
├── ssph/         # SSPH results
├── gsph/         # GSPH results
├── disph/        # DISPH results
└── gdisph/       # GDISPH results
```

### Phase 3: Visualization
**Purpose**: Compare methods visually

```bash
make -f lane_emden/Makefile.lane_emden lane_emden_compare_viz      # Plots
make -f lane_emden/Makefile.lane_emden lane_emden_compare_animate  # Animations
```

**Output**:
```
lane_emden/results/polytrope_n1.5_3d/comparison/
├── plots/                    # Comparison plots
└── animations/               # Comparison animations
    └── comparison_animation.gif
```

## 🎯 Available Targets

### Relaxation Targets
| Target | Description |
|--------|-------------|
| `lane_emden_relax` | Run relaxation only (default: 40000 steps) |
| `lane_emden_relax RELAX_STEPS=10000` | Override relaxation steps |
| `lane_emden_resume` | Resume from latest checkpoint |

### Multi-Method Comparison Targets
| Target | Description |
|--------|-------------|
| `lane_emden_compare_run` | Run all 4 SPH methods |
| `lane_emden_compare_viz` | Generate comparison plots |
| `lane_emden_compare_animate` | Generate comparison animation |
| `lane_emden_compare_all` | Run methods + visualizations |
| `lane_emden_compare_clean` | Clean comparison results (keeps relaxed ICs) |

### Complete Workflow Targets
| Target | Description |
|--------|-------------|
| `lane_emden_workflow` | Relaxation + Multi-Method + Visualizations |

### Management Targets
| Target | Description |
|--------|-------------|
| `lane_emden_list` | List available presets |
| `lane_emden_clean` | Clean all results |
| `lane_emden_help` | Show help |

## 📝 Examples

### Example 1: Standard Workflow (Recommended)
```bash
# Complete workflow with 10,000 relaxation steps
make -f lane_emden/Makefile.lane_emden lane_emden_workflow RELAX_STEPS=10000
```

### Example 2: Quick Test (Short Relaxation)
```bash
# Fast test with minimal relaxation
make -f lane_emden/Makefile.lane_emden lane_emden_workflow RELAX_STEPS=1000
```

### Example 3: High-Quality Run (Long Relaxation)
```bash
# Production run with thorough relaxation
make -f lane_emden/Makefile.lane_emden lane_emden_workflow RELAX_STEPS=50000
```

### Example 4: Different Polytrope (n=1.5, 2D)
```bash
# 2D version (requires DIM=2 build)
make -f lane_emden/Makefile.lane_emden lane_emden_workflow \
  PRESET=polytrope_n1_5_2d \
  RELAX_STEPS=10000
```

### Example 5: Only Run Specific Methods
```bash
# Run relaxation once
make -f lane_emden/Makefile.lane_emden lane_emden_relax RELAX_STEPS=10000

# Then manually run only SSPH and GSPH by editing the Makefile target
# Or wait for individual method targets to be implemented
```

## 🔍 Checking Results

### Directory Structure After Complete Workflow
```
lane_emden/results/polytrope_n1.5_3d/
├── snapshot_final.csv              # Relaxed initial conditions
├── checkpoints/                     # Relaxation checkpoints
├── comparison/
│   ├── ssph/
│   │   ├── snapshot_0000.csv
│   │   ├── snapshot_0001.csv
│   │   └── ...
│   ├── gsph/
│   │   └── ...
│   ├── disph/
│   │   └── ...
│   ├── gdisph/
│   │   └── ...
│   ├── plots/
│   │   ├── density_comparison_t0.00.png
│   │   ├── density_comparison_t1.00.png
│   │   └── ...
│   └── animations/
│       └── comparison_animation.gif
```

### Viewing Results
1. **Check snapshots**: Each method saves time snapshots in its directory
2. **View plots**: Open `comparison/plots/*.png` files
3. **Watch animations**: Open `comparison/animations/comparison_animation.gif`

## ⚙️ Configuration Details

### Default Settings (from preset)
- **Polytropic index**: n = 1.5 (γ = 5/3, monatomic ideal gas)
- **Particles**: 5400
- **Relaxation steps**: 40000 (override with `RELAX_STEPS`)
- **Simulation time**: 5.0 (dimensionless units)
- **Output frequency**: Every 0.05 time units
- **Kernel**: Wendland
- **Gravity**: Disabled during relaxation, enabled for hydrodynamic runs

### Changing Settings
You can modify the preset JSON files directly:
```bash
# Edit preset
vi lane_emden/config/presets/polytrope_n1_5_3d.json

# Or use config manager
python3 lane_emden/scripts/lane_emden_config_manager.py update \
  --preset polytrope_n1_5_3d \
  --relax-steps 20000
```

## 🐛 Troubleshooting

### Error: "Lane-Emden requires DIM=3 build"
**Solution**: Rebuild for 3D
```bash
cd build
cmake -DSPH_DIM=3 ..
make -j8
cd ..
```

### Error: "No relaxed initial conditions found"
**Solution**: Run relaxation first
```bash
make -f lane_emden/Makefile.lane_emden lane_emden_relax RELAX_STEPS=10000
```

### Error: "Missing results directories"
**Solution**: Run comparison simulations first
```bash
make -f lane_emden/Makefile.lane_emden lane_emden_compare_run
```

### Relaxation Not Converging
**Solution**: Increase relaxation steps
```bash
make -f lane_emden/Makefile.lane_emden lane_emden_relax RELAX_STEPS=50000
```

## 📚 Physics Background

### What is a Lane-Emden Polytrope?
A self-gravitating sphere in hydrostatic equilibrium with polytropic equation of state:
- P = K ρ^γ
- Where γ = 1 + 1/n (n = polytropic index)

### Why n=1.5 (γ=5/3)?
- Represents monatomic ideal gas
- Most relevant for stellar interiors
- Main sequence stars, stellar structure calculations

### Why Relaxation?
- Initial particle distribution from analytical solution has small errors
- Relaxation removes transients and achieves numerical equilibrium
- Provides consistent initial conditions for all SPH methods

### Why Compare Methods?
- Different SPH formulations have different conservation properties
- SSPH: Simple, fast, but has issues with discontinuities
- GSPH: Better shock capturing, exact momentum conservation
- DISPH: Better with density gradients
- GDISPH: Hybrid approach combining advantages

## 🎓 Scientific Use Cases

1. **Method Validation**: Compare SPH methods against known analytical solution
2. **Stability Analysis**: Test which methods maintain equilibrium best
3. **Convergence Studies**: Vary particle number and relaxation steps
4. **Shock Testing**: Add perturbations to trigger collapse/oscillations
5. **Kernel Comparison**: Test different kernels (Wendland, Cubic Spline, etc.)

## 📖 See Also

- `lane_emden/README.md` - General Lane-Emden documentation
- `lane_emden/WORKFLOW.md` - Detailed workflow guide
- `sample/khi/Makefile.khi` - Similar multi-method comparison for KHI
- `docs/GDISPH_IMPLEMENTATION.md` - GDISPH technical details
