# Lane-Emden Architecture Design

## Overview
This document describes the new scalable architecture for Lane-Emden initial conditions supporting multiple configurations for production simulations.

## Directory Structure

```
lane_emden/
├── config/
│   ├── presets/                    # Production-ready configurations
│   │   ├── polytrope_n0_2d.json   # Incompressible (γ→∞)
│   │   ├── polytrope_n0_3d.json
│   │   ├── polytrope_n0.5_2d.json # γ = 3.0
│   │   ├── polytrope_n0.5_3d.json
│   │   ├── polytrope_n1_2d.json   # γ = 2.0 (isothermal limit)
│   │   ├── polytrope_n1_3d.json
│   │   ├── polytrope_n1.5_2d.json # γ = 1.667 (monatomic gas)
│   │   ├── polytrope_n1.5_3d.json
│   │   ├── polytrope_n3_2d.json   # γ = 1.333 (radiation-dominated)
│   │   ├── polytrope_n3_3d.json
│   │   ├── polytrope_n5_2d.json   # γ = 1.2
│   │   └── polytrope_n5_3d.json
│   ├── templates/                  # Base configuration templates
│   │   ├── base_2d.json
│   │   └── base_3d.json
│   └── relaxation_defaults.json   # Default relaxation settings
│
├── data/
│   ├── numerical_solutions/        # Lane-Emden equation solutions
│   │   ├── 2d/
│   │   │   ├── n0.0.dat
│   │   │   ├── n0.5.dat
│   │   │   ├── n1.0.dat
│   │   │   ├── n1.5.dat
│   │   │   ├── n3.0.dat
│   │   │   └── n5.0.dat
│   │   └── 3d/
│   │       ├── n0.0.dat
│   │       ├── n0.5.dat
│   │       ├── n1.0.dat
│   │       ├── n1.5.dat
│   │       ├── n3.0.dat
│   │       └── n5.0.dat
│   └── plots/                      # Visualization of solutions
│       ├── n0.0_comparison.png
│       ├── n1.5_comparison.png
│       └── all_polytropes.png
│
├── scripts/
│   ├── generators/                 # Particle distribution generators
│   │   ├── generate_2d_profile.py
│   │   └── generate_3d_profile.py
│   ├── analysis/
│   │   ├── analyze_relaxation.py
│   │   └── verify_equilibrium.py
│   ├── visualization/
│   │   ├── plot_initial_state.py
│   │   └── animate_relaxation.py
│   └── config_manager.py           # Unified configuration manager
│
├── results/                        # Organized by preset name
│   ├── polytrope_n1.5_3d/
│   │   ├── logs/
│   │   ├── snapshots/
│   │   └── animations/
│   └── polytrope_n3_2d/
│       ├── logs/
│       ├── snapshots/
│       └── animations/
│
├── docs/
│   ├── polytrope_physics.md        # Physics reference
│   ├── preset_guide.md             # Which preset to use when
│   └── configuration_schema.md     # JSON schema documentation
│
└── README.md                       # Quick start guide
```

## Configuration Schema

### Preset Configuration File Structure

Each preset in `config/presets/` follows this schema:

```json
{
  "name": "polytrope_n1.5_3d",
  "description": "3D n=1.5 polytrope (monatomic ideal gas)",
  "dimension": 3,
  
  "physics": {
    "polytropic_index": 1.5,
    "adiabatic_index": 1.66666666666666666666666666666666667,
    "total_mass": 1.0,
    "radius": 1.0,
    "gravitational_constant": 1.0
  },
  
  "initial_conditions": {
    "data_file": "lane_emden/data/numerical_solutions/3d/n1.5.dat",
    "particle_count": 5400,
    "distribution": "fibonacci_sphere"
  },
  
  "relaxation": {
    "enabled": true,
    "steps": 100,
    "convergence_threshold": 1e-6,
    "output_frequency": 10
  },
  
  "simulation": {
    "end_time": 5.0,
    "output_time": 0.05,
    "output_directory": "lane_emden/results/polytrope_n1.5_3d"
  },
  
  "numerical": {
    "neighbor_number": 50,
    "kernel": "wendland",
    "sph_type": "disph",
    "use_gravity": true,
    "periodic": false,
    "iterative_smoothing_length": false
  },
  
  "artificial_viscosity": {
    "alpha": 1.0,
    "use_balsara_switch": true,
    "use_time_dependent_av": true
  }
}
```

## Usage Examples

### Running a Specific Preset

```bash
# Run specific preset
make lane_emden PRESET=polytrope_n1.5_3d

# Run with relaxation only
make lane_emden_relax_only PRESET=polytrope_n3_2d

# Run with custom relaxation steps
make lane_emden PRESET=polytrope_n1.5_3d RELAX_STEPS=200
```

### Creating Custom Configuration

```bash
# Generate new preset from template
python3 lane_emden/scripts/config_manager.py create \
  --name custom_n2_3d \
  --dimension 3 \
  --polytrope 2.0 \
  --particles 10000 \
  --relax-steps 150

# Apply changes to preset
python3 lane_emden/scripts/config_manager.py update \
  --preset polytrope_n1.5_3d \
  --relax-steps 200 \
  --end-time 10.0
```

## Preset Selection Guide

| Preset | n | γ | Physics | Use Case |
|--------|---|---|---------|----------|
| `polytrope_n0` | 0 | ∞ | Incompressible | Dense stellar cores |
| `polytrope_n0.5` | 0.5 | 3.0 | Very stiff | White dwarfs |
| `polytrope_n1` | 1 | 2.0 | Isothermal limit | Molecular clouds |
| `polytrope_n1.5` | 1.5 | 5/3 | Monatomic ideal gas | **Most stars** |
| `polytrope_n3` | 3 | 4/3 | Radiation pressure | Massive stars, AGN |
| `polytrope_n5` | 5 | 1.2 | Degenerate gas | Low-mass stars |

**2D vs 3D:**
- **2D**: Faster testing, cylindrical symmetry problems
- **3D**: Production runs, realistic stellar models

## Migration from Old System

The old system had:
```
sample/lane_emden/lane_emden.json
config/relaxation_config.json
```

New system organizes everything under `lane_emden/` with:
- Multiple presets instead of single config
- Separation of physics parameters from numerical settings
- Clear naming convention: `polytrope_n{index}_{2d|3d}`

### Backward Compatibility

The old `sample/lane_emden/` will be symlinked to `lane_emden/results/legacy/` for compatibility.

## Implementation Phases

### Phase 1: Core Infrastructure (Immediate)
1. Create directory structure
2. Move existing files to new locations
3. Create base templates
4. Create config_manager.py

### Phase 2: Preset Library (Week 1)
1. Generate all polytrope data files (n=0, 0.5, 1, 1.5, 3, 5)
2. Create all preset configurations
3. Test each preset

### Phase 3: Automation (Week 2)
1. Update Makefile with PRESET variable
2. Create visualization scripts
3. Documentation and examples

### Phase 4: Advanced Features (Future)
1. Auto-convergence detection
2. Parameter optimization
3. Batch processing multiple presets
4. CI/CD testing of all presets
