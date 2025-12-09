# Lane-Emden Initial Conditions

**⚠️ SINGLE SOURCE OF TRUTH**: All configurations are managed through **presets** in `config/presets/`. 

The C++ config file (`sample/lane_emden/lane_emden.json`) is **auto-generated** from presets - don't edit it directly!

## 🚀 Quick Start (New Workflow)

### Option 1: Using Makefile (Easiest)
```bash
# Run relaxation with 10000 steps
make lane_emden_relax_only RELAX_STEPS=10000

# Run full simulation (relaxation + evolution)
make lane_emden RELAX_STEPS=10000

# Use different preset
make lane_emden PRESET=polytrope_n1_5_2d RELAX_STEPS=5000

# List available presets
make lane_emden_list

# Get help
make lane_emden_help
```

### Option 2: Using Config Manager Directly
```bash
# 1. Update preset (Single Source of Truth)
python3 lane_emden_config_manager.py update --preset polytrope_n1_5_3d --relax-steps 10000

# 2. Generate C++ config from preset
python3 lane_emden_config_manager.py apply --preset polytrope_n1_5_3d

# 3. Run simulation
./build/sph lane_emden
```

**Note**: `lane_emden_config_manager.py` is a symlink to `lane_emden/scripts/lane_emden_config_manager.py`

**📖 Full workflow guide**: See [WORKFLOW.md](WORKFLOW.md)

## 📋 Configuration Management

### List Presets
```bash
python3 lane_emden_config_manager.py list
```

### Update Preset Parameters
```bash
python3 lane_emden_config_manager.py update --preset polytrope_n1_5_3d --relax-steps 5000
python3 lane_emden_config_manager.py update --preset polytrope_n1_5_3d --end-time 10.0
```

### Apply Preset (Generate C++ Config)
```bash
python3 lane_emden_config_manager.py apply --preset polytrope_n1_5_3d
```

### Create New Preset
```bash
python3 lane_emden_config_manager.py create \
  --name polytrope_n3_3d \
  --polytrope 3.0 \
  --dimension 3 \
  --relax-steps 1000
```

## Available Presets

Currently available:
- `polytrope_n1.5_3d` - 3D monatomic ideal gas (γ=5/3) - **Most common**
- `polytrope_n1.5_2d` - 2D version for testing

Coming soon:
- `polytrope_n0_3d` - Incompressible fluid
- `polytrope_n1_3d` - Isothermal (γ=2)
- `polytrope_n3_3d` - Radiation-dominated (γ=4/3)
- `polytrope_n5_3d` - Degenerate gas (γ=1.2)

See `docs/LANE_EMDEN_ARCHITECTURE.md` for full list and physics details.

## Directory Structure

```
lane_emden/
├── config/
│   ├── presets/           # Production-ready configurations
│   ├── templates/         # Base templates for creating new configs
│   └── relaxation_defaults.json
├── data/
│   ├── numerical_solutions/  # Lane-Emden equation solutions
│   │   ├── 2d/
│   │   └── 3d/
│   └── plots/
├── scripts/
│   ├── generators/        # Particle distribution generators
│   ├── analysis/          # Analysis tools
│   ├── visualization/     # Plotting and animation
│   └── config_manager.py  # Configuration management
├── results/               # Organized by preset name
└── docs/                  # Documentation
```

## Creating Custom Configurations

### Using Config Manager

```bash
# Create new preset from template
python3 lane_emden/scripts/config_manager.py create \
  --name custom_n3_3d \
  --dimension 3 \
  --polytrope 3.0 \
  --gamma 1.333333 \
  --particles 10000

# Update existing preset
python3 lane_emden/scripts/config_manager.py update \
  --preset polytrope_n1.5_3d \
  --relax-steps 200 \
  --end-time 10.0

# List all presets
python3 lane_emden/scripts/config_manager.py list
```

### Manual Configuration

Edit JSON files in `config/presets/`:

```json
{
  "name": "custom_config",
  "dimension": 3,
  "physics": {
    "polytropic_index": 1.5,
    "adiabatic_index": 1.667
  },
  "relaxation": {
    "steps": 100
  }
}
```

## Polytrope Physics Reference

| n | γ | Description | Physical Systems |
|---|---|-------------|------------------|
| 0 | ∞ | Incompressible | Dense stellar cores |
| 0.5 | 3.0 | Very stiff EOS | White dwarfs (non-relativistic) |
| 1.0 | 2.0 | Isothermal limit | Molecular clouds |
| **1.5** | **5/3** | **Monatomic gas** | **Main sequence stars** ⭐ |
| 3.0 | 4/3 | Radiation pressure | Massive stars, AGN |
| 5.0 | 1.2 | Degenerate gas | Low-mass stars, brown dwarfs |

**Relation**: γ = 1 + 1/n (for n > 0)

## Migration Notes

**Old System** (deprecated):
- Config: `config/relaxation_config.json`
- Sample: `sample/lane_emden/lane_emden.json`
- Results: `sample/lane_emden/results/`

**New System**:
- Presets: `lane_emden/config/presets/`
- Results: `lane_emden/results/{preset_name}/`
- Backward compatible via symlinks

## Documentation

- `docs/LANE_EMDEN_ARCHITECTURE.md` - Full architecture design
- `docs/CHECKPOINT_SYSTEM.md` - Checkpoint/resume system for long relaxation runs
- `docs/RELAXATION_ONLY_MODE.md` - Relaxation-only mode documentation
- `WORKFLOW.md` - Single Source of Truth workflow guide

## Resume from Checkpoint

```bash
# Resume from specific checkpoint file
make -f lane_emden/Makefile.lane_emden lane_emden_resume \
  CHECKPOINT=lane_emden/results/polytrope_n1.5_3d/snapshot_0009.csv

# See docs/CHECKPOINT_SYSTEM.md for full documentation
```

## Examples

### Example 1: Standard Star Simulation
```bash
make lane_emden PRESET=polytrope_n1.5_3d
```

### Example 2: Quick Relaxation Test
```bash
make lane_emden_relax_only PRESET=polytrope_n1.5_3d RELAX_STEPS=20
```

### Example 3: Radiation-Dominated System
```bash
# Coming soon
make lane_emden PRESET=polytrope_n3_3d
```

## Contributing

When adding new polytrope indices:

1. Generate numerical solution: `python3 scripts/generators/solve_lane_emden.py --n 2.0`
2. Create preset config: `python3 scripts/config_manager.py create --polytrope 2.0`
3. Test: `make lane_emden_relax_only PRESET=polytrope_n2.0_3d`
4. Document physics in `docs/polytrope_physics.md`
