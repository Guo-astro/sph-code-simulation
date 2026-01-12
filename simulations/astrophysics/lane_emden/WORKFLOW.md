# Lane-Emden Single Source of Truth Workflow

## 🎯 Philosophy: Presets are the Single Source of Truth

All Lane-Emden configurations are managed through **presets** in `lane_emden/config/presets/`. The C++ runtime config in `lane_emden/lane_emden.json` is **generated** from presets, not edited directly.

## 📁 Directory Structure

```
lane_emden/
├── config/
│   ├── presets/          ← SINGLE SOURCE OF TRUTH
│   │   ├── polytrope_n1_5_3d.json
│   │   └── polytrope_n1_5_2d.json
│   ├── templates/        ← Templates for creating new presets
│   └── relaxation_defaults.json
├── data/
│   └── numerical_solutions/
│       ├── 2d/n1.5.dat
│       └── 3d/n1.5.dat
├── scripts/
├── results/
└── docs/

lane_emden/
└── lane_emden.json       ← GENERATED from presets (DO NOT EDIT)
```

## 🔧 Workflow

### 1. List Available Presets
```bash
python3 lane_emden_config_manager.py list
python3 lane_emden_config_manager.py list --verbose  # Detailed info
```

### 2. Update a Preset (Single Source of Truth)
```bash
# Update relaxation steps
python3 lane_emden_config_manager.py update --preset polytrope_n1_5_3d --relax-steps 10000

# Update simulation time
python3 lane_emden_config_manager.py update --preset polytrope_n1_5_3d --end-time 10.0

# Update multiple parameters
python3 lane_emden_config_manager.py update --preset polytrope_n1_5_3d \
    --relax-steps 5000 \
    --end-time 8.0 \
    --particles 10000
```

### 3. Apply Preset to Generate C++ Config
```bash
# Generate lane_emden/lane_emden.json from preset
python3 lane_emden_config_manager.py apply --preset polytrope_n1_5_3d
```

### 4. Run Simulation
```bash
./build/sph lane_emden
```

### 5. Create New Preset
```bash
# Create n=3.0 polytrope (radiation-dominated)
python3 lane_emden_config_manager.py create \
    --name polytrope_n3_3d \
    --polytrope 3.0 \
    --dimension 3 \
    --relax-steps 1000

# Apply it
python3 lane_emden_config_manager.py apply --preset polytrope_n3_3d
```

### 6. Validate Preset
```bash
python3 lane_emden_config_manager.py validate --preset polytrope_n1_5_3d
```

## 📊 Example: Complete Workflow

```bash
# 1. Update preset for high-resolution relaxation
python3 lane_emden_config_manager.py update \
    --preset polytrope_n1_5_3d \
    --relax-steps 10000

# 2. Generate C++ config from preset
python3 lane_emden_config_manager.py apply --preset polytrope_n1_5_3d

# 3. Run simulation
./build/sph lane_emden > lane_emden/results/run.log 2>&1 &

# 4. Monitor progress
tail -f lane_emden/results/run.log
```

## 🚨 Important Rules

1. **NEVER** edit `lane_emden/lane_emden.json` directly
2. **ALWAYS** update presets in `lane_emden/config/presets/`
3. **ALWAYS** run `apply` command after updating a preset
4. Presets are version-controlled; generated configs are not

## 🎓 Available Polytropes

| n   | γ     | Physical Meaning                  | Use Case               |
|-----|-------|-----------------------------------|------------------------|
| 0.5 | 3.0   | Very stiff EOS                    | White dwarfs           |
| 1.0 | 2.0   | Isothermal limit                  | Molecular clouds       |
| 1.5 | 5/3   | Monatomic ideal gas               | Main sequence stars    |
| 2.0 | 1.5   | Intermediate                      | -                      |
| 3.0 | 4/3   | Radiation-dominated               | Massive stars          |
| 5.0 | 1.2   | Degenerate gas                    | Low-mass stars         |

## 🔍 Quick Reference

```bash
# List presets
python3 lane_emden_config_manager.py list

# Update preset
python3 lane_emden_config_manager.py update --preset <name> --relax-steps <N>

# Apply preset
python3 lane_emden_config_manager.py apply --preset <name>

# Create preset
python3 lane_emden_config_manager.py create --name <name> --polytrope <n> --dimension <2|3>

# Validate preset
python3 lane_emden_config_manager.py validate --preset <name>
```

## 📝 Notes

- The `apply` command converts rich preset format to C++-compatible format
- Presets support physics metadata, documentation, and validation
- Multiple presets can coexist; switch between them with `apply`
- Old workflow (direct JSON editing) is deprecated but still works for backward compatibility
