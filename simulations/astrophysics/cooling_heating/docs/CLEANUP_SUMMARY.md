# Folder Structure Cleanup Summary

**Date:** December 1, 2024  
**Purpose:** Reorganize cooling_heating test to match standard codebase structure

---

## Problem Statement

The `cooling_heating/` directory had accumulated significant clutter:
- **17 markdown files scattered in root** (should be in `docs/`)
- **2 Makefiles**: `Makefile.cooling` and `Makefile.cooling_heating` (redundant)
- **Duplicate config files**: `cnm_relaxation.json` in root AND `config/`
- **No consistent structure** compared to other test cases (shock_tube, sedov)

---

## Standard Structure Pattern

Based on existing test cases (`shock_tube/`, `sedov/`), the standard pattern is:

```
sample/<test>/
├── Makefile.<test>      # SINGLE unified Makefile
├── README.md            # Main documentation
├── <test>.json         # Active config (auto-generated)
├── config/
│   └── presets/        # Configuration presets
│       └── *.json
├── scripts/            # Python scripts
│   └── *.py
├── results/            # Simulation outputs
│   └── <run_name>/
└── docs/               # Documentation
    └── *.md
```

---

## Actions Taken

### 1. Created Subdirectories
```bash
mkdir -p docs/ config/presets/
```

### 2. Moved Documentation (17 files)
```bash
mv *.md docs/
```

Moved files:
- ANIMATION_GUIDE.md
- CNM_EQUILIBRIUM_PHYSICS.md
- CNM_TEST_NEXT_STEPS.md
- CNM_TEST_SUMMARY.md
- COOLING_HEATING_SUMMARY.md
- ISM_COOLING_HEATING_IMPLEMENTATION.md
- ISM_INTEGRATION_SUMMARY.md
- KOYAMA_INUTSUKA_COOLING_IMPLEMENTATION.md
- LOW_SNAPSHOT_COUNT_DIAGNOSIS.md
- NEXT_STEPS.md
- PARAMETER_MAPPING.md
- README.md
- TEST_CASE_SUMMARY.md
- THERMAL_CFL_FIX.md (CRITICAL - thermal timescale fix)
- THERMAL_TIMESCALE_IMPLEMENTATION.md
- And 2 others

### 3. Moved Config Presets
```bash
mv cnm_relaxation.json config/presets/
rm -f config/ism_cooling_1d.json  # duplicate removed
```

### 4. Merged Makefiles

**Before:**
- `Makefile.cooling` (169 lines) - CNM test targets
- `Makefile.cooling_heating` (201 lines) - Python solver targets

**After:**
- `Makefile.cooling_heating` (305 lines) - **UNIFIED**
  - Contains ALL targets from both Makefiles
  - Fixed duplicate `cooling_clean` target (merged)
  - CNM targets: `cnm_run`, `cnm_plot`, `cnm_animate`, `cnm_all`, `cnm_clean`
  - Python targets: `cooling_reproduce`, `cooling_run`, `cooling_compare`

```bash
rm -f Makefile.cooling  # removed redundant file
```

---

## Final Structure

```
simulations/astrophysics/cooling_heating/
├── Makefile.cooling_heating    ← SINGLE unified Makefile ✓
├── ism_cooling_1d.json        ← Active config (auto-generated)
├── config/
│   └── presets/
│       └── cnm_relaxation.json ← CNM preset ✓
├── scripts/                    ← Python scripts
│   ├── animate_cnm_relaxation.py
│   ├── plot_cnm_relaxation.py
│   └── (13 other chemistry scripts)
├── results/                    ← Simulation outputs
│   ├── cnm_relaxation/
│   └── ism_cooling_1d/
└── docs/                       ← ALL documentation ✓
    └── (15 markdown files)
```

---

## Verification

### No Makefile Warnings
```bash
$ make -f simulations/astrophysics/cooling_heating/Makefile.cooling_heating cnm_help
# ✓ No warnings (duplicate cooling_clean fixed)
```

### Structure Matches Standard
```
STANDARD (shock_tube):          YOUR (cooling_heating):
├── Makefile.shock_tube         ├── Makefile.cooling_heating ✓
├── README.md                   ├── docs/                     ✓
├── config/presets/             ├── config/presets/           ✓
├── scripts/                    ├── scripts/                  ✓
└── results/                    └── results/                  ✓
```

---

## Key Benefits

1. **Professional Organization**: Matches codebase conventions
2. **Single Source of Truth**: ONE Makefile, no redundancy
3. **Clean Root**: Only essential files visible
4. **Easy Navigation**: Logical subdirectory structure
5. **Maintainability**: Clear separation of docs/configs/scripts/results

---

## Related Documentation

- **THERMAL_CFL_FIX.md** - Critical thermal timescale fix (dtMax: 1.0 → 0.01)
- **CNM_TEST_SUMMARY.md** - Physics and setup for K&I (2000) Model C10
- **ANIMATION_GUIDE.md** - How to use animation scripts

---

## Usage

```bash
# CNM thermal relaxation test
make -f simulations/astrophysics/cooling_heating/Makefile.cooling_heating cnm_help
make -f simulations/astrophysics/cooling_heating/Makefile.cooling_heating cnm_all

# Python cooling solver reproduction
make -f simulations/astrophysics/cooling_heating/Makefile.cooling_heating cooling_help
make -f simulations/astrophysics/cooling_heating/Makefile.cooling_heating cooling_reproduce
```
