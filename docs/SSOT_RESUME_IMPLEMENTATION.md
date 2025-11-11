# SSOT Implementation Summary - Checkpoint Resume

## Changes Made

### 1. Simplified Resume Configuration (JSON)

**File**: `lane_emden/config/presets/polytrope_n1_5_3d_resume.json`

**Before** (67% duplication, 55 lines):
```json
{
  "physics": { /* duplicated from checkpoint */ },
  "initial_conditions": { /* ignored anyway */ },
  "numerical": { /* duplicated from checkpoint */ },
  "artificial_viscosity": { /* duplicated from checkpoint */ },
  "checkpoint": { "resumeFile": "..." },
  "simulation": { ... }
}
```

**After** (SSOT compliant, 18 lines - 67% reduction):
```json
{
  "name": "polytrope_n1_5_3d_resume",
  "description": "Resume... - SSOT: All physics/particle state loaded from checkpoint",
  "checkpoint": {
    "enabled": true,
    "resumeFile": "/path/to/snapshot_0009.csv",
    "saveInterval": 1000,
    ...
  },
  "simulation": {
    "end_time": 5.0,
    "output_time": 0.05,
    ...
  }
}
```

### 2. C++ Code Changes

#### A. `solver.cpp` - Make Physics Parameters Optional

**Modified**: `Solver::read_parameterfile()`

```cpp
// Check if we're in resume mode
bool is_resume_mode = input.get<bool>("checkpoint.autoResume", false) || 
                      !input.get<std::string>("checkpoint.resumeFile", "").empty();

if (is_resume_mode) {
    std::cout << "Resume mode detected - physics parameters will be loaded from checkpoint" << std::endl;
}

// Use defaults for physics parameters (will be overridden by checkpoint)
m_param->physics.gamma = input.get<real>("gamma", 1.6666666666666667);  // Default
m_param->physics.neighbor_number = input.get<int>("neighborNumber", 50); // Default
// ... other parameters with defaults ...
```

**Why**: Allows minimal config files without requiring gamma, G, etc. when resuming.

#### B. `output_manager.hpp/cpp` - Return Physics Metadata

**Modified**: `OutputManager::load_for_resume()` signature

```cpp
// Before:
bool load_for_resume(const std::string& filepath,
                     std::shared_ptr<Simulation> sim,
                     CheckpointMetadata& checkpoint_meta);

// After:
bool load_for_resume(const std::string& filepath,
                     std::shared_ptr<Simulation> sim,
                     CheckpointMetadata& checkpoint_meta,
                     OutputMetadata* output_meta = nullptr);  // NEW: returns physics
```

**Implementation**:
```cpp
// Populate checkpoint metadata
checkpoint_meta = metadata.checkpoint_data;

// Return output metadata if requested (for physics parameters - SSOT)
if (output_meta != nullptr) {
    *output_meta = metadata;
}
```

**Why**: Provides access to physics parameters stored in checkpoint for SSOT override.

#### C. `solver.cpp` - Override Parameters from Checkpoint

**Modified**: `Solver::initialize()`

```cpp
// Check if we should resume from checkpoint
CheckpointMetadata resume_metadata;
OutputMetadata checkpoint_physics;  // NEW: Will contain physics parameters
bool resumed = false;

if(m_resume_from_checkpoint && !m_checkpoint_file.empty()) {
    if(m_output_manager->load_for_resume(m_checkpoint_file, m_sim, 
                                         resume_metadata, &checkpoint_physics)) {
        resumed = true;
        
        // SSOT: Override physics parameters from checkpoint metadata
        std::cout << "Applying physics parameters from checkpoint (SSOT):" << std::endl;
        m_param->physics.gamma = checkpoint_physics.gamma;
        m_param->gravity.constant = checkpoint_physics.gravitational_constant;
        m_param->gravity.is_valid = checkpoint_physics.use_gravity;
        m_param->type = checkpoint_physics.sph_type;
        m_param->kernel = checkpoint_physics.kernel_type;
        m_param->physics.neighbor_number = checkpoint_physics.neighbor_number;
        m_param->av.use_balsara_switch = checkpoint_physics.use_balsara;
        m_param->av.use_time_dependent_av = checkpoint_physics.use_time_dependent_av;
        // ... log each parameter ...
    }
}
```

**Why**: Ensures checkpoint file is authoritative source for all physics parameters (SSOT).

### 3. Documentation

Created comprehensive docs explaining SSOT principle:

- `lane_emden/config/presets/SSOT_DESIGN.md` - Detailed design rationale
- `lane_emden/config/presets/RESUME_INSTRUCTIONS.md` - Updated with SSOT notes
- `lane_emden/RESUME_QUICK_START.md` - Quick reference with SSOT explanation

### 4. Makefile Integration

**File**: `lane_emden/Makefile.lane_emden`

Added `lane_emden_resume` target that works with minimal configs.

## Benefits

### 1. Eliminated Duplication
- **Before**: 55 lines of config (physics duplicated)
- **After**: 18 lines of config (only runtime control)
- **Reduction**: 67% smaller, zero duplication

### 2. Eliminated Conflict Risk
- **Before**: Config γ=1.5, Checkpoint γ=1.667 → Which wins?
- **After**: Only checkpoint has γ → No conflict possible

### 3. Clear Separation of Concerns
- **Checkpoint file**: State (SSOT for physics + particles)
- **Config file**: Control (how long to run, where to save)

### 4. Easier Maintenance
- Don't need to keep config in sync with checkpoint
- Can't accidentally use wrong physics parameters

## How It Works

```
┌─────────────────────────────────────────────────────────────┐
│ 1. Read minimal config file                                 │
│    - checkpoint.resumeFile = "snapshot_0009.csv"            │
│    - simulation.end_time = 5.0                              │
│    - Use DEFAULT values for gamma, G, etc.                  │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 2. Detect resume mode                                       │
│    - resumeFile is not empty                                │
│    - Set is_resume_mode = true                              │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 3. Load checkpoint                                           │
│    - Read particles from CSV/HDF5                           │
│    - Read metadata (physics params)                         │
│    - Return OutputMetadata with γ, G, SPH type, etc.       │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 4. Override parameters from checkpoint (SSOT)               │
│    - m_param->gamma = checkpoint_physics.gamma              │
│    - m_param->G = checkpoint_physics.gravitational_constant │
│    - m_param->type = checkpoint_physics.sph_type            │
│    - ... all physics from checkpoint ...                    │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 5. Continue simulation                                      │
│    - Use physics from checkpoint (authoritative)            │
│    - Use runtime control from config (end_time, etc.)       │
│    - Save results to output_directory from config           │
└─────────────────────────────────────────────────────────────┘
```

## Verification

### Console Output When Resuming

You should now see:
```
Resume mode detected - physics parameters will be loaded from checkpoint
...
=== Attempting to resume from checkpoint ===
Loading checkpoint from: snapshot_0009.csv
...
Applying physics parameters from checkpoint (SSOT):
  - gamma: 1.66667
  - G: 1 (gravity disabled)
  - SPH type: 0
  - Kernel: 1
  - Neighbor number: 50
  - Balsara switch: enabled
  - Time-dependent AV: enabled
=== Successfully resumed from checkpoint ===
```

### What Gets Overridden

From checkpoint metadata (SSOT):
- ✅ gamma (adiabatic index)
- ✅ G (gravitational constant)
- ✅ use_gravity flag
- ✅ SPH type (SSPH/DISPH/GSPH)
- ✅ Kernel type
- ✅ Neighbor number
- ✅ Balsara switch flag
- ✅ Time-dependent AV flag
- ✅ All particle data

From config file (runtime control):
- ✅ end_time
- ✅ output_time
- ✅ output_directory
- ✅ checkpoint settings
- ✅ resumeFile path

## Testing

1. **Build**: `cd build && make -j8` ✅ Success
2. **Resume**: `./build/sph lane_emden/config/presets/polytrope_n1_5_3d_resume.json`
3. **Verify**: Check console output shows "Applying physics parameters from checkpoint"

## Files Modified

- ✅ `src/solver.cpp` - Resume mode detection + SSOT parameter override
- ✅ `src/output_manager.cpp` - Return physics metadata
- ✅ `include/output_manager.hpp` - Updated signature
- ✅ `lane_emden/config/presets/polytrope_n1_5_3d_resume.json` - Minimal config
- ✅ `lane_emden/Makefile.lane_emden` - Fixed CONFIG_MANAGER path
- ✅ Documentation files (SSOT_DESIGN.md, RESUME_INSTRUCTIONS.md, etc.)

## Backward Compatibility

Old-style resume configs with physics parameters will still work:
- Config provides gamma=1.5
- Checkpoint provides gamma=1.667
- Result: **Checkpoint wins** (gamma=1.667 used)
- Console shows override message

This ensures SSOT even with legacy configs.
