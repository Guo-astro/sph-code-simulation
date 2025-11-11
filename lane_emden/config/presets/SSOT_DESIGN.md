# SSOT Design for Resume Configuration

## Problem: Redundant Configuration (Anti-Pattern)

**BAD - Violates DRY and SSOT:**
```json
{
  "physics": {
    "polytropic_index": 1.5,        // ❌ Duplicated from checkpoint
    "adiabatic_index": 1.666...,    // ❌ Duplicated from checkpoint
    "gravitational_constant": 1.0    // ❌ Duplicated from checkpoint
  },
  "initial_conditions": { ... },     // ❌ Ignored anyway
  "numerical": { ... },               // ❌ Duplicated from checkpoint
  "checkpoint": {
    "resumeFile": "snapshot.csv"     // ✅ This is the ONLY source of truth
  }
}
```

**Issues:**
1. **Duplication** - Same data in two places (config + checkpoint)
2. **Conflict Risk** - What if config values differ from checkpoint?
3. **Maintenance** - Must keep config in sync with checkpoint
4. **Confusion** - Which is authoritative?

## Solution: Minimal Resume Config (SSOT)

**GOOD - Follows SSOT:**
```json
{
  "name": "resume_config",
  "description": "Minimal resume - SSOT: checkpoint is authoritative",
  "checkpoint": {
    "enabled": true,
    "resumeFile": "/path/to/snapshot.csv",  // ✅ Source of truth
    "saveInterval": 1000,                    // ✅ Runtime control
    "directory": "results/checkpoints"       // ✅ Runtime control
  },
  "simulation": {
    "end_time": 10.0,                        // ✅ Runtime control
    "output_time": 0.1,                      // ✅ Runtime control
    "output_directory": "results"            // ✅ Runtime control
  }
}
```

**Benefits:**
1. ✅ **Single Source of Truth** - Checkpoint file is authoritative
2. ✅ **No Duplication** - Physics/particle data only in checkpoint
3. ✅ **No Conflicts** - Impossible to have mismatched parameters
4. ✅ **Clarity** - Clear separation: checkpoint=state, config=control

## What Comes From Where

### From Checkpoint File (SSOT for State)
- ✅ All particle data (position, velocity, density, energy, etc.)
- ✅ Physics parameters (γ, G, K, polytropic index, etc.)
- ✅ Numerical settings (kernel type, SPH type, neighbor number)
- ✅ Simulation state (time, step number)
- ✅ Lane-Emden relaxation parameters (α, ρ_c, R, M)

### From Resume Config (Runtime Control Only)
- ✅ Where to resume from (`resumeFile`)
- ✅ How long to continue (`end_time`)
- ✅ Output frequency (`output_time`)
- ✅ Where to save results (`output_directory`)
- ✅ Checkpoint settings (`saveInterval`, `maxCheckpoints`)

## Implementation in Code

The `load_for_resume()` function in `output_manager.cpp` reads:
```cpp
// Read all state from checkpoint file
metadata = read_metadata_from_checkpoint(filepath);

// Restore everything from checkpoint
sim->particles = metadata.particles;          // Particle data
sim->time = metadata.time;                    // Time
params->gamma = metadata.gamma;               // Physics
params->G = metadata.gravitational_constant;  // Physics
// ... all other state from checkpoint ...
```

The config file only controls:
```cpp
// From config: runtime behavior
end_time = config["simulation"]["end_time"];
output_interval = config["simulation"]["output_time"];
checkpoint_interval = config["checkpoint"]["saveInterval"];
```

## Real-World Analogy

**Resume from checkpoint** is like **loading a saved game**:
- 🎮 **Save file** = Complete game state (SSOT)
- ⚙️ **Settings menu** = Controls (graphics, volume, difficulty)
- ❌ **Don't duplicate** game state in settings menu

Same principle here:
- 📁 **Checkpoint file** = Complete simulation state (SSOT)
- ⚙️ **Resume config** = Runtime controls (how long, where to save)
- ❌ **Don't duplicate** physics/state in config

## Migration from Old Style

If you have old resume configs with physics parameters:

**Before (verbose, violates SSOT):**
```json
{
  "physics": { /* 50 lines */ },
  "initial_conditions": { /* 30 lines */ },
  "numerical": { /* 40 lines */ },
  "checkpoint": { "resumeFile": "snapshot.csv" }
}
```

**After (minimal, follows SSOT):**
```json
{
  "checkpoint": { "resumeFile": "snapshot.csv" },
  "simulation": { "end_time": 10.0, "output_time": 0.1 }
}
```

**Result:** 90% reduction in config size, zero duplication, zero conflict risk.

## Validation

To verify SSOT compliance, check:
1. ✅ Resume config has NO physics parameters
2. ✅ Resume config has NO particle data
3. ✅ Resume config has NO initial conditions
4. ✅ Resume config ONLY has checkpoint path + runtime control
5. ✅ All state comes from checkpoint file

## References

- **DRY**: Don't Repeat Yourself - Every piece of knowledge must have a single, unambiguous, authoritative representation
- **SSOT**: Single Source of Truth - Data should exist in exactly one place
- **Separation of Concerns**: State (checkpoint) vs. Control (config)
