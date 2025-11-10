# Checkpoint System Documentation

## Overview

The SPH simulation checkpoint system allows you to save and resume simulations at any point, particularly useful for long-running relaxation processes. This system saves the complete particle state, simulation time, and metadata, enabling seamless continuation after interruptions.

## Key Features

- ✅ **Auto-save during relaxation** - Save checkpoints at regular intervals during relaxation
- ✅ **Resume from any checkpoint** - Continue relaxation or simulation from where you left off
- ✅ **Minimal overhead** - Binary format for fast I/O
- ✅ **Metadata tracking** - JSON metadata files for human readability
- ✅ **Automatic cleanup** - Keep only N most recent checkpoints
- ✅ **Works for both relaxation and evolution** - Unified checkpoint format

## Quick Start

### 1. Enable Checkpoints in Configuration

Add a `checkpoint` section to your JSON configuration file:

```json
{
  "checkpoint": {
    "enabled": true,
    "saveInterval": 500,
    "directory": "checkpoints",
    "saveOnExit": true,
    "maxCheckpoints": 5,
    "autoResume": false,
    "resumeFile": ""
  }
}
```

### 2. Run Simulation with Checkpoints

```bash
# Start a new relaxation with checkpoints enabled
./build/sph lane_emden
```

The simulation will automatically save checkpoints to the specified directory every 500 steps.

### 3. Resume from Checkpoint

**Option A: Auto-resume from latest checkpoint**
```json
{
  "checkpoint": {
    "enabled": true,
    "autoResume": true
  }
}
```

**Option B: Resume from specific checkpoint**
```json
{
  "checkpoint": {
    "enabled": true,
    "resumeFile": "checkpoints/checkpoint_002500.chk"
  }
}
```

Then run:
```bash
./build/sph lane_emden
```

## Configuration Options

### `checkpoint` Section

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `enabled` | bool | `false` | Enable/disable checkpoint system |
| `saveInterval` | int | `100` | Save checkpoint every N steps/iterations |
| `directory` | string | `"checkpoints"` | Directory to store checkpoint files |
| `saveOnExit` | bool | `true` | Always save checkpoint when simulation ends |
| `maxCheckpoints` | int | `5` | Maximum number of checkpoints to keep (0 = unlimited) |
| `autoResume` | bool | `false` | Automatically resume from latest checkpoint if available |
| `resumeFile` | string | `""` | Specific checkpoint file to resume from |

## Checkpoint File Format

### Files Generated

For each checkpoint, two files are created:

1. **Binary checkpoint file** (`.chk`)
   - Contains complete particle data (positions, velocities, densities, etc.)
   - Optimized for fast I/O
   - Example: `checkpoint_002500.chk`

2. **Metadata file** (`.json`)
   - Human-readable checkpoint information
   - Contains step number, time, particle count, etc.
   - Example: `checkpoint_002500.json`

### Metadata Example

```json
{
  "version": 1,
  "time": 0.245,
  "step": 0,
  "particle_num": 5400,
  "is_relaxation": true,
  "relaxation_step": 2500,
  "relaxation_total_steps": 30000,
  "accumulated_time": 0.245,
  "config_hash": "1.6666666666666667_50_1",
  "preset_name": "lane_emden"
}
```

## Use Cases

### Use Case 1: Long Relaxation Runs

**Problem**: Relaxation with 30,000 steps takes hours and may be interrupted.

**Solution**: Enable checkpoints with moderate save interval

```json
{
  "relaxation": {
    "enabled": true,
    "steps": 30000
  },
  "checkpoint": {
    "enabled": true,
    "saveInterval": 1000,
    "maxCheckpoints": 10
  }
}
```

**Workflow**:
1. Start relaxation: `./build/sph lane_emden`
2. If interrupted (Ctrl+C, system crash, etc.), checkpoints are saved
3. Resume by setting `"autoResume": true` and run again
4. Relaxation continues from last checkpoint

### Use Case 2: Iterative Relaxation

**Problem**: You want to check relaxation progress and decide whether to continue.

**Solution**: Run relaxation in batches with checkpoint resume

```json
{
  "relaxation": {
    "steps": 5000
  },
  "checkpoint": {
    "enabled": true,
    "saveOnExit": true,
    "autoResume": true
  }
}
```

**Workflow**:
1. Run 5000 steps: `./build/sph lane_emden`
2. Analyze results
3. If not converged, run again - automatically resumes from step 5000
4. Repeat until satisfied

### Use Case 3: Testing Different Relaxation Parameters

**Problem**: Want to try different relaxation parameters from the same initial state.

**Solution**: Save checkpoint at start, resume from it with different configs

```bash
# Run 1000 steps to generate checkpoint
./build/sph lane_emden  # with checkpoint.enabled=true

# Copy checkpoint
cp checkpoints/checkpoint_001000.chk my_base_state.chk

# Try different damping factors by resuming from same checkpoint
# In config: "resumeFile": "my_base_state.chk"
```

### Use Case 4: Production Runs with Safety

**Problem**: Long evolution simulation that must not lose progress.

**Solution**: Frequent checkpoints with automatic cleanup

```json
{
  "checkpoint": {
    "enabled": true,
    "saveInterval": 50,
    "maxCheckpoints": 20,
    "saveOnExit": true,
    "autoResume": true
  }
}
```

## Command-Line Examples

### Example 1: Start Fresh Relaxation with Checkpoints

```bash
# Edit lane_emden.json to enable checkpoints
./build/sph lane_emden
```

Output:
```
Checkpoint system enabled:
  - Save interval: every 500 steps
  - Directory: checkpoints
  - Max checkpoints: 5
  
=== Starting Relaxation Phase (30000 steps) ===
Progress: [=====...] 10% | Step 3000/30000
✓ Checkpoint saved: checkpoint_003000.chk
  - Particles: 5400
  - Time: 0.156
  - Relaxation step: 3000/30000
```

### Example 2: Resume After Interruption

```bash
# Set autoResume: true in config
./build/sph lane_emden
```

Output:
```
=== Attempting to resume from checkpoint ===
✓ Checkpoint loaded: checkpoints/checkpoint_015000.chk
  - Particles: 5400
  - Time: 0.782
  - Resuming relaxation at step: 15000/30000
=== Successfully resumed from checkpoint ===

=== Starting Relaxation Phase (30000 steps) ===
Resuming from step 15000 (time=0.782)
Progress: [========...] 50% | Step 15000/30000
```

### Example 3: Resume from Specific Checkpoint

```bash
# In JSON config:
# "checkpoint": {
#   "resumeFile": "checkpoints/checkpoint_010000.chk"
# }

./build/sph lane_emden
```

## Checkpoint Management

### List Available Checkpoints

```bash
ls -lh checkpoints/
```

Output:
```
checkpoint_001000.chk    5.2M
checkpoint_001000.json   320B
checkpoint_002000.chk    5.2M
checkpoint_002000.json   320B
...
```

### View Checkpoint Metadata

```bash
cat checkpoints/checkpoint_002500.json
```

### Delete Old Checkpoints Manually

```bash
# Keep only last 3 checkpoints
cd checkpoints
ls -t checkpoint_*.chk | tail -n +4 | xargs rm
ls -t checkpoint_*.json | tail -n +4 | xargs rm
```

### Copy Checkpoint for Backup

```bash
# Backup checkpoint before critical change
cp checkpoints/checkpoint_025000.chk ~/backup/critical_state.chk
cp checkpoints/checkpoint_025000.json ~/backup/critical_state.json
```

## Troubleshooting

### Q: Checkpoint not loading

**A**: Check that:
- Checkpoint files exist in the specified directory
- `checkpoint.enabled` is `true`
- Either `autoResume` is `true` or `resumeFile` is specified
- File paths are correct (absolute or relative to run directory)

### Q: "Invalid checkpoint file format" error

**A**: Checkpoint file is corrupted or from incompatible version. Start fresh.

### Q: Simulation starts from beginning despite checkpoint

**A**: Verify:
```bash
# Check if autoResume is enabled
grep autoResume lane_emden.json

# Check if checkpoint exists
ls -lh checkpoints/checkpoint_*.chk | tail -5
```

### Q: Too many checkpoint files filling disk

**A**: Reduce `maxCheckpoints` or increase `saveInterval`:
```json
{
  "checkpoint": {
    "maxCheckpoints": 3,
    "saveInterval": 2000
  }
}
```

### Q: Want to restart from scratch but autoResume is on

**A**: Either:
1. Set `"autoResume": false` in config, or
2. Delete checkpoint directory: `rm -rf checkpoints/`

## Technical Details

### What's Saved in Checkpoints?

For each particle:
- Position (`pos[3]`)
- Velocity (`vel[3]`)
- Acceleration (`acc[3]`)
- Mass (`mass`)
- Density (`dens`)
- Pressure (`pres`)
- Internal energy (`ene`)
- Smoothing length (`smth`)
- Sound speed (`sound`)
- Artificial viscosity (`alpha`)
- Balsara switch (`balsara`)
- Particle ID (`id`)

Global state:
- Simulation time
- Step/iteration count
- Relaxation progress (if applicable)

### File Size

Approximate checkpoint file sizes:
- **5,400 particles**: ~5 MB per checkpoint
- **50,000 particles**: ~50 MB per checkpoint
- **500,000 particles**: ~500 MB per checkpoint

### Performance Impact

- **Saving**: ~0.1-0.5 seconds for typical runs (5K-10K particles)
- **Loading**: ~0.1-0.3 seconds
- **Overhead**: Negligible if `saveInterval` is reasonable (≥100 steps)

### Compatibility

Checkpoints are compatible across:
- ✅ Same compiler (same architecture)
- ✅ Same `DIM` setting (2D vs 3D)
- ✅ Same `real` type (float vs double)

Checkpoints are **NOT** compatible across:
- ❌ Different architectures (x86 vs ARM)
- ❌ Different `DIM` settings
- ❌ Different `real` precision

## Integration with Lane-Emden Workflow

### Preset Configuration

Create a preset with checkpoint enabled:

```json
{
  "name": "polytrope_n1_5_long_relax",
  "relaxation": {
    "enabled": true,
    "steps": 50000
  },
  "checkpoint": {
    "enabled": true,
    "saveInterval": 2500,
    "directory": "lane_emden/results/polytrope_n1.5_3d/checkpoints",
    "autoResume": true,
    "maxCheckpoints": 10
  }
}
```

### Makefile Integration

You can add Makefile targets for checkpoint management:

```makefile
.PHONY: lane_emden_checkpoint lane_emden_resume

lane_emden_checkpoint:
	@echo "Running Lane-Emden with checkpoints..."
	./build/sph lane_emden

lane_emden_resume:
	@echo "Resuming Lane-Emden from checkpoint..."
	# Ensure autoResume is enabled in config
	./build/sph lane_emden

lane_emden_clean_checkpoints:
	@echo "Cleaning checkpoints..."
	rm -rf lane_emden/results/*/checkpoints/
```

## Best Practices

1. **Choose appropriate save interval**
   - Too frequent: Performance impact, disk space waste
   - Too infrequent: Lose too much progress on failure
   - Recommended: 500-2000 steps for relaxation

2. **Set max_checkpoints wisely**
   - Keep 3-5 recent checkpoints for safety
   - More if testing/debugging
   - Consider disk space vs. safety tradeoff

3. **Use descriptive directories**
   - Include run name: `"directory": "runs/run_2025_01_10/checkpoints"`
   - Separate checkpoints by configuration

4. **Backup critical checkpoints**
   - Copy checkpoints before major parameter changes
   - Archive successful relaxation states

5. **Monitor checkpoint files**
   - Periodically check checkpoint directory size
   - Verify checkpoints are being created
   - Test resume capability early in long runs

## Future Enhancements

Planned features:
- [ ] Compression support (gzip/lz4)
- [ ] Parallel I/O for large particle counts
- [ ] Checkpoint validation/checksum
- [ ] Cloud storage integration
- [ ] Automatic checkpoint on signal (SIGTERM/SIGINT)
- [ ] Restart with different parameters (adaptive)

## See Also

- [RELAXATION_ONLY_MODE.md](RELAXATION_ONLY_MODE.md) - Relaxation-only mode
- [WORKFLOW.md](../WORKFLOW.md) - Lane-Emden workflow
- [LANE_EMDEN_ARCHITECTURE.md](LANE_EMDEN_ARCHITECTURE.md) - Architecture overview

---

**Created**: 2025-11-10  
**Last Updated**: 2025-11-10  
**Author**: SPH Development Team
