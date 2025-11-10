# Checkpoint System Quick Reference

## Basic Configuration

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

## Common Workflows

### 1. Enable Checkpoints (First Time)

```bash
# Edit your JSON config, add checkpoint section
# Then run:
./build/sph lane_emden
```

### 2. Resume from Latest Checkpoint

```json
"checkpoint": {
  "enabled": true,
  "autoResume": true
}
```

```bash
./build/sph lane_emden
```

### 3. Resume from Specific Checkpoint

```json
"checkpoint": {
  "enabled": true,
  "resumeFile": "checkpoints/checkpoint_010000.chk"
}
```

```bash
./build/sph lane_emden
```

### 4. Start Fresh (Ignore Checkpoints)

```json
"checkpoint": {
  "enabled": true,
  "autoResume": false,
  "resumeFile": ""
}
```

Or delete checkpoint directory:
```bash
rm -rf checkpoints/
./build/sph lane_emden
```

## Configuration Parameters

| Parameter | Type | Default | What It Does |
|-----------|------|---------|--------------|
| `enabled` | bool | `false` | Turn checkpoint system on/off |
| `saveInterval` | int | `100` | Save every N steps |
| `directory` | string | `"checkpoints"` | Where to save checkpoints |
| `saveOnExit` | bool | `true` | Save when simulation finishes |
| `maxCheckpoints` | int | `5` | Keep only N most recent (0=all) |
| `autoResume` | bool | `false` | Auto-load latest checkpoint |
| `resumeFile` | string | `""` | Load specific checkpoint |

## Checkpoint Files

### What Gets Created

```
checkpoints/
├── checkpoint_000500.chk     # Binary particle data
├── checkpoint_000500.json    # Metadata (human-readable)
├── checkpoint_001000.chk
├── checkpoint_001000.json
└── ...
```

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
  "accumulated_time": 0.245
}
```

## Typical Use Cases

### Long Relaxation (30K+ steps)

```json
{
  "relaxation": {
    "steps": 30000
  },
  "checkpoint": {
    "enabled": true,
    "saveInterval": 1000,
    "maxCheckpoints": 10,
    "autoResume": true
  }
}
```

**Benefit**: Can interrupt and resume anytime without losing progress

### Testing/Development

```json
{
  "checkpoint": {
    "enabled": true,
    "saveInterval": 100,
    "maxCheckpoints": 3,
    "autoResume": false
  }
}
```

**Benefit**: Quick checkpoints for debugging, manual resume control

### Production Runs

```json
{
  "checkpoint": {
    "enabled": true,
    "saveInterval": 500,
    "maxCheckpoints": 20,
    "saveOnExit": true,
    "autoResume": true
  }
}
```

**Benefit**: Safety net for long runs, automatic recovery

## Command Examples

### View Checkpoints

```bash
ls -lh checkpoints/
cat checkpoints/checkpoint_002500.json
```

### Resume After Interrupt

```bash
# Method 1: Edit config to set autoResume=true
./build/sph lane_emden

# Method 2: Set resumeFile to specific checkpoint
./build/sph lane_emden
```

### Clean Up Old Checkpoints

```bash
# Automatic (via maxCheckpoints setting in config)
# Or manual:
rm checkpoints/checkpoint_00*.chk
rm checkpoints/checkpoint_00*.json
```

### Backup Important Checkpoint

```bash
cp checkpoints/checkpoint_020000.chk ~/backup/good_state.chk
cp checkpoints/checkpoint_020000.json ~/backup/good_state.json
```

## Troubleshooting

| Problem | Solution |
|---------|----------|
| Checkpoints not being created | Check `enabled: true` and write permissions |
| Resume not working | Verify `autoResume: true` or `resumeFile` is set |
| "Invalid checkpoint" error | Checkpoint corrupted, delete and restart |
| Too many files | Reduce `maxCheckpoints` or increase `saveInterval` |
| Disk space full | Lower `maxCheckpoints`, increase `saveInterval` |
| Starts from beginning | Check `autoResume` and checkpoint directory |

## Integration with Lane-Emden

### Example Preset with Checkpoints

`lane_emden/config/presets/polytrope_n1_5_checkpoint.json`:

```json
{
  "name": "polytrope_n1_5_3d_checkpoint",
  "relaxation": {
    "enabled": true,
    "steps": 30000
  },
  "checkpoint": {
    "enabled": true,
    "saveInterval": 1000,
    "directory": "lane_emden/results/polytrope_n1.5_3d/checkpoints",
    "autoResume": true,
    "maxCheckpoints": 10
  }
}
```

### Apply and Run

```bash
# Apply preset
python3 lane_emden_config_manager.py apply --preset polytrope_n1_5_checkpoint

# Run with checkpoints
./build/sph lane_emden

# If interrupted, just run again - will auto-resume!
./build/sph lane_emden
```

## Performance Notes

- **Overhead**: ~0.2 seconds per checkpoint (5K particles)
- **File size**: ~1 MB per 1K particles
- **Recommended interval**: 500-2000 steps for relaxation

## What's Saved?

Per particle:
- Position, velocity, acceleration
- Mass, density, pressure, energy
- Smoothing length, sound speed
- Artificial viscosity parameters
- Particle ID

Global state:
- Simulation time
- Step/iteration count
- Relaxation progress

## Quick Tips

1. **Start with defaults** - The default settings work well for most cases
2. **Use autoResume for long runs** - No need to manually find latest checkpoint
3. **Backup critical states** - Copy checkpoint before risky changes
4. **Monitor disk space** - Checkpoints can accumulate quickly
5. **Test resume early** - Verify checkpoint/resume works before long run

---

For full details, see [CHECKPOINT_SYSTEM.md](CHECKPOINT_SYSTEM.md)
