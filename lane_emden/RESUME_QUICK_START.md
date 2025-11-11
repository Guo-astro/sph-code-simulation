# Quick Start: Resume Lane-Emden Simulation

## SSOT Design Principle

The resume configuration follows **Single Source of Truth (SSOT)**:
- ✅ **Checkpoint file** = source of truth for physics & particle state
- ✅ **Config file** = only runtime control (how long to run, where to save)
- ❌ **No duplication** = physics parameters NOT in resume config

All physics parameters (γ, G, K, masses, etc.) and particle data are loaded from the checkpoint file metadata.

## Using the Makefile Command (Easiest)

### Resume from your checkpoint file:

```bash
make -f lane_emden/Makefile.lane_emden lane_emden_resume \
  CHECKPOINT=lane_emden/results/polytrope_n1.5_3d/snapshot_0009.csv
```

### Alternative with short path (relative to preset's output directory):

```bash
make -f lane_emden/Makefile.lane_emden lane_emden_resume \
  PRESET=polytrope_n1_5_3d \
  CHECKPOINT=snapshot_0009.csv
```

### Resume from a checkpoint in the checkpoints subdirectory:

```bash
make -f lane_emden/Makefile.lane_emden lane_emden_resume \
  CHECKPOINT=lane_emden/results/polytrope_n1.5_3d/checkpoints/checkpoint_0010.csv
```

## What the Command Does

1. **Loads preset configuration** (default: `polytrope_n1_5_3d`)
2. **Configures checkpoint resume** with your specified file
3. **Runs the simulation** starting from the checkpoint state
4. **Continues calculation** until completion

## Examples

### Example 1: Resume from snapshot_0009.csv
```bash
cd /Users/guo/Downloads/sphcode
make -f lane_emden/Makefile.lane_emden lane_emden_resume \
  CHECKPOINT=lane_emden/results/polytrope_n1.5_3d/snapshot_0009.csv
```

### Example 2: Resume with different preset
```bash
make -f lane_emden/Makefile.lane_emden lane_emden_resume \
  PRESET=polytrope_n3_3d \
  CHECKPOINT=lane_emden/results/polytrope_n3_3d/snapshot_0015.csv
```

### Example 3: Resume using short filename (auto-detected from preset)
```bash
make -f lane_emden/Makefile.lane_emden lane_emden_resume \
  PRESET=polytrope_n1_5_3d \
  CHECKPOINT=snapshot_0009.csv
```

## Expected Output

When you run the command, you should see:

```
==========================================
Lane-Emden Resume from Checkpoint
==========================================
Preset: polytrope_n1_5_3d
Checkpoint file: lane_emden/results/polytrope_n1.5_3d/snapshot_0009.csv

Loading preset configuration...
Configuring resume from checkpoint...
✓ Configured to resume from: lane_emden/results/polytrope_n1.5_3d/snapshot_0009.csv
Output directory:
  📁 lane_emden/results/polytrope_n1.5_3d

Running simulation from checkpoint...
...
=== Attempting to resume from checkpoint ===
Loading checkpoint from: ...
Successfully loaded 5400 particles from checkpoint
Resume from step XXX, time XXX
=== Successfully resumed from checkpoint ===
...
```

## Troubleshooting

### Error: "Checkpoint file not found"
- Check that the path is correct
- Use `ls -lh <path>` to verify the file exists
- Try using the full path instead of relative path

### Error: "CHECKPOINT file not specified"
- You must provide `CHECKPOINT=<file>` parameter
- See examples above

### Simulation starts from beginning
- Make sure the checkpoint file path is correct
- Check that the file has valid metadata headers
- Verify the file is not corrupted

## All Lane-Emden Commands

To see all available commands:
```bash
make -f lane_emden/Makefile.lane_emden lane_emden_help
```

Key commands:
- `lane_emden` - Run simulation from scratch
- `lane_emden_resume` - Resume from checkpoint (NEW!)
- `lane_emden_relax_only` - Run only relaxation phase
- `lane_emden_animate` - Generate animation from results
- `lane_emden_list` - List all presets
- `lane_emden_help` - Show help
