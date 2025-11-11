# How to Resume from snapshot_0009.csv

## Quick Start

To resume your simulation from `snapshot_0009.csv`, use the provided configuration file:

```bash
cd /Users/guo/Downloads/sphcode
./build/sph lane_emden/config/presets/polytrope_n1_5_3d_resume.json
```

## What Will Happen

1. **The simulation will load** all particle data and physics parameters from `snapshot_0009.csv`
2. **Resume from that state** preserving:
   - Particle positions, velocities, densities
   - All physics parameters (γ, G, K, masses, etc.)
   - Simulation time and step counter
   - All state variables
3. **Continue calculation** until completion or the next checkpoint

## SSOT Principle (Single Source of Truth)

The resume configuration is **minimal by design**:
- **Physics parameters** come from the checkpoint file (SSOT)
- **Particle data** comes from the checkpoint file (SSOT)
- **Only control parameters** are specified (where to save, how long to run)

This ensures no conflicts between config file and checkpoint state.

## Key Configuration Settings

The resume configuration (`polytrope_n1_5_3d_resume.json`) is minimal - following SSOT:

```json
{
  "checkpoint": {
    "autoResume": false,
    "resumeFile": "/path/to/snapshot_0009.csv"
  },
  "simulation": {
    "end_time": 5.0,
    "output_time": 0.05,
    "output_directory": "lane_emden/results/polytrope_n1.5_3d"
  }
}
```

**What's NOT included (loaded from checkpoint instead):**
- Physics parameters (γ, G, K, masses) - from checkpoint metadata
- Particle data (positions, velocities, densities) - from checkpoint
- Numerical settings (kernel, SPH type, neighbors) - from checkpoint
- Initial conditions - irrelevant when resuming

**What IS included (runtime control):**
- Checkpoint configuration (where to resume from)
- Simulation control (how long to run, output settings)
- Output directory (where to save new results)

## Alternative: Using Relative Path

If you prefer a relative path, edit the config to use:

```json
"resumeFile": "lane_emden/results/polytrope_n1.5_3d/snapshot_0009.csv"
```

## Verification

When you run the simulation, you should see:

```
=== Attempting to resume from checkpoint ===
Loading checkpoint from: /Users/guo/Downloads/sphcode/lane_emden/results/polytrope_n1.5_3d/snapshot_0009.csv
Successfully loaded XXXX particles from checkpoint
Resume from step XXX, time XXX
=== Successfully resumed from checkpoint ===
```

## Important Notes

1. **SSOT Principle** - All physics parameters come from the checkpoint file, not the config
2. **Minimal Configuration** - Only specify runtime control (end time, output settings, checkpoint settings)
3. **Checkpoint is Authoritative** - The checkpoint file contains the complete state; config just controls what to do with it
4. **File format** - The code supports `.csv`, `.h5`, and `.vtk` checkpoint files

## Troubleshooting

### "Failed to load checkpoint"
- Check that the file path is correct
- Verify the file exists and is readable
- Ensure the CSV file has valid metadata headers

### "Unsupported checkpoint file format"
- Make sure the file extension is `.csv`, `.h5`, or `.vtk`

### Simulation starts from beginning
- Check that `resumeFile` is not empty
- Verify `autoResume` is set correctly
- Look for error messages in the console output

## Customizing the Resume

You can modify `polytrope_n1_5_3d_resume.json` to:

- **Change output directory**: Modify `simulation.output_directory`
- **Adjust checkpoint frequency**: Change `checkpoint.saveInterval`
- **Extend simulation**: Increase `simulation.end_time`
- **Change output frequency**: Modify `simulation.output_time`

## Example: Resume and Extend Simulation

To resume and run longer:

1. Edit `polytrope_n1_5_3d_resume.json`
2. Change `simulation.end_time` to desired value (e.g., `10.0`)
3. Run: `./build/sph lane_emden/config/presets/polytrope_n1_5_3d_resume.json`

The simulation will continue from snapshot_0009 until the new end time.
