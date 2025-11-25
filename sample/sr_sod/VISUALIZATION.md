# SR-GSPH Sod Shock Tube Visualization

## Single Source of Truth (SSOT)

All visualization is now handled by a single script: **`scripts/sr_sod_visualize.py`**

### Usage

```bash
# From sample/sr_sod directory:
python scripts/sr_sod_visualize.py results/<test_name> [options]

# Examples:
python scripts/sr_sod_visualize.py results/sod_kitajima --all
python scripts/sr_sod_visualize.py results/blast_wave --animate --format mp4
python scripts/sr_sod_visualize.py results/strong_blast --plot
```

### Options

| Option | Description |
|--------|-------------|
| `--animate, -a` | Create animation |
| `--plot, -p` | Create static plots (final state + evolution) |
| `--all` | Create all visualizations (default) |
| `--format {gif,mp4}` | Animation format (default: gif) |
| `--fps FPS` | Frames per second (default: 12) |

### Output Files

- `sr_sod_animation_<solver>.gif` or `.mp4` - Time evolution animation
- `sr_sod_final_state.png` - Final state plot
- `sr_sod_evolution.png` - Multi-time overlay plot

### Features

- **Physical Units**: All axes labeled with proper units (c, L/c, n₀, P₀)
- **Time Info Box**: Shows current time `t`, output timestep `Δt_out`, and frame number
- **Auto-Detection**: Automatically detects test type and solver from config/directory
- **Consistent Styling**: Uniform appearance across all visualizations

### Physical Units (Code Units)

In SR-GSPH simulations:
- Speed of light: c = 1
- Time: measured in L/c (light-crossing time)
- Velocity: measured in units of c
- Density: rest-frame baryon number density n
- Pressure: in code units P₀

### Deprecated Scripts

Old scripts are archived in `scripts/deprecated/`. See migration guide there.
