# 2D Neutron Star Merger Sample

Binary neutron star collision simulation in 2D using Special Relativistic Godunov SPH (SR-GSPH).

Based on the setup described in **Shibata & Hotokezaka (2019)**.

## Physics

- **Equation of State**: Γ=2 polytrope (P = K·ρ^γ), approximating nuclear matter
- **Density Profile**: Lane-Emden n=1 (corresponding to γ=2) gives the NS structure
- **Initial Setup**: Two identical neutron stars placed symmetrically, moving toward each other
- **Relativistic Treatment**: Uses SR-GSPH with Gaussian kernel

## Quick Start

```bash
# Build the code (from repository root)
make -j4 DIM=2

# Run quick test simulation (~3 seconds)
make -f simulations/astrophysics/ns_merger/Makefile.ns_merger ns_merger_quick

# Create animation
make -f simulations/astrophysics/ns_merger/Makefile.ns_merger ns_merger_animate

# Or do both at once
make -f simulations/astrophysics/ns_merger/Makefile.ns_merger ns_merger_all
```

## Available Targets

Run `make -f simulations/astrophysics/ns_merger/Makefile.ns_merger ns_merger_help` for full target list.

| Target | Description |
|--------|-------------|
| `ns_merger_quick` | Quick test (~1300 particles, 5 time units) |
| `ns_merger_run` | Full simulation (~5400 particles, 50 time units) |
| `ns_merger_animate` | Create animation GIF (quick test) |
| `ns_merger_grid` | Create evolution grid plot |
| `ns_merger_all` | Run quick test + all visualizations |
| `ns_merger_full_all` | Run full simulation + animation |

## Configuration

Two preset configurations are provided:

### Quick Test (`ns_merger_quick_test.json`)
- Particles: ~1300
- Star radius: 1.0
- Separation: 4.0
- Collision velocity: 0.1c
- End time: 5
- Computation time: ~3 seconds

### Full Simulation (`ns_merger_2d.json`)
- Particles: ~5400
- Star radius: 1.2
- Separation: 6.0
- Collision velocity: 0.15c
- End time: 50
- Computation time: ~1-2 minutes

## Configuration Parameters

The `ns_merger` section in the JSON config controls:

```json
"ns_merger": {
    "radius": 1.2,           // NS radius (code units)
    "central_density": 2.8,  // Central density ρ_c
    "separation": 6.0,       // Initial separation between star centers
    "collision_velocity": 0.15, // Each star's velocity toward collision (in c)
    "n_radial": 30,          // Radial particle resolution
    "polytropic_index": 1.0  // n=1 for γ=2 polytrope
}
```

## Output

Results are saved to `simulations/astrophysics/ns_merger/results/`:
- `test_quick/` - Quick test results
- `ns_merger_2d/` - Full simulation results

Each run generates:
- Snapshot CSV files (`snapshot_XXXX.csv`)
- Energy log (`energy.dat`)
- Simulation log (`.log`)

## Visualization

The visualization script supports:

```bash
# Plot final snapshot
python3 simulations/astrophysics/ns_merger/scripts/ns_merger_visualize.py \
    --results-dir simulations/astrophysics/ns_merger/results/test_quick --plot

# Create animation
python3 simulations/astrophysics/ns_merger/scripts/ns_merger_visualize.py \
    --results-dir simulations/astrophysics/ns_merger/results/test_quick --animate

# Evolution grid (multiple timesteps)
python3 simulations/astrophysics/ns_merger/scripts/ns_merger_visualize.py \
    --results-dir simulations/astrophysics/ns_merger/results/test_quick --grid

# Density profile evolution
python3 simulations/astrophysics/ns_merger/scripts/ns_merger_visualize.py \
    --results-dir simulations/astrophysics/ns_merger/results/test_quick --evolution
```

## Physical Background

### Lane-Emden Profile

For a polytropic equation of state P = K·ρ^(1+1/n), the hydrostatic equilibrium structure is given by the Lane-Emden equation. For n=1 (γ=2), the density profile is:

```
ρ(r) = ρ_c · sin(πr/R) / (πr/R)
```

where ρ_c is the central density and R is the stellar radius.

### Particle Placement

Particles are placed in concentric rings in 2D:
1. Each ring at radius r has circumference 2πr
2. Number of particles per ring scales with r
3. Total mass per ring is constant (for uniform mass particles)
4. This maintains good neighbor distribution

### SR-GSPH Requirements

Special Relativistic Godunov SPH requires:
- Gaussian kernel (mandatory)
- Baryon number `nu` for each particle (set equal to mass)
- Smoothing length `sml` initialization
- Relativistic energy/velocity handling

## Reference

Shibata, M. & Hotokezaka, K. (2019). "Merger and Mass Ejection of Neutron Star Binaries". Annual Review of Nuclear and Particle Science, 69, 41-64.

## Directory Structure

```
simulations/astrophysics/ns_merger/
├── Makefile.ns_merger
├── README.md
├── config/
│   └── presets/
│       ├── ns_merger_2d.json      # Full simulation config
│       └── ns_merger_quick_test.json  # Quick test config
├── results/
│   ├── test_quick/                # Quick test output
│   │   ├── snapshot_*.csv
│   │   ├── energy.dat
│   │   └── ns_merger.gif
│   └── ns_merger_2d/              # Full simulation output
└── scripts/
    └── ns_merger_visualize.py     # Visualization tool
```
