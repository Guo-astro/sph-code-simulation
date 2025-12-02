# IMBH-Cloud Configuration Presets

This directory contains organized JSON configuration files for the IMBH-Cloud tidal disruption simulations.

## Folder Structure

```
config/presets/
├── README.md              # This file
├── relaxation/            # Initial condition generation
│   ├── ic_relax_10k.json     # Testing: 10k particles
│   └── ic_relax_200k.json    # Production: 200k particles
├── hydrostatic/           # Self-gravity validation
│   ├── ic_verify_10k_hydrostatic.json
│   └── ic_verify_200k_hydrostatic.json
├── simulation/            # Main IMBH tidal disruption runs
│   ├── simulation_10k_b3pc_cool.json     # 10k, b=3pc, WITH cooling
│   ├── simulation_10k_b3pc_nocool.json   # 10k, b=3pc, NO cooling
│   ├── simulation_10k_b1p5pc_optimal.json
│   ├── simulation_200k_b3pc_cool.json    # 200k, b=3pc, WITH cooling
│   └── simulation_200k_b6pc.json
└── continuation/          # Resume from snapshots
    └── continue_10k_b3pc_nocool.json     # Continue from t=2.0 to t=4.0
```

## Naming Conventions

### Relaxation Configs (`relaxation/`)
- Pattern: `ic_relax_<particles>.json`
- Example: `ic_relax_10k.json` → ~10,000 particles

### Hydrostatic Configs (`hydrostatic/`)
- Pattern: `ic_verify_<particles>_hydrostatic.json`
- Example: `ic_verify_10k_hydrostatic.json`

### Simulation Configs (`simulation/`)
- Pattern: `simulation_<particles>_<impact>_<thermal>.json`
- Components:
  - `<particles>`: `10k`, `200k`, etc.
  - `<impact>`: `b3pc` (3 parsec), `b6pc` (6 parsec), `b1p5pc` (1.5 parsec)
  - `<thermal>`: `cool` (with K&I cooling) or `nocool` (adiabatic)
- Example: `simulation_10k_b3pc_nocool.json`

### Continuation Configs (`continuation/`)
- Pattern: `continue_<base_config>.json`
- Example: `continue_10k_b3pc_nocool.json` → extends the nocool run

## How to Add a New Simulation Case

1. **Create the JSON config** in `simulation/`:
   ```bash
   cp simulation/simulation_10k_b3pc_cool.json simulation/simulation_<YOUR_NAME>.json
   # Edit the new file with your parameters
   ```

2. **Register in Makefile** (one line change):
   ```makefile
   # In Makefile.imbh_cloud, find SIMULATION_CASES and add your case:
   SIMULATION_CASES := 10k_b3pc_cool 10k_b3pc_nocool ... <YOUR_NAME>
   ```

3. **Run**:
   ```bash
   make -f sample/imbh_cloud/Makefile.imbh_cloud sim_<YOUR_NAME>
   make -f sample/imbh_cloud/Makefile.imbh_cloud sim_<YOUR_NAME>_viz
   ```

## Key Configuration Parameters

| Parameter | Description | Typical Values |
|-----------|-------------|----------------|
| `endTime` | Simulation end time (code units) | 2.0 = ~2 Myr |
| `outputTime` | Snapshot interval | 0.02 = ~100 snapshots for t=2 |
| `imbh_parameters.M_BH` | Black hole mass (M_☉) | 100000 = 10⁵ M_☉ |
| `imbh_parameters.impact_parameter_b` | Impact parameter (pc) | 3.0 |
| `thermal.enable_cooling` | K&I cooling | true/false |
| `resetTimeOnResume` | For continuation | false = continue time |

## Continuation Runs

To continue a simulation from a snapshot:

1. Create a continuation config in `continuation/`:
   - Set `resumeFromSnapshot` to the snapshot path
   - Set `resetTimeOnResume: false` (continue time, don't restart at t=0)
   - Set `endTime` to the new end time

2. Or use the generic Makefile target:
   ```bash
   make imbh_continue SNAPSHOT=<path> END_TIME=<time>
   ```
