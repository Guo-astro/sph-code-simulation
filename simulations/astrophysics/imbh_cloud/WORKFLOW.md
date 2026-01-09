# IMBH-Cloud Simulation Workflow

## Quick Start

```bash
make relax       # Step 1: Position particles (no gravity)
make hydro       # Step 2: Test stability with self-gravity
make interact    # Step 3: Add IMBH and run flyby
```

Or run all:
```bash
make run
```

## Three-Phase Workflow

### Phase 1: Relaxation (no self-gravity)

**Goal**: Position particles to match density profile

```bash
make relax
```

**What happens**:
- Create cloud with ghost envelope for pressure boundary
- Run GLASS relaxation to fix particle spacing
- Run main relaxation: SPH forces - analytical equilibrium forces
- Velocities zeroed each step (quasi-static positioning)
- Output: particles positioned correctly, ghost particles removed

**Self-gravity**: OFF (only SPH pressure forces)

**Output**: `results/be_sphere_relaxation/`

### Phase 2: Hydrostatic Test (self-gravity ON, no IMBH)

**Goal**: Verify cloud is stable with self-gravity

```bash
make hydro
```

**What happens**:
- Load relaxed cloud from Phase 1
- Enable self-gravity
- Run for ~1 Myr without IMBH
- Check if cloud remains stable

**Success criteria**:
- Cloud radius stays approximately constant
- Velocity dispersion remains small (<0.5 km/s)
- No collapse or expansion

**If unstable**: Need more relaxation steps or check pressure matching

**Output**: `results/be_sphere_hydrostatic/`

### Phase 3: IMBH Flyby (self-gravity ON, IMBH ON)

**Goal**: Simulate tidal interaction

```bash
make interact
```

**What happens**:
- Load stable cloud from Phase 2 (or Phase 1 if Phase 2 passed)
- Add IMBH at 10 pc
- Run tidal interaction simulation

**Output**: `results/be_sphere_imbh/`

## Why Three Phases?

```
Phase 1 (Relaxation):
  - Only positions particles using analytical equilibrium
  - Self-gravity is OFF
  - Does NOT verify true hydrostatic balance

Phase 2 (Hydrostatic Test):
  - Turns ON self-gravity
  - Checks: pressure gradient = gravitational force?
  - If cloud expands/collapses → relaxation failed

Phase 3 (IMBH Flyby):
  - Only run after Phase 2 confirms stability
  - IMBH tidal forces now act on stable cloud
```

## Commands

| Command | Description |
|---------|-------------|
| `make relax` | Phase 1: Relaxation |
| `make hydro` | Phase 2: Hydrostatic test |
| `make interact` | Phase 3: IMBH flyby |
| `make run` | Run all three phases |
| `make clean` | Delete results |

## Physical Parameters

Cloud:
- Mass: 40 Msun
- Temperature: 15 K
- Central density: ~200 cm^-3
- Profile: rho = rho_c / (1 + (r/r_c)^2)

IMBH:
- Mass: 10^5 Msun
- Initial distance: 10 pc
