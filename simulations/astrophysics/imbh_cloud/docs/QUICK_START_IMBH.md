# Quick Start: IMBH Simulation from Relaxed Initial Condition

## Overview
This guide shows how to run IMBH-cloud tidal disruption simulations starting from the relaxed Lane-Emden sphere you've already generated.

## Prerequisites

✅ **You have completed:**
1. Relaxation run: `make imbh_relax_2k` → Generated `snapshot_0032.csv`
2. Built SPH code: `cd build && make` → Executable at `build/sph`

## Quick Start (5 minutes)

### Option 1: Using Pre-configured Preset

```bash
# Run IMBH simulation with relaxed 2k particle IC
./build/sph simulations/astrophysics/imbh_cloud/config/presets/imbh_from_relaxed_2k.json

# Monitor progress
tail -f simulations/astrophysics/imbh_cloud/results/imbh_relaxed_2k_b3pc/run.log
```

**What this does:**
- ✓ Loads relaxed cloud from `snapshot_0032.csv`
- ✓ Initializes IMBH at x=-20pc with v=10km/s toward cloud
- ✓ Impact parameter b=3pc (strong disruption)
- ✓ Simulates 2.0 Myr (≈4 cloud crossing times)
- ✓ Outputs 100 snapshots to `results/imbh_relaxed_2k_b3pc/`

### Option 2: Parameter Variations

**Different impact parameters:**

```bash
# Copy template and modify
cp simulations/astrophysics/imbh_cloud/config/presets/imbh_from_relaxed_2k.json \
   simulations/astrophysics/imbh_cloud/config/my_imbh_b4pc.json

# Edit my_imbh_b4pc.json:
#   - Change "impact_parameter_b": 4.0
#   - Change output_directory to "results/imbh_relaxed_2k_b4pc"

# Run
./build/sph simulations/astrophysics/imbh_cloud/config/my_imbh_b4pc.json
```

**Different velocities:**

```bash
# For v = 5 km/s:
#   - Change "BH_initial_velocity": [5.0, 0.0, 0.0]
#   - Adjust end_time to 3.0 or 4.0 Myr (slower encounter)

# For v = 20 km/s:
#   - Change "BH_initial_velocity": [20.0, 0.0, 0.0]
#   - Can reduce end_time to 1.0 Myr (faster encounter)
```

## Expected Results

### Timeline

```
t = 0.0 Myr:  Cloud at origin (relaxed equilibrium)
              BH at x = -20 pc, velocity = +10 km/s

t ≈ 1.5 Myr:  BH reaches pericenter (closest approach at y = b)
              Tidal forces strongest → cloud compression

t ≈ 2.0 Myr:  Post-encounter phase
              Tidal tails formed, possible accretion onto BH
```

### Output Files

```
simulations/astrophysics/imbh_cloud/results/imbh_relaxed_2k_b3pc/
├── snapshot_0000.csv  # t = 0.00 Myr (initial state from relaxed IC)
├── snapshot_0001.csv  # t = 0.02 Myr
├── ...
├── snapshot_0100.csv  # t = 2.00 Myr (final state)
├── energy.dat         # Energy conservation tracking
└── run.log           # Simulation log
```

## Verification Checks

### 1. Energy Conservation

```bash
# Plot energy vs time
python3 << 'EOF'
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('simulations/astrophysics/imbh_cloud/results/imbh_relaxed_2k_b3pc/energy.dat', 
                  comments='#')
time = data[:, 0]
E_kin = data[:, 1]
E_thermal = data[:, 2]
E_pot = data[:, 3]
E_total = E_kin + E_thermal + E_pot

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot(time, E_kin, label='Kinetic')
plt.plot(time, E_thermal, label='Thermal')
plt.plot(time, E_pot, label='Potential')
plt.plot(time, E_total, 'k--', label='Total', linewidth=2)
plt.xlabel('Time [Myr]')
plt.ylabel('Energy [code units]')
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(1, 2, 2)
dE = (E_total - E_total[0]) / abs(E_total[0]) * 100
plt.plot(time, dE)
plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
plt.axhline(y=1, color='r', linestyle='--', alpha=0.3, label='1% threshold')
plt.axhline(y=-1, color='r', linestyle='--', alpha=0.3)
plt.xlabel('Time [Myr]')
plt.ylabel('ΔE/E₀ [%]')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('energy_conservation.png', dpi=150)
print("✓ Energy plot saved to energy_conservation.png")
EOF
```

**Good conservation:** ΔE/E₀ < 1% over entire simulation

### 2. Cloud Stability Before Encounter

```bash
# Check first few snapshots - cloud should remain stable
python3 << 'EOF'
import pandas as pd
import numpy as np

# Load initial and t=0.1 Myr snapshots
snap0 = pd.read_csv('simulations/astrophysics/imbh_cloud/results/imbh_relaxed_2k_b3pc/snapshot_0000.csv',
                    comment='#')
snap5 = pd.read_csv('simulations/astrophysics/imbh_cloud/results/imbh_relaxed_2k_b3pc/snapshot_0005.csv',
                    comment='#')

# Check if cloud center-of-mass remains near origin
r0 = np.sqrt(snap0['x']**2 + snap0['y']**2 + snap0['z']**2)
r5 = np.sqrt(snap5['x']**2 + snap5['y']**2 + snap5['z']**2)

print(f"Initial cloud radius: <r> = {r0.mean():.3f} pc, σ = {r0.std():.3f} pc")
print(f"At t=0.1 Myr:         <r> = {r5.mean():.3f} pc, σ = {r5.std():.3f} pc")
print(f"Change: Δr = {abs(r5.mean() - r0.mean()):.4f} pc")

if abs(r5.mean() - r0.mean()) < 0.1:
    print("✓ Cloud remains stable before BH encounter")
else:
    print("⚠ Warning: Cloud may be expanding/contracting")
EOF
```

### 3. Tidal Disruption Diagnostics

```bash
# Check final snapshot for tidal features
python3 << 'EOF'
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load final snapshot
snap = pd.read_csv('simulations/astrophysics/imbh_cloud/results/imbh_relaxed_2k_b3pc/snapshot_0100.csv',
                   comment='#')

# Spatial distribution
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# XY projection
axes[0].scatter(snap['x'], snap['y'], c=snap['rho'], s=1, cmap='viridis', alpha=0.6)
axes[0].set_xlabel('x [pc]')
axes[0].set_ylabel('y [pc]')
axes[0].set_title('XY Projection (density-colored)')
axes[0].set_aspect('equal')
axes[0].grid(True, alpha=0.3)

# XZ projection
axes[1].scatter(snap['x'], snap['z'], c=snap['rho'], s=1, cmap='viridis', alpha=0.6)
axes[1].set_xlabel('x [pc]')
axes[1].set_ylabel('z [pc]')
axes[1].set_title('XZ Projection')
axes[1].set_aspect('equal')
axes[1].grid(True, alpha=0.3)

# Density histogram
axes[2].hist(np.log10(snap['rho']), bins=50, alpha=0.7, edgecolor='black')
axes[2].set_xlabel('log₁₀(ρ) [code units]')
axes[2].set_ylabel('Particle count')
axes[2].set_title('Density Distribution')
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('tidal_disruption_final.png', dpi=150)
print("✓ Final state plot saved to tidal_disruption_final.png")
print(f"  Particles: {len(snap)}")
print(f"  Density range: {snap['rho'].min():.2e} to {snap['rho'].max():.2e}")
print(f"  Spatial extent: x ∈ [{snap['x'].min():.1f}, {snap['x'].max():.1f}] pc")
EOF
```

## Troubleshooting

### "Error: IMBH external force module not found"

**Solution:** The IMBH point-mass external force needs to be implemented. Check:

```bash
# Search for existing implementation
grep -r "PointMassBH\|IMBHForce" include/ src/

# If not found, you'll need to implement it
# See RESEARCH_SETUP.md section "Physics Modules Required"
```

### "Checkpoint file not found"

**Solution:** Verify the relaxed snapshot exists:

```bash
ls -lh simulations/astrophysics/imbh_cloud/results/lane_emden_2k_relax/snapshot_0032.csv

# If missing, run relaxation first:
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_relax_2k
```

### "Simulation crashes immediately"

**Check:**
1. ✓ Is `use_gravity: true` enabled?
2. ✓ Does relaxed snapshot have proper format (CSV with header)?
3. ✓ Is gamma = 1.6666... (5/3) consistent with relaxation?

```bash
# Verify snapshot header
head -30 simulations/astrophysics/imbh_cloud/results/lane_emden_2k_relax/snapshot_0032.csv
```

## Next Steps

1. **Run baseline case** (b=3pc, v=10km/s)
2. **Parameter scan:** Vary b = 3, 4, 5, 6 pc
3. **Velocity scan:** Vary v = 0, 5, 10, 15, 20 km/s
4. **Higher resolution:** Generate and use 200k particle relaxed IC
5. **Implement visualization:** Create animations of tidal disruption

For detailed physics and resolution requirements, see:
- `simulations/astrophysics/imbh_cloud/docs/RESEARCH_SETUP.md` - Full research setup guide
- `simulations/astrophysics/imbh_cloud/README.md` - Project overview
- `simulations/astrophysics/imbh_cloud/QUICK_START.md` - Quick reference

## Contact & Support

- **Issue tracker:** Report problems with simulation setup
- **Documentation:** See `docs/` directory for physics details
- **Examples:** Check `sample/` for other test cases
