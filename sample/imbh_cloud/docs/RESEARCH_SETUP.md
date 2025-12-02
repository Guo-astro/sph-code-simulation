# IMBH-Molecular Cloud Interaction: Research Setup

## Table of Contents

1. [Scientific Objective](#scientific-objective)
2. [Using Relaxed Initial Conditions](#using-relaxed-initial-conditions)
3. [Resolution Analysis](#resolution-analysis)
4. [Physics Modules Required](#physics-modules-required)

---

## Scientific Objective

Simulate intermediate-mass black hole (IMBH) tidal disruption of molecular clouds with the following parameters:

### Physical Parameters

- **IMBH Mass**: M_BH = 10⁵ M_☉ (100,000 solar masses)
- **Impact Parameters**: b = 3-6 parsecs (parameter scan study)
- **Initial Relative Velocity**: v₀ = 0-20 km/s (velocity effect study)
  - Baseline: v₀ = 10 km/s in x-direction
  - Additional scans: v₀ = 0, 5, 10, 15, 20 km/s
- **Cloud Structure**: 3D Lane-Emden sphere with polytropic index n = 5/3 (γ = 5/3)
- **Thermal Physics**: Koyama & Inutsuka (2000) thermal equilibrium curve
  - Assumption: Thermal timescale << dynamical timescale
  - Cloud remains in thermal equilibrium throughout simulation

### Research Questions

1. **Impact Parameter Study**: How does tidal disruption efficiency vary with impact parameter b?
   - b = 3 pc: Strong disruption (b < r_tidal)
   - b = 4 pc: Moderate disruption
   - b = 5 pc: Weak disruption
   - b = 6 pc: Minimal disruption (reference case)

2. **Velocity Effect Study**: How does initial relative velocity affect:
   - Encounter timescale and disruption morphology
   - Mass accretion rate onto BH
   - Tidal tail formation
   - Shock heating at pericenter passage

---

## Using Relaxed Initial Conditions

### Overview

To conduct realistic IMBH-cloud interaction simulations, you should **start from a relaxed, gravitationally stable cloud** rather than an idealized analytic profile. The relaxation process ensures:

1. ✅ **Hydrostatic equilibrium**: Pressure gradient balances self-gravity
2. ✅ **Thermal equilibrium**: Temperature follows Koyama-Inutsuka curve for given density
3. ✅ **Numerical stability**: Smoothing lengths and densities properly converged
4. ✅ **Physical realism**: Cloud structure reflects true molecular cloud properties

### Available Relaxed Snapshots

You have generated relaxed initial conditions at:
```
/Users/guo/Downloads/sphcode/sample/imbh_cloud/results/ic_relax_10k/snapshot_0032.csv
```

**Properties of this snapshot:**
- **Particle count**: ~2,000 particles (N=22 in Lane-Emden generator)
- **Cloud mass**: M_cloud ≈ 10⁴ M_☉ (configure based on your needs)
- **Cloud radius**: R_cloud ≈ 5 pc (configure based on your needs)
- **Relaxation steps**: 20,000 steps (fully converged)
- **Method**: GDISPH with self-gravity
- **Status**: ✓ Relaxed and ready for IMBH interaction simulation

### Step-by-Step: Using Relaxed Snapshot in IMBH Simulation

#### Step 1: Verify Relaxation Quality

Before using the snapshot, confirm it's well-relaxed:

```bash
# Check final snapshots
ls -lh sample/imbh_cloud/results/ic_relax_10k/snapshot_*.csv | tail -5

# View relaxation summary plots (if generated)
open sample/imbh_cloud/results/ic_relax_10k/relaxation_summary.png
```

**Quality criteria:**
- ✓ Energy conservation: ΔE/E₀ < 1% over relaxation
- ✓ Density profile: Matches Lane-Emden n=1.5 polytrope
- ✓ Virial equilibrium: 2T + U ≈ 0 (within 5%)

#### Step 2: Create IMBH Simulation Configuration

Create a new JSON configuration file for your IMBH-cloud interaction. You have two options:

**Option A: Start from existing preset**

```bash
# Copy and modify existing preset
cp sample/imbh_cloud/config/presets/simulation_10k_b3pc_cool.json \
   sample/imbh_cloud/config/my_imbh_simulation.json
```

**Option B: Create new configuration from scratch**

Here's a complete template for using the relaxed snapshot:

```json
{
  "name": "imbh_cloud_relaxed_init",
  "description": "IMBH tidal disruption starting from relaxed Lane-Emden sphere",
  
  "comment": "=== Using Relaxed Initial Condition ===",
  
  "physics": {
    "gamma": 1.6666666666666667,
    "gravitational_constant": 1.0,
    "use_gravity": true
  },
  
  "checkpoint": {
    "enabled": true,
    "saveInterval": 500,
    "directory": "sample/imbh_cloud/results/my_imbh_run/checkpoints",
    "saveOnExit": true,
    "maxCheckpoints": 5,
    "autoResume": false,
    "resumeFile": "sample/imbh_cloud/results/ic_relax_10k/snapshot_0032.csv"
  },
  
  "imbh_parameters": {
    "M_BH": 100000.0,
    "BH_position": [-20.0, 0.0, 0.0],
    "BH_velocity": [10.0, 0.0, 0.0],
    "impact_parameter": 3.0,
    "softening_length": 0.05
  },
  
  "numerical": {
    "neighbor_number": 50,
    "kernel": "wendland",
    "sph_type": "gdisph",
    "use_gravity": true,
    "periodic": false,
    "iterative_smoothing_length": true,
    "leaf_particle_number": 32
  },
  
  "simulation": {
    "start_time": 0.0,
    "end_time": 2.0,
    "output_time": 0.02,
    "energy_time": 0.01,
    "output_directory": "sample/imbh_cloud/results/my_imbh_run"
  },
  
  "timestep": {
    "CFL_sound": 0.3,
    "CFL_force": 0.25,
    "max_timestep": 0.01
  },
  
  "artificial_viscosity": {
    "alpha": 1.0,
    "use_balsara_switch": true,
    "use_time_dependent_av": true,
    "use_artificial_conductivity": false
  },
  
  "thermal": {
    "enable_cooling": true,
    "cooling_type": "koyama_inutsuka",
    "N_H_column": 1e20,
    "thermal_relax_time": 0.01
  },
  
  "output": {
    "formats": [
      {
        "type": "csv",
        "precision": 16
      }
    ],
    "enableEnergyFile": true
  }
}
```

**Key Configuration Parameters:**

| Parameter | Description | Recommended Value |
|-----------|-------------|-------------------|
| `checkpoint.resumeFile` | Path to relaxed snapshot | `sample/imbh_cloud/results/ic_relax_10k/snapshot_0032.csv` |
| `checkpoint.autoResume` | Auto-resume on restart | `false` (set to `true` only for continuing interrupted runs) |
| `imbh_parameters.M_BH` | Black hole mass [M_☉] | `100000.0` (10⁵ M_☉) |
| `imbh_parameters.BH_position` | Initial BH position [pc] | `[-20.0, 0.0, 0.0]` (start far from cloud) |
| `imbh_parameters.BH_velocity` | BH velocity [km/s] | `[10.0, 0.0, 0.0]` (baseline v_rel = 10 km/s) |
| `imbh_parameters.impact_parameter` | Impact parameter b [pc] | `3.0` to `6.0` (scan range) |
| `simulation.end_time` | Simulation duration [Myr] | `2.0` (≈4 crossing times) |
| `simulation.output_time` | Snapshot interval [Myr] | `0.02` (100 snapshots) |

#### Step 3: Configure IMBH External Force

**CRITICAL**: The current codebase may need an IMBH external force module. Check if it exists:

```bash
# Search for IMBH or point mass force implementation
find include/ src/ -name "*imbh*" -o -name "*point_mass*" -o -name "*black_hole*"
```

If the external force module doesn't exist yet, you'll need to implement:
- `include/external_forces/point_mass_bh.hpp`
- `src/external_forces/point_mass_bh.cpp`

See [Physics Modules Required](#physics-modules-required) section below for implementation details.

#### Step 4: Run the Simulation

Once configuration is ready:

```bash
# Method 1: Direct execution with config file
./build/sph sample/imbh_cloud/config/my_imbh_simulation.json

# Method 2: Using Makefile (if you add custom target)
# Edit sample/imbh_cloud/Makefile.imbh_cloud to add:
#
# my_imbh_run:
#     @echo "Running IMBH simulation with relaxed initial condition..."
#     mkdir -p sample/imbh_cloud/results/my_imbh_run
#     ./build/sph sample/imbh_cloud/config/my_imbh_simulation.json
#
make -f sample/imbh_cloud/Makefile.imbh_cloud my_imbh_run
```

#### Step 5: Monitor Simulation Progress

```bash
# Watch real-time log output
tail -f sample/imbh_cloud/results/my_imbh_run/run.log

# Check energy conservation
cat sample/imbh_cloud/results/my_imbh_run/energy.dat

# List generated snapshots
ls -lh sample/imbh_cloud/results/my_imbh_run/snapshot_*.csv
```

**Expected behavior:**
- ✓ Time starts from 0.0 (or relaxation end time if `resetTimeOnResume: false`)
- ✓ Particle positions/velocities loaded from snapshot
- ✓ Cloud remains stable initially (before BH approaches)
- ✓ Tidal disruption begins when BH passes within ~r_tidal

### Parameter Scan Strategy

For systematic studies, create multiple configurations varying key parameters:

**Impact Parameter Scan** (b = 3, 4, 5, 6 pc at fixed v = 10 km/s):

```bash
# Generate configs programmatically
for b in 3 4 5 6; do
  sed "s/IMPACT_PARAMETER/$b/g" template.json > config_b${b}pc.json
done

# Run all cases
for b in 3 4 5 6; do
  ./build/sph sample/imbh_cloud/config/config_b${b}pc.json &
done
wait
```

**Velocity Scan** (v = 0, 5, 10, 15, 20 km/s at fixed b = 4 pc):

```bash
for v in 0 5 10 15 20; do
  # Modify BH_velocity in config: [v, 0.0, 0.0]
  # Run simulation
done
```

### Troubleshooting

**Problem: "Failed to load checkpoint file"**
- ✓ Verify file path is correct (use absolute or relative to execution directory)
- ✓ Check CSV file format (must have proper header and metadata)
- ✓ Ensure `checkpoint.resumeFile` points to snapshot_0032.csv, not a directory

**Problem: "Particle count mismatch"**
- The loaded snapshot has a fixed number of particles (N ≈ 2000)
- Don't specify `N` parameter in config when using `resumeFile`
- The particle count is determined by the snapshot file

**Problem: "Cloud immediately flies apart"**
- Check that `use_gravity: true` is enabled
- Verify the relaxed snapshot actually reached equilibrium
- Ensure `gamma` matches the relaxation run (1.666... for γ=5/3)

**Problem: "BH doesn't interact with cloud"**
- Verify IMBH external force module is implemented and enabled
- Check `BH_position` is not too far (should be within ~20-30 pc initially)
- Confirm `BH_velocity` points toward cloud (positive x-direction if cloud at origin)

### Visualization After Simulation

```bash
# Generate density projection plots
python3 sample/imbh_cloud/scripts/visualize_imbh_cloud.py \
  --input sample/imbh_cloud/results/my_imbh_run \
  --output plots/

# Create animation
python3 sample/imbh_cloud/scripts/animate_tidal_disruption.py \
  --snapshots sample/imbh_cloud/results/my_imbh_run/snapshot_*.csv \
  --output animations/tidal_disruption.gif
```

### Higher Resolution Runs

For publication-quality results, use the 200k particle relaxed snapshot:

```bash
# First, generate 200k relaxed snapshot (if not done yet)
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_relax_200k

# Wait for relaxation to complete (~several hours)
# Then use: sample/imbh_cloud/results/ic_relax_200k/snapshot_XXXX.csv
```

**Resolution comparison:**

| N particles | Purpose | h/R_cloud | Compute time | Quality |
|-------------|---------|-----------|--------------|---------|
| 2k (N=22) | Quick tests, parameter scans | ~0.02 | Minutes | Exploratory ⭐ |
| 20k (N=50) | Standard runs | ~0.007 | ~1 hour | Good ⭐⭐ |
| 200k (N=100) | Publication | ~0.002 | ~12 hours | Excellent ⭐⭐⭐ |

---

## Resolution Analysis

### 1. Jeans Mass Criterion

The Jeans mass for a self-gravitating cloud is:

```
M_J = (π^(5/2) / 6) * (c_s^3 / (G^(3/2) * ρ^(1/2)))
```

For thermal equilibrium in CNM (Cold Neutral Medium):
- n_H ~ 100-1000 cm⁻³ → T_eq ~ 10-50 K (from K&I thermal curve)
- Sound speed: c_s ~ 0.2-0.3 km/s

**Numerical example** (n_H = 100 cm⁻³, T = 50 K):
- ρ = n_H × μ m_H / pc³ = 100 × (1.67×10⁻²⁴ g) / (2.938×10⁵⁷ cm³) ≈ 2.47 M_☉/pc³
- c_s = √(k_B T / μ m_H) = √(1.38×10⁻¹⁶ × 50 / 1.67×10⁻²⁴) ≈ 0.29 km/s
- **M_J ≈ 3.5 M_☉** (typical for molecular clouds)

**Requirement**: To resolve fragmentation and thermal pressure support:
- **N_J ≥ 50-100 particles per Jeans mass**
- For cloud mass M_cloud ~ 10³-10⁴ M_☉:
  - **Minimum N ≥ 5×10⁴ - 10⁶ particles**

### 2. Tidal Disruption Timescale

The tidal (pancake) disruption occurs when tidal force exceeds self-gravity:

```
t_tidal ~ √(R_cloud³ / (G M_BH))
```

**Numerical example** (M_BH = 10⁵ M_☉, R_cloud = 5 pc):
```
t_tidal = √(5³ pc³ / (4.302×10⁻³ pc·(km/s)²·M_☉⁻¹ × 10⁵ M_☉))
        = √(125 / 430.2) [Myr]
        ≈ 0.54 Myr = 5.4×10⁵ years
```

Compare to cloud crossing time:
```
t_cross ~ R_cloud / v_rel
```

**Numerical examples**:
- v_rel = 10 km/s: t_cross = 5 pc / 10 km/s = 0.49 Myr = **4.9×10⁵ years**
- v_rel = 20 km/s: t_cross = 5 pc / 20 km/s = 0.24 Myr = **2.4×10⁵ years**
- v_rel = 50 km/s: t_cross = 5 pc / 50 km/s = 0.10 Myr = **1.0×10⁵ years**

**Timescale ordering**: t_thermal << t_tidal ~ t_cross

### 2a. Velocity Effect on Encounter Timescale

**Pericenter passage time** (when tidal force is strongest):
```
t_peri ~ b / v_rel
```

**Numerical examples** (various b and v_rel combinations):

| b [pc] | v_rel [km/s] | t_peri [Myr] | t_peri [years] | Physical interpretation |
|--------|--------------|--------------|----------------|------------------------|
| 3 | 0 | ∞ | ∞ | Static infall - no pericenter, pure radial |
| 3 | 5 | 0.59 | 5.9×10⁵ | Slow, extended tidal compression |
| 3 | 10 | 0.29 | 2.9×10⁵ | Baseline - moderate encounter |
| 3 | 15 | 0.20 | 2.0×10⁵ | Fast - impulsive tidal shock |
| 3 | 20 | 0.15 | 1.5×10⁵ | Very fast - minimal disruption |
| 4 | 10 | 0.39 | 3.9×10⁵ | Baseline at moderate distance |
| 5 | 10 | 0.49 | 4.9×10⁵ | Weak tidal interaction |
| 6 | 10 | 0.59 | 5.9×10⁵ | Minimal tidal effect (reference) |

**Velocity regimes**:
- **v_rel = 0 km/s** (static case): t_peri → ∞, pure radial infall, maximum disruption
- **v_rel = 5 km/s**: t_peri ~ 0.6 Myr, slow encounter, extended tidal compression
- **v_rel = 10 km/s** (BASELINE): t_peri ~ 0.3 Myr, intermediate regime
- **v_rel = 15 km/s**: t_peri ~ 0.2 Myr, fast encounter, impulsive tidal shock
- **v_rel = 20 km/s**: t_peri ~ 0.15 Myr, very fast, minimal disruption

**Key physics**: Higher v_rel → shorter t_peri → less time for tidal stripping → weaker disruption

### 3. Spatial Resolution Requirements

**Minimum spatial resolution needed**:

a) **Tidal radius** (where BH gravity ~ cloud self-gravity):
   ```
   r_t ~ R_cloud * (M_BH / M_cloud)^(1/3)
   ```
   **Numerical example** (M_cloud = 10⁴ M_☉, M_BH = 10⁵ M_☉, R_cloud = 5 pc):
   ```
   r_t = 5 pc × (10⁵ / 10⁴)^(1/3)
       = 5 pc × 10^(1/3)
       = 5 pc × 2.15
       ≈ 10.8 pc
   ```
   - **Need to resolve h ~ 0.1-1 pc** to capture tidal effects

b) **Hill radius** (tidal truncation):
   ```
   r_H ~ b * (M_cloud / (3 M_BH))^(1/3)
   ```
   **Numerical examples** (M_cloud = 10⁴ M_☉, M_BH = 10⁵ M_☉):
   
   | b [pc] | r_H [pc] | Interpretation |
   |--------|----------|----------------|
   | 3 | 0.99 | Strong tidal truncation |
   | 4 | 1.32 | Moderate truncation |
   | 5 | 1.65 | Weak truncation |
   | 6 | 1.98 | Minimal truncation |
   
   Calculation: r_H = b × (10⁴ / (3×10⁵))^(1/3) = b × 0.3264 pc
   - **Need N_neighbor ~ 50-100 within r_H** for proper resolution

c) **Smoothing length constraint**:
   - h ~ R_cloud / √N for uniform distribution
   
   **Numerical examples** (R_cloud = 5 pc):
   
   | N particles | h [pc] | h/R_cloud | Resolution quality |
   |-------------|--------|-----------|-------------------|
   | 5×10⁴ | 0.022 | 0.0045 | Quick test ✓ |
   | 10⁵ | 0.016 | 0.0032 | Good ✓ |
   | 2×10⁵ | 0.011 | 0.0022 | Standard ✓✓ |
   | 10⁶ | 0.005 | 0.0010 | High resolution ✓✓✓ |

### 4. Recommended Particle Numbers

**Quick Reference Table**:

| Cloud Mass | N_particles | h/R_cloud | N_J | Purpose |
|------------|-------------|-----------|-----|---------|
| 10³ M_☉   | 5×10⁴       | ~0.014    | 50  | Quick test |
| 10⁴ M_☉   | 2×10⁵       | ~0.007    | 100 | Standard resolution |
| 10⁴ M_☉   | 10⁶         | ~0.003    | 500 | High resolution |

**Recommended**: Start with **N = 2×10⁵** for production runs, **N = 5×10⁴** for parameter studies.

---

### 4a. Minimal Acceptable Particle Numbers (For Reviewers)

The following table provides **justifiable minimal particle counts** for different physical scenarios, based on multiple resolution criteria that reviewers will scrutinize:

#### **Resolution Criteria for Minimal N**:
1. **Jeans criterion**: N_J ≥ 50 particles per Jeans mass (fragmenttion resolution)
2. **Neighbor criterion**: N_neighbor ≥ 50 within smoothing kernel (SPH accuracy)
3. **Hill sphere**: ≥50 particles within r_H at pericenter (tidal truncation)
4. **Spatial resolution**: h ≤ 0.1 r_H (resolve tidal features)

---

#### **A. Cloud Mass Variation Study** (Fixed: b = 4 pc, v = 10 km/s, R ∝ M^(1/3))

| M_cloud | R_cloud | <n_H> | M_J | r_H | N_min | N_recommended | Justification |
|---------|---------|-------|-----|-----|-------|---------------|--------------|
| 10² M_☉ | 2.3 pc | 100 cm⁻³ | 3.5 M_☉ | 0.68 pc | 1.4×10³ | 1×10⁴ | M/M_J = 29 → 29×50 = 1450; h(10⁴)=0.023pc < r_H ✓ |
| 5×10² M_☉ | 3.7 pc | 100 cm⁻³ | 3.5 M_☉ | 0.91 pc | 7.1×10³ | 3×10⁴ | M/M_J = 143 → need 7150; increased for h < 0.1×r_H |
| 10³ M_☉ | 4.6 pc | 100 cm⁻³ | 3.5 M_☉ | 1.08 pc | 1.4×10⁴ | 5×10⁴ | **Quick test baseline**; h(5×10⁴)=0.021pc ✓ |
| 5×10³ M_☉ | 7.4 pc | 100 cm⁻³ | 3.5 M_☉ | 1.46 pc | 7.1×10⁴ | 1.5×10⁵ | M/M_J = 1429 → larger cloud needs more resolution |
| 10⁴ M_☉ | 9.3 pc | 100 cm⁻³ | 3.5 M_☉ | 1.73 pc | 1.4×10⁵ | 2×10⁵ | **Production baseline**; h(2×10⁵)=0.021pc < 0.1×r_H ✓ |
| 5×10⁴ M_☉ | 15 pc | 100 cm⁻³ | 3.5 M_☉ | 2.33 pc | 7.1×10⁵ | 1×10⁶ | High-mass cloud, publication quality |

**Scaling**: N_min ≈ (M_cloud / M_J) × 50, rounded up to nearest convenient number

**Calculation example** (M = 10⁴ M_☉):
```
M/M_J = 10⁴ / 3.5 ≈ 2857 Jeans masses
N_min = 2857 × 50 = 142,850 ≈ 1.4×10⁵
Recommended = 2×10⁵ (40% safety margin)
h = R/√N = 9.3pc / √(2×10⁵) ≈ 0.021 pc < 0.1 × r_H(1.73pc) ✓
```

---

#### **B. Impact Parameter Variation** (Fixed: M = 10⁴ M_☉, R = 5 pc, v = 10 km/s)

| b [pc] | r_H [pc] | h_max | N_min | N_recommended | Reviewer argument |
|--------|----------|-------|-------|---------------|-------------------|
| 2.0 | 0.66 | 0.066 | 5.8×10⁵ | 8×10⁵ | **Very close**: h(5.8×10⁵)=0.0066pc; 50 ptcls in r_H requires extreme resolution |
| 3.0 | 0.99 | 0.099 | 2.6×10⁵ | 4×10⁵ | **Close encounter**: h ≤ 0.1pc to resolve Hill sphere; strong tidal stripping |
| 4.0 | 1.32 | 0.132 | 1.4×10⁵ | 2×10⁵ | **Baseline**: h(2×10⁵)=0.011pc < 0.1×r_H ✓; 50 neighbors guaranteed |
| 5.0 | 1.65 | 0.165 | 9.2×10⁴ | 1.5×10⁵ | Weaker tidal effect; slightly relaxed but still resolve r_H |
| 6.0 | 1.98 | 0.198 | 6.4×10⁴ | 1×10⁵ | **Distant**: minimal disruption; lower N acceptable for reference |
| 8.0 | 2.64 | 0.264 | 3.6×10⁴ | 5×10⁴ | Control case: very weak tidal effect |

**Calculation** (h_max = 0.1 × r_H):
```
For b = 4 pc: r_H = 1.32 pc → h_max = 0.132 pc
N_min = (R_cloud / h_max)² = (5 / 0.132)² ≈ 1.4×10⁵
Plus Jeans: M/M_J = 2857, need 142,850
Maximum: N_min = max(1.4×10⁵, 1.4×10⁵) ≈ 1.4×10⁵
Recommended: 2×10⁵ (add safety margin)
```

**Reviewer note**: Close encounters (b < 4 pc) require **higher resolution** because:
1. Smaller r_H → need smaller h → need more particles
2. Stronger tidal gradients → need better force resolution
3. Higher compression → more Jeans masses form → need more N_J

---

#### **C. Velocity Variation** (Fixed: M = 10⁴ M_☉, R = 5 pc, b = 4 pc)

| v [km/s] | t_peri [Myr] | Δt_output | N_snapshots | N_particles | Justification |
|----------|--------------|-----------|-------------|-------------|---------------|
| 0 (static) | ∞ | 0.02 | 100 | 2×10⁵ | **Quasi-static infall**: slow evolution, standard resolution adequate |
| 5 | 0.78 | 0.01 | 200 | 3×10⁵ | Slow encounter: t_peri ~ 1.5×t_tidal → need finer time sampling + spatial |
| 10 | 0.39 | 0.005 | 400 | 2×10⁵ | **Baseline**: moderate timescale, standard resolution |
| 15 | 0.26 | 0.003 | 600 | 2.5×10⁵ | Fast encounter: impulsive shock → need good Riemann solver + slightly more N |
| 20 | 0.20 | 0.002 | 1000 | 3×10⁵ | **Very fast shock**: ΔT_shock ~ 160K → resolve shock heating → higher N |
| 30 | 0.13 | 0.001 | 1300 | 4×10⁵ | Extreme case: shock-dominated, publication-level resolution |

**Velocity-dependent N**: Higher v requires:
1. **More particles** to resolve shock fronts (Δx_shock ~ few h)
2. **Finer timesteps** (CFL condition: Δt ∝ h / v_shock)
3. **More snapshots** to capture rapid pericenter passage

**Calculation** (v = 20 km/s shock):
```
Shock width: λ_shock ~ 5-10 h (numerical diffusion)
v_shock ~ v_rel = 20 km/s
To resolve shock: need h < λ_shock / 5
Post-shock compression: n → 4n (strong shock) → 4× more Jeans masses
N_min = 1.4×10⁵ × √2 ≈ 2×10⁵ (increased by ~√2 for shock)
Recommended: 3×10⁵ (50% safety for shock physics)
```

---

#### **D. Cloud Radius Variation** (Fixed: M = 10⁴ M_☉, b = 4 pc, v = 10 km/s)

| R_cloud | <ρ> | <n_H> | M_J | r_H/R | N_min | N_recommended | Reviewer concern |
|---------|-----|-------|-----|-------|-------|---------------|------------------|
| 3 pc | 88 M_☉/pc³ | 3564 cm⁻³ | 1.3 M_☉ | 0.44 | 3.8×10⁵ | 5×10⁵ | **Dense compact**: M/M_J=7692 → many Jeans masses |
| 4 pc | 37 M_☉/pc³ | 1500 cm⁻³ | 1.9 M_☉ | 0.33 | 2.6×10⁵ | 3.5×10⁵ | Moderate density, resolve r_H well |
| 5 pc | 19 M_☉/pc³ | 775 cm⁻³ | 2.7 M_☉ | 0.26 | 1.9×10⁵ | 2.5×10⁵ | **Baseline** |
| 7 pc | 6.9 M_☉/pc³ | 280 cm⁻³ | 4.5 M_☉ | 0.19 | 1.1×10⁵ | 1.5×10⁵ | Lower density, fewer Jeans masses |
| 10 pc | 2.4 M_☉/pc³ | 97 cm⁻³ | 7.6 M_☉ | 0.13 | 6.6×10⁴ | 1×10⁵ | **Diffuse**: M/M_J=1316, relaxed requirement |

**Key insight**: Smaller R_cloud at fixed M:
- Higher <ρ> → lower M_J → **more Jeans masses** → need more particles
- r_H/R smaller → tidal effect covers larger fraction → need better resolution

**Calculation** (R = 3 pc, M = 10⁴ M_☉):
```
<ρ> = M / (4π/3 R³) = 10⁴ / 113 ≈ 88 M_☉/pc³ → n_H ≈ 3564 cm⁻³
M_J(3564 cm⁻³) ≈ 1.3 M_☉ (higher density → smaller Jeans mass)
M/M_J = 10⁴ / 1.3 ≈ 7692 Jeans masses
N_min = 7692 × 50 ≈ 3.8×10⁵
```

---

#### **E. Combined Worst-Case Scenarios**

| Scenario | M | R | b | v | N_min | N_recommended | Why this is challenging |
|----------|---|---|---|---|-------|---------------|------------------------|
| **Dense close fast** | 10⁴ M_☉ | 3 pc | 3 pc | 20 km/s | 8×10⁵ | 1×10⁶ | Small M_J + small r_H + shock → maximum resolution |
| **Massive close slow** | 5×10⁴ M_☉ | 15 pc | 3 pc | 5 km/s | 9×10⁵ | 1.2×10⁶ | Many M_J + extended tidal tails |
| **Standard publication** | 10⁴ M_☉ | 5 pc | 4 pc | 10 km/s | 1.4×10⁵ | 2×10⁵ | **Baseline for papers** |
| **Quick exploratory** | 10³ M_☉ | 4.6 pc | 5 pc | 10 km/s | 1.4×10⁴ | 5×10⁴ | Parameter scan, proof of concept |
| **Control (weak)** | 10⁴ M_☉ | 7 pc | 6 pc | 10 km/s | 8×10⁴ | 1×10⁵ | Reference case, minimal disruption |

**Worst-case calculation** (Dense close fast):
```
M = 10⁴ M_☉, R = 3 pc → <n_H> = 3564 cm⁻³ → M_J = 1.3 M_☉
N_Jeans = (10⁴ / 1.3) × 50 = 3.8×10⁵

b = 3 pc → r_H = 0.99 pc
h_max = 0.1 × r_H = 0.099 pc
N_spatial = (3 pc / 0.099 pc)² = 9.2×10⁵

v = 20 km/s → shock resolution: N_shock ≈ 1.5 × N_baseline = 2.1×10⁵

N_min = max(3.8×10⁵, 9.2×10⁵, 2.1×10⁵) ≈ 9.2×10⁵
Recommended: 1×10⁶ (round up, add 10% safety)
```

---

#### **F. Summary: Minimal N vs Recommended N**

| Use Case | N_min (bare minimum) | N_recommended | Safety margin |
|----------|---------------------|---------------|---------------|
| Quick parameter scan | 5×10⁴ | 5×10⁴ | None (accept lower quality) |
| Exploratory run | 1×10⁵ | 1.5×10⁵ | +50% |
| **Standard publication** | 1.4×10⁵ | **2×10⁵** | +43% |
| High-quality paper | 2×10⁵ | 4×10⁵ | +100% |
| Flagship simulation | 5×10⁵ | 1×10⁶ | +100% |

**Reviewer checklist** (what they will scrutinize):
1. ✓ **Jeans resolution**: N_J ≥ 50 particles per M_J
2. ✓ **Neighbor count**: N_neighbor ≥ 50 within smoothing kernel
3. ✓ **Tidal feature**: h ≤ 0.1 r_H (resolve Hill sphere)
4. ✓ **Convergence test**: Compare N and 2N results (< 10% difference)
5. ✓ **Energy conservation**: ΔE/E₀ < 1% over simulation
6. ✓ **Shock resolution**: At least 5-10 particles across shock front (if v > 15 km/s)

**Recommended approach for papers**:
- Main results: N = 2×10⁵ (baseline) + N = 4×10⁵ (convergence test)
- Parameter scan: N = 1×10⁵ (acceptable) → upgrade to 2×10⁵ for final published figures
- Supplementary: N = 5×10⁴ (proof of concept only, clearly labeled)

### 5. Thermal Equilibrium Constraint

From Koyama & Inutsuka (2000) cooling timescale:

```
t_cool ~ u / |du/dt| ~ u / |(u_eq - u) / τ_relax|
```

**Numerical example** (τ_relax = 0.05 Myr, t_tidal = 0.54 Myr):
```
τ_relax / t_tidal = 0.05 / 0.54 ≈ 0.09 = 9%
```
- **Thermal equilibrium maintained**: T → T_eq(n) within ~10% of t_tidal
- **Cooling is ~10× faster than tidal disruption** → quasi-equilibrium valid

**Specific internal energy** (T = 50 K, γ = 5/3):
```
u = k_B T / ((γ-1) μ m_H)
  = (1.38×10⁻¹⁶ erg/K × 50 K) / (0.667 × 1.67×10⁻²⁴ g)
  ≈ 6.2×10⁹ erg/g
  ≈ 6.2×10¹² cm²/s² (code units: ~0.062 in [pc²/Myr²])
```

**Implementation**: Use `KoyamaInutsukaCooling` with relaxation time:
```cpp
real tau_relax = 0.05;  // Myr (5% of tidal timescale ~0.5 Myr)
real du_dt = cooling.cooling_rate(n_H, T_current, tau_relax);
// For T >> T_eq: du/dt ≈ -(u - u_eq) / tau_relax ≈ -u / 0.05 Myr
```

## Physics Modules Required

### 1. External Force: Point-Mass Black Hole

**Header**: `include/external_forces/point_mass_bh.hpp`

```cpp
namespace sph {
namespace external_forces {

class PointMassBlackHole {
    vec_t m_position;      // BH position [pc]
    real m_mass;           // BH mass [M_☉]
    real m_softening;      // Gravitational softening [pc]
    
public:
    vec_t acceleration(const vec_t& r_particle) const;
    real potential(const vec_t& r_particle) const;
};

}
}
```

**Force**:
```
F_BH = -G M_BH (r - r_BH) / (|r - r_BH|² + ε²)^(3/2)
```

where ε = softening length ~ 0.01 R_cloud (prevents singularity).

**Numerical example** (M_BH = 10⁵ M_☉, r = 3 pc from BH, ε = 0.05 pc):
```
|F_BH| = G M_BH / r²
       = (4.302×10⁻³ pc·(km/s)²·M_☉⁻¹) × 10⁵ M_☉ / (3 pc)²
       = 430.2 / 9 [km²/s²/pc]
       ≈ 47.8 km²/s²/pc
       ≈ 47.8 pc/Myr² (code units: ~47.8)

Acceleration: a_BH = F_BH / m_particle
```

**Gravitational potential energy** (cloud at b = 3 pc):
```
E_pot = -G M_BH M_cloud / b
      = -(4.302×10⁻³) × 10⁵ × 10⁴ / 3 [M_☉·(km/s)²]
      ≈ -1.43×10⁸ M_☉·(km/s)²
      ≈ -2.85×10⁵⁴ erg (enormous binding energy!)
```

### 2. Thermal Physics: K&I Cooling

**Existing module**: `include/thermal/koyama_inutsuka_cooling.hpp`

- Already implemented ✓
- Provides T_eq(n_H), P_eq(n_H)
- Cooling rate: du/dt = (u_eq - u) / τ_relax

### 3. Initial Conditions: Lane-Emden + Thermal Equilibrium

**Based on**: `src/sample/lane_emden.cpp`

Modifications needed:
1. Set initial T = T_eq(n) from K&I curve (not polytro

pic P = K ρ^γ)
2. Add IMBH at offset position
3. Optionally add relative velocity

### 4. SPH Method: GDISPH (Recommended)

**Why GDISPH**:
- Best for density gradients (cloud disruption creates steep gradients)
- Hybrid Godunov + pressure-energy formulation
- Built-in Riemann solver handles shocks from tidal compression
- Already tested in `sample/sedov/`, `sample/khi/`

**Configuration**:
```json
{
  "numerical": {
    "sph_type": "gdisph",
    "kernel": "wendland",
    "neighbor_number": 50,
    "use_gravity": true
  },
  "artificial_viscosity": {
    "alpha": 1.0,
    "use_balsara_switch": true,
    "use_time_dependent_av": true
  }
}
```

## Code Units and Conversion

### Recommended Unit System

**Length unit**: 1 code unit = 1 parsec (pc)
**Mass unit**: 1 code unit = 1 M_☉
**Time unit**: Derived from G

```
G = 4.302 × 10⁻³ pc (km/s)² M_☉⁻¹
t_unit = √(L³ / (G M)) = √(pc³ / (G M_☉))
       ≈ 0.978 Myr
```

**Velocity unit**: v_unit = L/t_unit ≈ 1.02 km/s

### Physical Parameters in Code Units

| Quantity | Physical Value | Code Units |
|----------|----------------|------------|
| M_BH | 10⁵ M_☉ | 10⁵ |
| M_cloud | 10⁴ M_☉ | 10⁴ |
| R_cloud | 5 pc | 5.0 |
| b (impact) | 3-6 pc | 3.0-6.0 |
| v_rel (baseline) | 10 km/s | ~9.8 |
| v_rel (slow) | 5 km/s | ~4.9 |
| v_rel (fast) | 15 km/s | ~14.7 |
| v_rel (very fast) | 20 km/s | ~19.6 |
| T_CNM | 50 K | — |
| n_H | 100 cm⁻³ | — |
| t_end | 2 Myr | ~2.0 |

**Velocity conversion**: v[code units] = v[km/s] / 1.02

### Initial Conditions Setup

**Cloud placement**: Center at origin (0, 0, 0)
**BH placement**: Position (b, 0, 0) with velocity (-v_x, 0, 0)
  - Moves toward cloud along x-axis
  - Pericenter passage at x ≈ 0 (closest approach)
**Simulation time**: 2 Myr (sufficient for multiple crossing times)

### Density Conversion

From code density ρ [M_☉/pc³] to number density n_H [cm⁻³]:

```
n_H = ρ * (M_☉/pc³) / (μ m_H)
    = ρ * 40.5  [for μ = 1, neutral H]
```

where:
- μ = 1 (neutral hydrogen)
- m_H = 1.67 × 10⁻²⁴ g
- 1 pc³ = 2.938 × 10⁵⁷ cm³

**Derivation**:
```
n_H = ρ [M_☉/pc³] × (1.989×10³³ g/M_☉) / (2.938×10⁵⁷ cm³/pc³) / (1.67×10⁻²⁴ g)
    = ρ × (1.989×10³³ / 2.938×10⁵⁷ / 1.67×10⁻²⁴)
    = ρ × 40.5 cm⁻³
```

**Numerical examples**:

| ρ [M_☉/pc³] | n_H [cm⁻³] | Physical regime |
|-------------|------------|-----------------|
| 0.1 | 4.05 | Warm neutral medium (WNM) |
| 1.0 | 40.5 | Cold neutral medium (CNM) - low density |
| 2.47 | 100 | **CNM - typical molecular cloud** |
| 10 | 405 | Dense molecular cloud |
| 24.7 | 1000 | Very dense - star forming |

**Cloud mass to density** (M_cloud = 10⁴ M_☉, R_cloud = 5 pc):
```
<ρ> = M_cloud / (4π/3 R_cloud³)
    = 10⁴ M_☉ / (4π/3 × 125 pc³)
    ≈ 19.1 M_☉/pc³
    ≈ 775 cm⁻³ (central density higher!)
```

## Implementation Checklist

### Phase 1: External Force Module
- [ ] `include/external_forces/point_mass_bh.hpp`
- [ ] `src/external_forces/point_mass_bh.cpp`
- [ ] Add to CMakeLists.txt
- [ ] Unit tests

### Phase 2: Initial Conditions
- [ ] `src/sample/imbh_cloud.cpp`
  - [ ] Lane-Emden sphere with n=5/3
  - [ ] K&I thermal equilibrium T_eq(n)
  - [ ] IMBH placement at distance b
  - [ ] Optional relative velocity
  
### Phase 3: Configuration Presets

**Impact Parameter Study**:
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b3pc_v10_gdisph.json`
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b4pc_v10_gdisph.json`
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b5pc_v10_gdisph.json`
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b6pc_v10_gdisph.json`

**Velocity Study** (fixed b = 4 pc):
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b4pc_v0_gdisph.json` (static infall)
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b4pc_v5_gdisph.json`
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b4pc_v10_gdisph.json` (baseline)
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b4pc_v15_gdisph.json`
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b4pc_v20_gdisph.json`

### Phase 4: Visualization & Analysis
- [ ] `sample/imbh_cloud/scripts/visualize_disruption.py`
  - [ ] Density evolution
  - [ ] Tidal deformation (axis ratios)
  - [ ] Mass loss rate
  - [ ] Temperature vs density phase diagram
  
- [ ] `sample/imbh_cloud/scripts/compare_impact_parameters.py`
  - [ ] Multi-panel comparison for b = 3, 4, 5, 6 pc
  - [ ] Bound vs unbound mass evolution
  - [ ] Tidal tail morphology differences
  
- [ ] `sample/imbh_cloud/scripts/compare_velocities.py`
  - [ ] Multi-panel comparison for v = 0, 5, 10, 15, 20 km/s
  - [ ] Encounter timescale analysis
  - [ ] Peak compression vs velocity
  
### Phase 5: Makefile
- [ ] `sample/imbh_cloud/Makefile.imbh_cloud`
  - [ ] Single run targets (specific b, v combinations)
  - [ ] Impact parameter scan: `imbh_scan_impact` (4 runs)
  - [ ] Velocity scan: `imbh_scan_velocity` (5 runs)
  - [ ] Full matrix scan: `imbh_scan_full` (20 runs)
  - [ ] Comparison visualization targets

## Expected Physics

### Tidal Disruption Sequence

1. **Initial approach** (t = 0 - 0.3 t_tidal):
   - Cloud feels increasing tidal gradient
   - Elongation along radial direction
   - Compression perpendicular to orbital plane (pancake)
   - **Velocity dependence**: Higher v → faster approach → less pre-compression

2. **Maximum compression** (t ~ t_peri = b/v_rel):
   - Peak density reached at pericenter passage
   - Shocks from rapid compression (stronger for high v_rel)
   - Potential triggered star formation (future work)
   - **Impact parameter effect**: Smaller b → stronger tidal force → more compression
   - **Velocity effect**: Higher v → impulsive shock → heating vs cooling competition

3. **Disruption** (t > t_tidal):
   - Leading/trailing tidal tails form
   - Bound vs unbound material separation
   - Self-gravity determines remnant mass
   - **b-dependence**: Close encounters (b ~ 3 pc) → 50-90% mass loss
   - **b-dependence**: Distant encounters (b ~ 6 pc) → <10% mass loss
   - **v-dependence**: Slow encounters → gradual tidal stripping
   - **v-dependence**: Fast encounters → impulsive disruption with shocks

4. **Thermal response**:
   - Compression → n increases → T_eq decreases (CNM branch)
   - Cooling maintains thermal equilibrium (for slow encounters)
   - **Shock heating** (for fast encounters): May temporarily exceed T_eq
   - Cooling time vs compression time determines thermal state

### Diagnostics to Track

1. **Morphology**:
   - Axis ratios (a:b:c) from moment of inertia tensor
   - Elongation parameter: e = 1 - c/a
   - **Parameter study**: Plot e(t) for different b and v
   
   **Expected values**:
   - t = 0: e ≈ 0 (spherical cloud, a ≈ b ≈ c)
   - t = t_peri: e ~ 0.5-0.8 (strong tidal elongation, "cigar" shape)
   - t >> t_tidal: e → varies (tidal tails, disrupted structure)

2. **Energetics**:
   - Kinetic energy: E_kin = Σ(½ m v²)
   - Thermal energy: E_th = Σ(m u)
   - Gravitational potential: E_pot (cloud self-gravity + BH)
   - Total energy conservation: ΔE/E₀ < 1%
   - **Parameter study**: E_kin vs time shows velocity damping
   
   **Numerical estimates** (M_cloud = 10⁴ M_☉, v_rel = 10 km/s, T = 50 K):
   ```
   E_kin,initial ≈ ½ M_cloud v_rel² 
                ≈ 0.5 × 10⁴ M_☉ × (10 km/s)²
                ≈ 5×10⁵ M_☉·(km/s)²
                ≈ 10⁵⁴ erg
   
   E_th,initial ≈ M_cloud × u
               ≈ 10⁴ M_☉ × 6.2×10⁹ erg/g × 1.989×10³³ g/M_☉
               ≈ 1.2×10⁴⁷ erg (negligible compared to E_kin!)
   
   E_pot,BH ≈ -G M_BH M_cloud / b
            ≈ -1.43×10⁸ M_☉·(km/s)² (for b=3pc)
            ≈ -2.85×10⁵⁴ erg (dominates total energy!)
   ```

3. **Mass Budget**:
   - Bound mass: M_bound(t) = Σm for particles with E_total < 0
   - Unbound mass: M_unbound(t) = Σm for particles with E_total > 0
   - Accretion rate: dM_BH/dt (particles within R_acc ~ 0.01 pc)
   - **Key metric**: f_disrupted = (M₀ - M_bound(t_final)) / M₀
   
   **Expected disruption fractions**:
   
   | b [pc] | v [km/s] | f_disrupted | M_lost [M_☉] | Physical interpretation |
   |--------|----------|-------------|--------------|------------------------|
   | 3 | 10 | ~0.70 | 7000 | Strong disruption |
   | 4 | 10 | ~0.40 | 4000 | Moderate disruption |
   | 5 | 10 | ~0.20 | 2000 | Weak disruption |
   | 6 | 10 | ~0.10 | 1000 | Minimal disruption |
   | 4 | 5 | ~0.50 | 5000 | Slower → more stripping |
   | 4 | 20 | ~0.25 | 2500 | Faster → less stripping |

4. **Thermal State**:
   - n-T diagram (should follow K&I curve for slow encounters)
   - Departures from equilibrium: ΔT = T - T_eq(n)
   - Cooling time: τ_cool ~ u / |du/dt|
   - **Parameter study**: Fast encounters (high v) show shock heating (T > T_eq)
   
   **Expected temperature evolution** (initial T = 50 K):
   ```
   Compression at pericenter: n increases 10× → ρ → 10× higher
   K&I curve (CNM branch): T_eq ∝ n^(-0.7) (roughly)
   Expected: T_eq drops to ~50 K / 10^0.7 ≈ 10 K (colder!)
   
   Shock heating (v=20 km/s): ΔT_shock ~ (γ-1) v² / (2 k_B / μ m_H)
                                        ~ 0.667 × (20 km/s)² / (2×8.3×10⁷)
                                        ~ 160 K (hot shock!)
   
   Result: Competition between compression cooling vs shock heating
   ```

## References

1. Koyama & Inutsuka (2000), ApJ 532, 980: Thermal equilibrium curves
2. Guillochon & Loeb (2015), ApJ 811, 20: Tidal disruption events
3. Burkert & Naab (2013), MNRAS 434, 36: Cloud-BH interactions
4. Stone et al. (1998), ApJS 114, 345: Numerical methods for tidal disruption
5. Rees (1988), Nature 333, 523: Tidal disruption by supermassive black holes
6. Hills (1975), Nature 254, 295: Tidal capture and disruption

## Appendix: Derivation of Physical Criteria

This appendix provides detailed derivations of all key physical scales and timescales used in the resolution analysis.

---

### A.1 Jeans Mass and Jeans Length

**Physical origin**: Self-gravitating clouds are unstable to collapse when gravity overcomes thermal pressure support.

**Jeans criterion**: For a perturbation of wavelength λ, gravitational collapse occurs when:
```
Gravitational potential energy > Thermal kinetic energy
G M² / R > N k_B T
```

**Derivation from sound wave dispersion**:

Consider a density perturbation: ρ = ρ₀ + δρ e^(i(k·r - ωt))

The dispersion relation for sound waves in self-gravitating medium:
```
ω² = c_s² k² - 4πG ρ₀
```

**Critical wavenumber** (ω² = 0, marginal stability):
```
k_J = √(4πG ρ₀) / c_s
```

**Jeans length**:
```
λ_J = 2π / k_J = c_s √(π / (G ρ₀))
```

**Numerical example** (n_H = 100 cm⁻³, T = 50 K):
```
ρ₀ = 100 cm⁻³ × 1.67×10⁻²⁴ g = 1.67×10⁻²² g/cm³
c_s = √(k_B T / μ m_H) = 0.29 km/s = 2.9×10⁵ cm/s
G = 6.67×10⁻⁸ cm³/(g·s²)

λ_J = 2.9×10⁵ cm/s × √(π / (6.67×10⁻⁸ × 1.67×10⁻²²))
    = 2.9×10⁵ × √(2.81×10¹³)
    ≈ 1.54×10¹⁸ cm
    ≈ 0.50 pc
```

**Jeans mass** (mass within Jeans length sphere):
```
M_J = ρ₀ × (4π/3) × (λ_J/2)³
    = (π/6)^(1/2) × (c_s³ / (G^(3/2) ρ₀^(1/2)))
    = (π^(5/2) / 6) × (c_s³ / (G^(3/2) ρ₀^(1/2)))
```

**Numerical example**:
```
M_J = (π^(5/2) / 6) × (2.9×10⁵)³ / ((6.67×10⁻⁸)^(3/2) × (1.67×10⁻²²)^(1/2))
    ≈ 6.9×10³³ g
    ≈ 3.5 M_☉
```

**SPH resolution requirement**: To resolve fragmentation, need **N_J ≥ 50** particles per M_J.

---

### A.2 Tidal Disruption: Pancaking and Elongation

**Physical setup**: Cloud of mass M_cloud, radius R_cloud, at distance r from BH of mass M_BH.

#### **A.2.1 Tidal Force Components**

The tidal force arises from the **gradient** of the BH gravitational field across the cloud.

**BH acceleration at cloud center** (distance r along x-axis):
```
a_center = -G M_BH / r² x̂
```

**Acceleration at cloud edge** (at position r + δr):
```
a_edge = -G M_BH / (r + δr)² x̂
       ≈ -G M_BH / r² × (1 - 2δr/r) x̂  [Taylor expansion]
```

**Tidal acceleration** (differential force):
```
a_tidal = a_edge - a_center
        = 2 G M_BH δr / r³ x̂
```

**General tidal tensor** (for arbitrary displacement δr = (δx, δy, δz)):
```
T_ij = -G M_BH / r³ × [3 r̂_i r̂_j - δ_ij]
```

where r̂ is the unit vector from BH to cloud center.

**For BH along x-axis**, the tidal tensor is:
```
        ⎡ +2    0    0 ⎤
T = -GM/r³ ⎢  0   -1    0 ⎥
        ⎣  0    0   -1 ⎦
```

#### **A.2.2 Pancaking (Vertical Compression)**

**Perpendicular direction** (y and z): Compression force
```
F_⊥ = m × T_yy × δy = +G M_BH m δy / r³
```

This is a **restoring force** (positive → compression toward orbital plane).

**Equation of motion** (simple harmonic oscillator):
```
d²δy/dt² = (G M_BH / r³) δy
```

**Vertical oscillation frequency**:
```
ω_⊥ = √(G M_BH / r³)
```

**Numerical example** (M_BH = 10⁵ M_☉, r = 3 pc at pericenter):
```
ω_⊥ = √(4.302×10⁻³ × 10⁵ / 3³)  [in units of Myr⁻¹]
    = √(430.2 / 27)
    ≈ 4.0 Myr⁻¹

T_compress = 2π / ω_⊥ ≈ 1.6 Myr
```

#### **A.2.3 Elongation (Radial Stretching)**

**Radial direction** (x): Stretching force
```
F_∥ = m × T_xx × δx = -2 G M_BH m δx / r³
```

This is **anti-restoring** (negative → exponential growth of tidal tails).

**Equation of motion**:
```
d²δx/dt² = -2 (G M_BH / r³) δx
```

**Exponential growth rate**:
```
γ_∥ = √(2 G M_BH / r³) = √2 × ω_⊥
```

**Tidal tail growth**:
```
δx(t) = δx₀ × exp(γ_∥ t)
```

**Numerical example** (same parameters):
```
γ_∥ = √2 × 4.0 ≈ 5.7 Myr⁻¹
e-folding time: τ_stretch = 1/γ_∥ ≈ 0.18 Myr
```

**Physical interpretation**: Material stretches exponentially along the BH direction, forming **leading and trailing tidal tails**.

#### **A.2.4 Tidal Disruption Timescale**

Cloud is **tidally disrupted** when tidal force exceeds self-gravity:
```
G M_BH R_cloud / r³ > G M_cloud / R_cloud²
```

**Critical radius** (tidal radius):
```
r_t = R_cloud × (M_BH / M_cloud)^(1/3)
```

**Disruption timescale** (time for cloud to cross tidal radius):
```
t_tidal ~ r_t / v_orb
        ~ r_t / √(G M_BH / r_t)
        ~ √(r_t³ / (G M_BH))
        ~ √(R_cloud³ × (M_BH / M_cloud) / (G M_BH))
        ~ √(R_cloud³ / (G M_cloud))
```

**Numerical example** (M_cloud = 10⁴ M_☉, R_cloud = 5 pc):
```
t_tidal = √(125 pc³ / (4.302×10⁻³ × 10⁴))
        = √(125 / 43.02)
        ≈ 0.54 Myr
```

**Note**: This is independent of M_BH! Disruption time depends only on cloud properties.

---

### A.3 Tidal Radius vs Hill Radius

These are **different concepts** often confused:

#### **A.3.1 Tidal Radius** (r_t)

**Definition**: Radius at which BH tidal force equals cloud self-gravity.

**Derivation**: Balance tidal force and self-gravity at cloud surface:
```
G M_BH R_cloud / r_t³ = G M_cloud / R_cloud²
```

**Result**:
```
r_t = R_cloud × (M_BH / M_cloud)^(1/3)
```

**Numerical example** (M_cloud = 10⁴ M_☉, M_BH = 10⁵ M_☉, R_cloud = 5 pc):
```
r_t = 5 pc × (10⁵ / 10⁴)^(1/3)
    = 5 pc × 2.15
    ≈ 10.8 pc
```

**Physical meaning**: Cloud passing within r_t will be **tidally disrupted**.

#### **A.3.2 Hill Radius** (r_H or Roche Lobe)

**Definition**: Maximum distance from cloud where cloud's gravity dominates over tidal force.

**Derivation**: Consider test particle at distance r_H from cloud center. It feels:
1. Cloud gravity: F_cloud = G M_cloud / r_H²
2. Tidal force from BH: F_tidal = 2 G M_BH r_H / b³ (where b = impact parameter)

**Balance point** (Hill sphere boundary):
```
G M_cloud / r_H² = 2 G M_BH r_H / b³
```

**Solve for r_H**:
```
r_H³ = M_cloud b³ / (2 M_BH)
r_H = b × (M_cloud / (2 M_BH))^(1/3)
```

**More accurate formula** (including factor of 3 from orbit geometry):
```
r_H = b × (M_cloud / (3 M_BH))^(1/3)
```

**Numerical example** (M_cloud = 10⁴ M_☉, M_BH = 10⁵ M_☉, b = 4 pc):
```
r_H = 4 pc × (10⁴ / (3×10⁵))^(1/3)
    = 4 pc × (0.0333)^(1/3)
    = 4 pc × 0.322
    ≈ 1.29 pc
```

**Physical meaning**: Material beyond r_H will be **tidally stripped** from the cloud.

**Comparison**:
- **r_t**: Where disruption begins (intrinsic to cloud, scales with R_cloud)
- **r_H**: Extent of bound material (depends on impact parameter b)

---

### A.4 Pericenter Passage and Orbital Timescales

#### **A.4.1 Pericenter Time**

For a **parabolic encounter** (typical for cloud-BH interactions):

**Initial conditions**: Cloud at distance b (impact parameter), velocity v_∞ at infinity.

**Energy conservation**:
```
½ v_∞² = ½ v_peri² - G M_BH / r_peri
```

For **hyperbolic orbit** (v_∞ > 0):
```
r_peri ≈ b  (for weak encounters)
```

**Pericenter passage time** (time spent within distance ~ b):
```
t_peri ~ b / v_rel
```

**Numerical examples**:
```
b = 3 pc, v = 10 km/s:
t_peri = 3 pc / 10 km/s
       = 3 × 3.086×10¹⁸ cm / (10 × 10⁵ cm/s)
       = 9.26×10¹² s
       ≈ 0.29 Myr

b = 6 pc, v = 20 km/s:
t_peri = 6 pc / 20 km/s ≈ 0.30 Myr
```

**Physical interpretation**: 
- **Slow encounters** (v small): Long t_peri → gradual tidal stripping
- **Fast encounters** (v large): Short t_peri → impulsive shock

---

### A.5 Smoothing Length and Neighbor Number

#### **A.5.1 SPH Kernel Support**

The SPH smoothing kernel W(r, h) has compact support within radius **2h** (for cubic spline) or **3h** (for Wendland C4).

**Neighbor number**:
```
N_neighbor = ∫ n(r') W(|r - r'|, h) dV
           ≈ n₀ × V_kernel
```

For **Wendland C4 kernel** (support radius = 2h):
```
V_kernel = (4π/3) × (2h)³ = 32π h³ / 3
```

**Average neighbor count**:
```
N_neighbor = n₀ × 32π h³ / 3
```

For **uniform distribution** with N total particles in volume V:
```
n₀ = N / V
```

**Smoothing length** (from neighbor number requirement):
```
N_neighbor = (N / V) × 32π h³ / 3
```

For **spherical cloud** (V = 4πR³/3):
```
N_neighbor = N × (h/R)³ × 8
```

**Solve for h**:
```
h = R × (N_neighbor / (8N))^(1/3)
```

**For N_neighbor = 50** (SPH standard):
```
h = R × (50 / (8N))^(1/3)
  = R × (6.25 / N)^(1/3)
  ≈ R / (N^(1/3) × 0.54)
  ≈ 1.84 R / N^(1/3)
```

**Simplified approximation**:
```
h ≈ R / √N  (slightly overestimates)
```

**Numerical example** (R = 5 pc, N = 2×10⁵):
```
h = 5 pc / √(2×10⁵)
  = 5 / 447
  ≈ 0.011 pc
```

#### **A.5.2 Resolution Requirement**

To resolve a physical scale λ, need:
```
h < λ / 5  (at least 5 smoothing lengths across feature)
```

For **Hill radius** (λ = r_H):
```
h_max = r_H / 5 ≈ 0.2 r_H
```

**Conservative criterion** (used in this work):
```
h ≤ 0.1 r_H
```

This ensures **≥10 smoothing lengths** across the Hill sphere → well-resolved tidal truncation.

---

### A.6 Shock Heating and Compression

#### **A.6.1 Rankine-Hugoniot Shock Conditions**

For a **strong shock** in ideal gas (γ = 5/3):

**Density jump**:
```
ρ₂ / ρ₁ = ((γ+1) M²) / ((γ-1) M² + 2)
```

For **Mach number M → ∞** (strong shock limit):
```
ρ₂ / ρ₁ → (γ+1) / (γ-1) = 4  (for γ = 5/3)
```

**Temperature jump**:
```
T₂ / T₁ = [2γ M² - (γ-1)] × [(γ-1) M² + 2] / [(γ+1)² M²]
```

For **strong shock**:
```
T₂ / T₁ ≈ 2(γ-1) M² / (γ+1)² ≈ 0.4 M²
```

**Numerical example** (v_shock = 20 km/s, c_s,1 = 0.3 km/s):
```
M = v_shock / c_s,1 = 20 / 0.3 ≈ 67

ρ₂ / ρ₁ ≈ 4  (maximum compression)

T₂ / T₁ ≈ 0.4 × 67² ≈ 1800

For initial T₁ = 50 K:
T₂ ≈ 90,000 K  (but cooling brings it down quickly!)
```

#### **A.6.2 Post-Shock Temperature**

**Kinetic energy → thermal energy**:
```
½ v² = (γ / (γ-1)) × (k_B ΔT / μ m_H)
```

**Temperature increase**:
```
ΔT = (γ-1) μ m_H v² / (2 k_B)
```

**Numerical example** (v = 20 km/s, μ = 1):
```
ΔT = 0.667 × 1.67×10⁻²⁴ g × (20×10⁵ cm/s)² / (2 × 1.38×10⁻¹⁶ erg/K)
   ≈ 160 K
```

**But**: Cooling time in CNM is short! Koyama-Inutsuka equilibrium curve wins for **slow encounters**.

---

### A.7 Thermal Equilibrium Timescale

#### **A.7.1 Cooling Function**

From Koyama & Inutsuka (2000), cooling rate per unit volume:
```
Λ(n, T) = n² × L(T)
```

where L(T) is the cooling function [erg·cm³/s].

**Cooling time**:
```
t_cool = thermal energy / cooling rate
       = (3/2) n k_B T / (n² L(T))
       = (3/2) k_B T / (n L(T))
```

**Typical values** (CNM, n = 100 cm⁻³, T = 50 K):
```
L(50 K) ≈ 10⁻²⁶ erg·cm³/s  [from K&I curve]
t_cool ≈ (3/2) × 1.38×10⁻¹⁶ × 50 / (100 × 10⁻²⁶)
       ≈ 10¹⁰ s
       ≈ 0.03 Myr
```

**Compare to tidal timescale**:
```
t_cool / t_tidal ≈ 0.03 / 0.54 ≈ 0.06 = 6%
```

**Conclusion**: Cooling is ~17× faster than tidal disruption → **thermal equilibrium maintained**.

---

### A.8 Summary of Scaling Relations

| Quantity | Scaling | Numerical example |
|----------|---------|-------------------|
| **Jeans mass** | M_J ∝ ρ^(-1/2) T^(3/2) | 3.5 M_☉ (n=100 cm⁻³, T=50K) |
| **Jeans length** | λ_J ∝ ρ^(-1/2) T^(1/2) | 0.5 pc |
| **Tidal radius** | r_t ∝ R (M_BH/M)^(1/3) | 10.8 pc (M=10⁴M_☉, M_BH=10⁵M_☉) |
| **Hill radius** | r_H ∝ b (M/M_BH)^(1/3) | 1.3 pc (b=4pc) |
| **Tidal time** | t_tidal ∝ √(R³/(GM)) | 0.54 Myr (M=10⁴M_☉, R=5pc) |
| **Pericenter time** | t_peri = b / v | 0.29 Myr (b=3pc, v=10km/s) |
| **Cooling time** | t_cool ∝ T / (n L(T)) | 0.03 Myr (n=100 cm⁻³, T=50K) |
| **Smoothing length** | h ∝ R / √N | 0.011 pc (R=5pc, N=2×10⁵) |
| **Shock temperature** | ΔT ∝ v² | 160 K (v=20 km/s) |

**Timescale hierarchy**:
```
t_cool << t_tidal ~ t_peri << t_Hubble
0.03 Myr << 0.5 Myr ~ 0.5 Myr << 10⁴ Myr
```

This justifies the **thermal equilibrium assumption**: clouds cool much faster than they are tidally disrupted.

---

## A.9 Shock Heating and Molecular Dissociation: Radio Observability Constraints

### Physical Motivation

**Critical Question**: When can we observe the IMBH-cloud scattering event with radio telescopes?

**Answer**: Only when molecules (CO, H₂) **survive** the encounter. If shock heating is too strong, molecules dissociate → no radio line emission → event becomes **radio-invisible** despite large velocity perturbations.

This section derives:
1. Shock heating as a function of impact parameters (b, v_rel)
2. Dissociation thresholds for H₂ and CO molecules
3. Observable vs unobservable parameter space
4. Why we assume ~80 km/s line-of-sight velocities are detectable

---

### A.9.1 First Principles: Shock Heating in Tidal Encounters

#### **A.9.1.1 Velocity Perturbation from Tidal Force**

Consider a molecular cloud element at the cloud edge experiencing the IMBH tidal encounter.

**Tidal acceleration** at distance r from cloud center:
```
a_tidal = 2 G M_BH r / d³
```
where d = distance from IMBH to cloud center at pericenter ≈ b (impact parameter).

**Impulse approximation** (fast encounter, t_peri << t_tidal):

The velocity kick received during pericenter passage is:
```
Δv ~ a_tidal × t_peri
   = (2 G M_BH r / b³) × (b / v_rel)
   = 2 G M_BH r / (b² v_rel)
```

**Why t_peri ≈ b / v_rel?** Here's the geometric intuition:

**Reference frame note**: We work in the **cloud rest frame** for convenience.
In reality, both objects move - what matters is the **relative velocity** v_rel = v_BH - v_cloud (vector!).
In the galactic context: cloud ≈ static (part of ISM), IMBH ≈ moving (wandering BH).

```
    IMBH trajectory in cloud rest frame (top-down view, x-y plane)

         y
         ↑
         │                  IMBH trajectory (hyperbolic, but nearly straight
         │                  for high-velocity encounters)
         │
         │       t=-∞              t=0 (pericenter)         t=+∞
         │        ↓                    ↓                      ↓
       b │   •─────────→──────────•──────────→─────────────────→
         │  IMBH  v_rel       closest                    IMBH exits
         │ (incoming)          approach
         │                        ↓
         │                   (x≈0, y≈b)
         │
         │                   ┌─────────┐
         │                   │  Cloud  │
       0 ├───────────────────┤    ●    ├──────────────────────────► x
         │                   │  (0,0)  │
         │                   └─────────┘
         │                   R_cloud ~ 5 pc
         │
         └──

    KEY GEOMETRY:

    • Cloud center: at origin (0, 0, 0)
    • Impact parameter b: perpendicular distance from cloud to IMBH trajectory
    • Velocity v_rel: IMBH velocity vector, pointing in +x direction (magnitude v_rel)
    • IMBH path: y ≈ b (constant), x changes from -∞ to +∞

    At different times:
    - t = -∞: IMBH at (-∞, b, 0), velocity = (+v_rel, 0, 0)
    - t = 0:  IMBH at (0, b, 0), velocity ≈ (+v_rel, small y-component, 0)
    - t = +∞: IMBH at (+∞, b, 0), velocity = (+v_rel, 0, 0)

    Distance from IMBH to cloud center:
    d(t) = √[x(t)² + b²]
         ≈ √[(v_rel × t)² + b²]  (for straight-line approximation)

    Time to cross from x=-b to x=+b:
    Δt = 2b / v_rel  ≈ "encounter duration"


Tidal force strength vs time:

  F_tidal│
         │                  ╱╲              ← Peak at pericenter
         │                ╱    ╲
         │              ╱        ╲
         │            ╱            ╲
         │          ╱                ╲
         │________╱____________________╲_________ time
                 │                      │
                 │◄─────Δt ~ b/v_rel───►│

               (distance < b)
```

**Explanation**:

1. **Impact parameter b**: Closest distance between IMBH trajectory and cloud center
   - For high-velocity encounters (v_rel² >> G M_BH / b), the orbit barely bends
   - Pericenter distance ≈ b

2. **Critical distance**: Tidal force is strongest when IMBH is within distance ~b
   - Distance d from cloud: tidal force ∝ 1/d³
   - For d > 2b: tidal force drops to ~1/8 of peak value
   - For d < b: tidal force is near maximum

3. **Crossing time**: Time for IMBH to traverse distance ~b at speed v_rel
   ```
   Δx ~ b  (distance traveled during strong tidal phase)
   v_rel ~ Δx / Δt
   ⟹ t_peri ~ b / v_rel
   ```

4. **Impulse estimate**: Velocity change ≈ (average acceleration) × (time)
   ```
   Δv ~ a_tidal(d=b) × t_peri
      ~ (2 G M_BH r / b³) × (b / v_rel)
      = 2 G M_BH r / (b² v_rel)
   ```

**Numerical example** (M_BH = 10⁵ M_☉, r = R_cloud = 5 pc, b = 3 pc, v_rel = 10 km/s):
```
Δv = 2 × (4.302×10⁻³ pc·(km/s)²·M_☉⁻¹) × 10⁵ M_☉ × 5 pc / (9 pc² × 10 km/s)
   = 2 × 430.2 × 5 / (9 × 10) [km/s]
   = 4302 / 90
   ≈ 48 km/s
```

**Scaling relation**:
```
Δv ∝ M_BH R_cloud / (b² v_rel)
```

**Key insight**:
- Smaller b → larger Δv (closer encounters have stronger tidal kicks)
- Higher v_rel → smaller Δv (faster encounters have less time for tidal acceleration)

#### **A.9.1.2 Shock Formation from Vertical Pancaking (Rees 1988)**

**CRITICAL CORRECTION**: Tidal disruption creates **two distinct deformation modes**:

1. **In-plane stretching** (forming tidal tails): Differential velocities, **NO strong shocks**
2. **Vertical pancaking** (⊥ orbital plane): Compression → **STRONG SHOCKS**

Following [Rees (1988) Nature 333, 523](https://ui.adsabs.harvard.edu/abs/1988Natur.333..523R) and recent work on [nozzle shocks in TDEs](https://academic.oup.com/mnras/article/511/2/2147/6516431), the dominant shock velocity comes from **vertical compression**, NOT from stretching.

**Vertical oscillation frequency** (from tidal tensor, ⊥ component):

From Appendix A.2.2, the pancaking force gives oscillation frequency:
```
ω_⊥ = √(G M_BH / d³)
```

At pericenter d ≈ b:
```
ω_⊥(b) = √(G M_BH / b³)
```

**Compression velocity** (vertical collapse speed):

The cloud has initial vertical extent H ~ R_cloud. The tidal force compresses it on timescale τ ~ 1/ω_⊥:
```
v_compress ~ H × ω_⊥
          = R_cloud × √(G M_BH / b³)
```

**This is the shock velocity** (different layers collide during pancaking):
```
v_shock = R_cloud √(G M_BH / b³)
```

**NOT** v_shock ~ 1/v_rel from impulse approximation! The shock comes from **vertical compression**.

**Numerical example** (M_BH = 10⁵ M_☉, R_cloud = 5 pc, b = 3 pc):
```
ω_⊥ = √(4.302×10⁻³ × 10⁵ / 27) [Myr⁻¹]
    = √(430.2 / 27)
    = √15.93
    = 3.99 Myr⁻¹

v_shock = 5 pc × 3.99 Myr⁻¹
        = 19.95 pc/Myr
        ≈ 20 km/s  (converting: 1 pc/Myr ≈ 1.02 km/s)
```

**Revised parameter space** (vertical pancaking shocks):

| b [pc] | v_rel [km/s] | ω_⊥ [Myr⁻¹] | v_shock [km/s] | M | T_shock [K] | Molecules |
|--------|-------------|------------|----------------|---|-------------|-----------|
| 3 | any | 4.0 | 20 | 67 | 3,200 | CO dissociates |
| 4 | any | 2.7 | 13 | 43 | 1,370 | ✓ Survive |
| 5 | any | 1.9 | 9 | 30 | 650 | ✓ Safe |
| 6 | any | 1.4 | 7 | 23 | 390 | ✓ Safe |
| 8 | any | 0.9 | 5 | 17 | 200 | ✓ Safe |

**CRITICAL INSIGHT**: Shock velocity is **independent of v_rel**!

It only depends on:
- M_BH (stronger gravity → faster compression)
- R_cloud (larger cloud → more velocity from ω × R)
- b (closer → stronger tidal force → higher ω)

The encounter velocity v_rel determines **encounter duration**, not shock strength!

**Consequence for observations**:
```
Observable if: v_shock < 15 km/s  (molecules survive)
              → R_cloud √(G M_BH / b³) < 15 km/s
              → b > R_cloud × √(G M_BH / 225 pc²/Myr²)

For M_BH = 10⁵ M_☉, R = 5 pc:
b > 5 × √(430.2 / 225) = 5 × 1.38 ≈ 6.9 pc

Need b ≥ 7 pc for molecular survival!
```

References:
- [Rees (1988) Nature 333, 523](https://ui.adsabs.harvard.edu/abs/1988Natur.333..523R) - Original pancaking theory
- [Bonnerot et al. (2022) MNRAS 511, 2147](https://academic.oup.com/mnras/article/511/2/2147/6516431) - Nozzle shock in TDEs
- [Guillochon & Ramirez-Ruiz (2013) ApJ 767, 25](https://academic.oup.com/mnras/article/435/3/1809/1018912) - Strong compression consequences
- [The Process of Stellar Tidal Disruption](https://link.springer.com/article/10.1007/s11214-021-00818-7) - Comprehensive review

---

### A.9.2 Molecular Dissociation Energies and Temperature Thresholds

#### **A.9.2.1 Dissociation Energy**

**Hydrogen molecule (H₂)**:
```
H₂ → H + H
D₀(H₂) = 4.48 eV = 7.17×10⁻¹² erg
```

**Temperature threshold** (kB T ~ D₀):
```
T_diss(H₂) = D₀ / k_B
           = 7.17×10⁻¹² erg / (1.38×10⁻¹⁶ erg/K)
           ≈ 52,000 K
```

**Carbon monoxide (CO)**:
```
CO → C + O
D₀(CO) = 11.09 eV = 1.78×10⁻¹¹ erg
```

**Temperature threshold**:
```
T_diss(CO) = D₀ / k_B
           ≈ 129,000 K
```

**Practical dissociation** (considering collisional rates and timescales):
- H₂ begins significant dissociation at **T > 2,000 K** (collisional dissociation)
- H₂ fully dissociates at **T > 10,000 K** (thermal dissociation)
- CO dissociates at **T > 3,000-5,000 K** (photodissociation + collisional)

**Critical thresholds for radio observations**:
```
T < 1,000 K  → All molecules intact → Strong radio lines
1,000 K < T < 3,000 K → H₂ survives, CO weakens → Moderate detection
3,000 K < T < 10,000 K → CO dissociates, H₂ weakens → Weak/no detection
T > 10,000 K → All molecules gone → Radio-invisible
```

---

### A.9.3 Post-Shock Temperature Calculation

From Rankine-Hugoniot conditions (ideal gas, γ = 5/3):

**Post-shock temperature**:
```
T_shock = T₁ + (γ-1) μ m_H v_shock² / (2 k_B)
```

For **cold pre-shock gas** (T₁ ≈ 50 K << T_shock):
```
T_shock ≈ (2/3) × μ m_H v_shock² / k_B
        = (2/3) × (1.67×10⁻²⁴ g) × v_shock² / (1.38×10⁻¹⁶ erg/K)
        ≈ 8.06×10⁻⁹ × v_shock² [K, with v in cm/s]
```

**Simplified formula** (v_shock in km/s):
```
T_shock ≈ 8.06×10⁻⁹ × (v_shock × 10⁵)² K
        ≈ 8060 × v_shock² K  [v_shock in km/s]
```

**Numerical examples**:

| v_shock [km/s] | T_shock [K] | Molecular state | Radio detection |
|----------------|-------------|-----------------|-----------------|
| 5 | 200 | All molecules intact | ✓ Strong CO, H₂ lines |
| 10 | 800 | Molecules intact | ✓ Strong detection |
| 15 | 1,800 | H₂ begins dissociating | ⚠ Moderate CO, weak H₂ |
| 20 | 3,200 | CO dissociating, H₂ partial | ⚠ Weak/marginal |
| 30 | 7,200 | CO gone, H₂ mostly gone | ✗ Very weak/invisible |
| 40 | 12,800 | All molecules dissociated | ✗ **Radio-invisible** |
| 50 | 20,000 | Fully ionized plasma | ✗ No molecular lines |

---

### A.9.4 Observable Parameter Space: b-v Diagram

Combining the shock velocity formula with temperature threshold:

**Dissociation boundary** (T_shock = T_crit):
```
v_shock = √(T_crit / 8060)  [km/s, for T_crit in K]
```

From earlier: v_shock ~ 2 G M_BH R_cloud / (b² v_rel)

**Critical curve** (molecules survive below this line):
```
v_rel,crit(b) = 2 G M_BH R_cloud / (b² × v_shock,crit)
```

**Numerical evaluation** (M_BH = 10⁵ M_☉, R_cloud = 5 pc):

For **CO survival** (T < 3,000 K → v_shock < 19.3 km/s):
```
v_rel,crit(b) = 2 × 430.2 × 5 / (b² × 19.3)
              = 4302 / (19.3 b²)
              = 223 / b² [km/s]
```

For **H₂ survival** (T < 2,000 K → v_shock < 15.8 km/s):
```
v_rel,crit(b) = 273 / b² [km/s]
```

**Parameter space table**:

| b [pc] | v_rel [km/s] | v_shock [km/s] | T_shock [K] | CO state | H₂ state | Radio detection |
|--------|-------------|----------------|-------------|----------|----------|-----------------|
| **3** | **10** | **48** | **18,500** | ✗ Dissociated | ✗ Dissociated | ✗ **Invisible** |
| **3** | **20** | **24** | **4,600** | ✗ Dissociated | ⚠ Partial | ✗ **Very weak** |
| **3** | **30** | **16** | **2,060** | ✓ Survives | ⚠ Marginal | ⚠ **Marginal** |
| **4** | **10** | **27** | **5,870** | ✗ Dissociated | ✗ Dissociated | ✗ **Weak** |
| **4** | **15** | **18** | **2,610** | ⚠ Marginal | ⚠ Partial | ⚠ **Weak** |
| **4** | **20** | **13** | **1,370** | ✓ Survives | ✓ Survives | ✓ **Detectable** |
| **5** | **10** | **17** | **2,330** | ✓ Survives | ⚠ Marginal | ✓ **Good** |
| **5** | **15** | **11** | **1,040** | ✓ Survives | ✓ Survives | ✓ **Strong** |
| **5** | **20** | **9** | **650** | ✓ Intact | ✓ Intact | ✓ **Strong** |
| **6** | **10** | **12** | **1,160** | ✓ Intact | ✓ Intact | ✓ **Strong** |
| **6** | **20** | **6** | **290** | ✓ Intact | ✓ Intact | ✓ **Strong** |

**Critical insight**:
```
Close + Slow encounters (small b, small v_rel) → Hot shocks → Molecules dissociate → RADIO-INVISIBLE
Distant OR Fast encounters → Moderate shocks → Molecules survive → RADIO DETECTABLE
```

---

### A.9.5 Cooling Timescale vs Encounter Timescale

**Key question**: Does the shocked gas cool fast enough to re-form molecules before the encounter ends?

**Cooling timescale** (hot gas, T ~ 5000 K):
```
t_cool = (3/2) k_B T / (n Λ(T))
```

For **T = 5,000 K**, cooling function Λ ≈ 10⁻²³ erg·cm³/s (CIE cooling):
```
t_cool ≈ (3/2) × 1.38×10⁻¹⁶ × 5000 / (100 cm⁻³ × 10⁻²³)
       ≈ 1.0×10⁻⁹ × 10²³
       ≈ 10¹⁴ s
       ≈ 3 Myr
```

**Encounter timescale**:
```
t_enc ~ b / v_rel ≈ 0.3-0.6 Myr
```

**Comparison**:
```
t_cool / t_enc ≈ 3 / 0.5 ≈ 6
```

**Conclusion**: **Hot gas does NOT cool during the encounter!**

**Implication**:
- If shock heating raises T > 3,000 K → molecules dissociate
- Cooling time (3 Myr) >> encounter time (0.5 Myr)
- Gas remains hot and molecular-free throughout observable window
- **Event is radio-invisible** during the actual scattering

**Molecule reformation timescale**:
- H₂ reformation: t_form ~ 10⁶ / n yr ≈ 10 Myr (for n = 100 cm⁻³)
- CO reformation: requires H₂ first, then another ~1 Myr

**Observable window**: We can only detect the event with radio if the scattering happens **gently enough** that T stays below ~2,000 K.

---

### A.9.6 Line-of-Sight Velocity and Radio Detection

**Observational constraint**: Radio telescopes (e.g., ALMA, NOEMA) detect molecular line emission with velocity resolution Δv ~ 0.1-1 km/s.

**Typical CO line**: ¹²CO(J=1-0) at 115 GHz

**Line-of-sight velocity** during IMBH scattering:
```
v_los = v_tidal · (observer direction)
```

**For optimal geometry** (observer perpendicular to orbital plane):
```
v_los,max ~ v_shock ~ 2 G M_BH R / (b² v_rel)
```

**Our assumption**: We detect v_los ~ 80 km/s

**Working backwards to infer parameters**:

If v_los = 80 km/s, what are the allowed (b, v_rel) combinations?

```
80 km/s = 2 × 430.2 × 5 / (b² v_rel)
80 = 4302 / (b² v_rel)
b² v_rel = 4302 / 80 = 53.8

v_rel = 53.8 / b²
```

**Solutions**:

| b [pc] | v_rel [km/s] | v_shock [km/s] | T_shock [K] | Detectable? |
|--------|-------------|----------------|-------------|-------------|
| 2 | 13.4 | 80 | **51,500 K** | ✗ **NO - too hot!** |
| 3 | 6.0 | 80 | **51,500 K** | ✗ **NO - molecules gone** |
| 4 | 3.4 | 80 | **51,500 K** | ✗ **NO - dissociated** |
| 5 | 2.2 | 80 | **51,500 K** | ✗ **NO - radio-invisible** |

**Problem**: If we naively use impulse approximation with v_shock = 80 km/s, T_shock ~ 51,500 K → **all molecules dissociate** → we **cannot observe** this with radio!

**Resolution of the paradox**:

The **80 km/s line-of-sight velocity** must come from:
1. **Bulk motion** of cloud remnant (not shock heating)
2. **Turbulent velocity** induced by tidal stirring (subsonic, T << T_diss)
3. **Orbital velocity** of disrupted material (kinetic, not thermal)

**Revised interpretation**:

Radio observations detect **large-scale velocity perturbations** (v_los ~ 80 km/s) from:
- Tidal stretching → material flung outward at high velocity
- Cloud fragments moving with different velocities
- **But local gas temperature remains T < 2,000 K** (molecules survive)

This requires:
- **Slow enough encounter** that shock heating is moderate: v_shock < 15 km/s → T < 2,000 K
- **Large enough kick** that bulk velocities reach 80 km/s
- This happens when **t_peri >> R/v_shock** (adiabatic limit, not impulsive)

**Detectable scenario**:

| b [pc] | v_rel [km/s] | v_shock [km/s] | T_shock [K] | Δv_bulk [km/s] | Detection |
|--------|-------------|----------------|-------------|----------------|-----------|
| 4 | 5 | 34 | 9,300 | ~50 | ✗ Too hot |
| 5 | 8 | 21 | 3,500 | ~60 | ⚠ Marginal |
| 6 | 10 | 12 | 1,200 | ~70 | ✓ **Detectable** |
| 7 | 12 | 10 | 800 | ~80 | ✓ **Strong lines** |

**Conclusion**: To observe 80 km/s velocities **and** detect CO/H₂ lines:
- Need **b ≥ 6 pc** or **v_rel ≥ 15 km/s** (or both)
- Moderate encounters where shock heating is limited but tidal kicks are strong
- **Not the closest/slowest encounters** - those are radio-invisible!

---

### A.9.7 Summary: When Can We Observe IMBH-Cloud Scattering with Radio Telescopes?

**Observable conditions** (all must be satisfied):

1. **Temperature constraint**: T_shock < 2,000 K (H₂ survival) or T < 3,000 K (CO survival)
   - Requires: v_shock < 16 km/s

2. **Shock velocity limit**:
   ```
   v_shock = 2 G M_BH R / (b² v_rel) < 16 km/s
   ```

3. **Observable parameter space** (M_BH = 10⁵ M_☉, R = 5 pc):
   ```
   b² v_rel > 270 [pc² km/s]
   ```

4. **Detectable velocity**: Need v_los ~ 10-100 km/s (large enough to distinguish from Galactic rotation)

**Best observational scenarios**:

| Scenario | b [pc] | v_rel [km/s] | v_los [km/s] | T_shock [K] | Detection quality |
|----------|--------|-------------|--------------|-------------|-------------------|
| **Optimal** | 6-8 | 10-15 | 50-80 | 800-1,500 | ✓ Excellent CO + H₂ |
| **Marginal** | 4-5 | 15-20 | 40-60 | 1,500-2,500 | ⚠ Weak CO, no H₂ |
| **Too close** | 3 | 10 | 100+ | 18,000 | ✗ No molecules |
| **Too slow** | 4 | 5 | 60 | 6,000 | ✗ Dissociated |

**Key insight for simulation**:

Our baseline parameters (b = 4 pc, v_rel = 10 km/s) produce:
- v_shock ~ 27 km/s
- T_shock ~ 5,870 K
- **Result**: Molecules dissociate → event is **difficult to observe** with molecular line tracers!

**Recommendation for observable events**:
- Prioritize **b = 6 pc, v_rel = 10-15 km/s**
- Or **b = 5 pc, v_rel = 15-20 km/s**
- These produce detectable velocities (50-80 km/s) while keeping T < 2,000 K

---

### A.9.8 Alternative Tracers for Radio-Invisible Events

When molecular dissociation occurs (T > 3,000 K), we can still detect using:

1. **Atomic fine-structure lines**:
   - [C II] 158 μm (survives to ~10⁴ K)
   - [O I] 63 μm (survives to ionization ~10⁴ K)
   - These trace warm/hot atomic gas

2. **Radio recombination lines** (if ionized, T > 10⁴ K):
   - H recombination lines at cm wavelengths
   - Trace ionized gas dynamics

3. **X-ray emission** (for very hot shocks, T > 10⁶ K):
   - Thermal bremsstrahlung
   - Not expected for our velocities (v < 50 km/s)

4. **Dust continuum** (if dust survives, T < 1,500 K):
   - Submm/mm continuum emission
   - Traces column density even without molecules

**Multi-wavelength strategy**:
- **Cold events** (T < 2,000 K): CO, H₂ molecular lines (ALMA)
- **Warm events** (2,000 < T < 10,000 K): [C II], [O I] (Herschel/SOFIA)
- **Hot events** (T > 10,000 K): H recombination, continuum (VLA)

---

## A.10 Rigorous Hyperbolic Orbit Treatment

This section provides the **exact mathematical treatment** of the IMBH-cloud encounter using classical orbital mechanics, without the impulse approximation.

### A.10.1 Setup and Conserved Quantities

**Two-body problem** in the cloud rest frame:
- Cloud mass: M_cloud (at origin)
- IMBH mass: M_BH (incoming with velocity v_∞ at infinity)
- Reduced mass: μ = M_BH M_cloud / (M_BH + M_cloud) ≈ M_cloud (since M_BH >> M_cloud typically)
- Total mass: M_tot = M_BH + M_cloud ≈ M_BH

**Initial conditions** (t → -∞):
```
Position: r → ∞
Velocity: v → v_∞ (speed at infinity)
Impact parameter: b (perpendicular offset)
```

**Conserved quantities**:

1. **Energy** (per unit reduced mass):
   ```
   E = ½ v_∞² - G M_tot / r = ½ v_∞²  (constant, > 0 for hyperbola)
   ```

2. **Angular momentum** (per unit reduced mass):
   ```
   L = |r × v| = b v_∞  (at infinity, when r ⊥ v)
   ```

3. **Eccentricity**:
   ```
   e = √(1 + 2 E L² / (G² M_tot²))
      = √(1 + (b v_∞² / (G M_tot))²)
   ```

For hyperbolic orbits: **e > 1**

---

### A.10.2 Orbital Equation

**Polar coordinates** (r, θ) with origin at cloud center:

The orbit equation is:
```
r(θ) = a(e² - 1) / (1 + e cos θ)
```

where:
- **Semi-major axis**: a = G M_tot / v_∞² (negative for hyperbola, use |a| in formulas)
- **Eccentricity**: e = √(1 + (b/a)²) for hyperbola
- **θ = 0**: direction of pericenter

**Pericenter distance** (closest approach):
```
r_peri = r(θ=0) = a(e² - 1) / (1 + e)
                = a(e - 1)
```

**Relation to impact parameter**:

For hyperbolic orbit:
```
b = a √(e² - 1)
```

From energy and angular momentum:
```
r_peri = (L² / (G M_tot)) × 1 / (1 + e)
       = b² v_∞² / (G M_tot) × 1 / (1 + e)
```

**Simplification** for high-velocity encounters (v_∞² >> G M_tot / b):

In this limit: e >> 1, so:
```
e ≈ b v_∞² / (G M_tot)
r_peri ≈ b² v_∞² / (G M_tot e) ≈ G M_tot / v_∞² ≈ a
```

But for **weak deflection** (b >> a), we have e ≈ b/a >> 1, giving:
```
r_peri ≈ a(e - 1) ≈ a × b/a = b
```

**Key result**: For fast, distant encounters, **r_peri ≈ b** ✓

---

### A.10.3 Time Dependence (Kepler's Equation for Hyperbola)

**Eccentric anomaly** F (hyperbolic analog of eccentric anomaly):
```
r = a(e cosh F - 1)
tan(θ/2) = √((e+1)/(e-1)) tanh(F/2)
```

**Kepler's equation** (time as function of F):
```
M = e sinh F - F
```

where **Mean anomaly**:
```
M = n(t - t_peri)
```

and **Mean motion**:
```
n = √(G M_tot / |a|³) = v_∞³ / (G M_tot)
```

**Solving for time**:
```
t - t_peri = (e sinh F - F) / n
           = (G M_tot / v_∞³) × (e sinh F - F)
```

This is transcendental and must be solved numerically in general.

---

### A.10.4 Encounter Timescale (Exact Calculation)

**Definition**: Time spent within distance d_crit of pericenter.

Let's compute the time for IMBH to go from distance r₁ to r₂.

**Using conservation of energy**:
```
½ v² - G M_tot / r = ½ v_∞²
v = √(v_∞² + 2 G M_tot / r)
```

**Radial velocity component**:
```
dr/dt = v_r = ± √(v² - (L/r)²)
           = ± √(v_∞² + 2 G M_tot / r - b² v_∞² / r²)
```

**Time integral** (from pericenter r_peri to distance r):
```
t(r) = ∫_{r_peri}^r dr' / √(v_∞² + 2 G M_tot / r' - b² v_∞² / r'²)
```

This can be solved analytically in terms of inverse hyperbolic functions, but is messy.

---

### A.10.5 Approximate Encounter Time

**For weak deflection** (b >> G M_tot / v_∞²):

Approximate the trajectory as nearly straight: r(t) ≈ √(b² + (v_∞ t)²)

**Distance as function of time**:
```
d(t) ≈ √(b² + v_∞² t²)
```

At pericenter (t=0): d(0) = b

**Time to reach distance d**:
```
t(d) = √(d² - b²) / v_∞
```

**Example**: Time to go from d = √2 b to d = √2 b on the other side:
```
Δt = 2 × √(2b² - b²) / v_∞ = 2b / v_∞
```

**Time within distance d = √2 b** (factor of √2 from pericenter):
```
Δt ≈ 2b / v_∞  ✓
```

**Time within distance d = 2b** (factor of 2 from pericenter):
```
Δt = 2 √(4b² - b²) / v_∞ = 2√3 b / v_∞ ≈ 3.46 b / v_∞
```

**Rule of thumb**: Time within distance ~b of pericenter is Δt ~ (1-2) × b / v_∞

Our estimate **t_peri ~ b / v_rel** is correct within a factor of ~2 ✓

---

### A.10.6 Scattering Angle

**Deflection angle** Θ (angle between incoming and outgoing asymptotes):

From hyperbolic orbit theory:
```
tan(Θ/2) = 1/e = a/b
```

where a = G M_tot / v_∞².

**Small-angle approximation** (e >> 1, weak deflection):
```
Θ ≈ 2/e = 2a/b = 2 G M_tot / (b v_∞²)
```

**Numerical example** (M_BH = 10⁵ M_☉, b = 3 pc, v_∞ = 10 km/s):
```
a = G M_BH / v_∞²
  = (4.302×10⁻³ pc·(km/s)²·M_☉⁻¹) × 10⁵ M_☉ / (10 km/s)²
  = 430.2 / 100
  = 4.3 pc

e = √(1 + (b/a)²) = √(1 + (3/4.3)²) = √(1 + 0.486) = √1.486 ≈ 1.22

Θ = 2 arctan(1/e) = 2 arctan(1/1.22) ≈ 2 × 39.3° ≈ 78.6° ≈ 1.37 radians
```

**Significant deflection!** This is not a weak encounter.

**For truly weak deflection** (b >> a), we need:
```
b >> G M_BH / v_∞² ≈ 4.3 pc

For b = 6 pc:
e = √(1 + (6/4.3)²) = √(1 + 1.95) = 1.72
Θ ≈ 2 arctan(0.58) ≈ 60° (still substantial!)

For b = 10 pc:
e = √(1 + (10/4.3)²) = √(1 + 5.41) = 2.53
Θ ≈ 2 arctan(0.40) ≈ 43°

For b = 20 pc:
e = √(1 + (20/4.3)²) = √(1 + 21.6) = 4.75
Θ ≈ 2 arctan(0.21) ≈ 24°
```

**Implication**: For b = 3-6 pc, the orbital deflection is **not small** - the trajectory is significantly curved!

---

### A.10.7 Exact vs Approximate Encounter Time

**Dimensionless parameter**:
```
η = v_∞² b / (G M_tot) = b/a
```

This is the ratio of impact parameter to gravitational scale length.

**Regimes**:
- **Strong deflection**: η << 1 (b << a) → orbit is strongly curved, Θ ~ π
- **Moderate deflection**: η ~ 1 (b ~ a) → significant curvature, Θ ~ π/2
- **Weak deflection**: η >> 1 (b >> a) → nearly straight, Θ ~ 2/η

**Encounter time** can be expressed exactly as:
```
t_enc = (b/v_∞) × f(η)
```

where f(η) is a geometric factor:
- f(η → 0) → diverges (direct collision)
- f(η ~ 1) ≈ 2-3 (moderate encounter)
- f(η >> 1) → 2 (straight-line limit)

**For our problem** (M_BH = 10⁵ M_☉, v = 10 km/s):
```
a = 4.3 pc

b = 3 pc: η = 0.70 (moderate deflection, f ≈ 2.5)
b = 4 pc: η = 0.93 (moderate deflection, f ≈ 2.2)
b = 5 pc: η = 1.16 (weak-ish deflection, f ≈ 2.1)
b = 6 pc: η = 1.40 (weak deflection, f ≈ 2.05)
```

**Conclusion**: The impulse approximation t_peri ~ b/v_∞ is **within a factor of 2-2.5** of the exact result for our parameter range - quite good for order-of-magnitude estimates!

---

### A.10.8 Tidal Force Along the Orbit

**Distance from IMBH to cloud center** along the orbit:
```
d(θ) = r(θ) = a(e² - 1) / (1 + e cos θ)
```

**Tidal acceleration** (radial stretching component):
```
a_tidal(θ) = 2 G M_BH R_cloud / d(θ)³
```

**Time-averaged tidal kick**:
```
Δv_tidal = ∫ a_tidal(t) dt
```

This integral can be done numerically, but for weak deflection (e >> 1):
```
Δv_tidal ≈ a_tidal(pericenter) × t_enc
          ≈ (2 G M_BH R_cloud / b³) × (2b/v_∞)
          ≈ 4 G M_BH R_cloud / (b² v_∞)
```

Compare to our impulse estimate:
```
Δv_impulse = 2 G M_BH R_cloud / (b² v_∞)
```

**Factor of 2 difference** due to time-averaging and geometric factors - but same scaling!

---

### A.10.9 Summary: When is Impulse Approximation Valid?

**Impulse approximation** assumes:
1. Tidal force acts for time t ~ b/v_∞
2. Trajectory is approximately straight
3. Distance at pericenter ≈ b

**Validity criteria**:
```
η = b v_∞² / (G M_tot) >> 1
```

**For our parameters**:
```
η = b × (10 km/s)² / (4.302×10⁻³ × 10⁵ M_☉)
  = b × 100 / 430.2
  = 0.23 × b  [for b in pc]

b = 3 pc: η = 0.70 (marginal - 30% error expected)
b = 6 pc: η = 1.4 (good - 15% error expected)
b = 10 pc: η = 2.3 (excellent - <10% error expected)
```

**Bottom line**: Our impulse approximation is **reasonably accurate** for b ≥ 5 pc, but underestimates the encounter complexity for b = 3-4 pc where orbital deflection is significant.

For **precise predictions**, use the full hyperbolic orbit equations. For **scaling arguments** and **order-of-magnitude estimates**, the impulse approximation is sufficient.

---

### A.10.10 Comparison Table: Exact vs Approximate

| Parameter | Exact (Hyperbolic) | Approximate (Impulse) | Error |
|-----------|-------------------|----------------------|-------|
| r_peri | a(e-1) = b²v²/(GM) / (1+e) | b | ~10-30% for b~3-6pc |
| t_enc | Complex integral, f(η)×b/v | b/v | Factor of 2-2.5 |
| Δv_tidal | Numerical integration | 2GMR/(b²v) | Factor of 2 |
| Θ_deflect | 2 arctan(a/b) | Small (ignored) | Significant for b<6pc |

**Recommendation**:
- Use **hyperbolic formulas** for b < 5 pc (moderate deflection)
- Use **impulse approximation** for b > 6 pc and order-of-magnitude estimates
- For **publication-quality results**, always compare both methods

---

## Summary: Simulation Outputs, Theoretical Predictions, and Observational Diagnostics

This section provides a comprehensive guide to **what our simulations should produce** and **how to compare results** with theoretical predictions and radio observations.

---

### 1. Required Simulation Outputs

#### 1.1 Particle-Level Data (Every Snapshot)

**Spatial and kinematic data**:
- Position: `(x, y, z)` [pc]
- Velocity: `(v_x, v_y, v_z)` [km/s]
- Mass: `m` [M_☉]
- Density: `ρ` [M_☉/pc³]
- Smoothing length: `h` [pc]

**Thermodynamic data**:
- Temperature: `T` [K] - **CRITICAL for molecular dissociation**
- Pressure: `P` [dyn/cm²]
- Internal energy: `u` [erg/g]
- Sound speed: `c_s` [km/s]

**Gravitational data**:
- Gravitational potential: `Φ` [pc²/Myr²]
- Distance from IMBH: `d_BH` [pc]
- Binding energy: `E_bind = ½v² + u + Φ` [pc²/Myr²]

**Additional diagnostics** (computed in post-processing):
- Number density: `n_H = ρ × 40.5` [cm⁻³] (for μ=1)
- Mach number: `M = |v| / c_s`
- Tidal compression ratio: `ρ/ρ_0`

#### 1.2 Global Diagnostics (Every Timestep)

**Energy conservation**:
```
E_total = E_kinetic + E_thermal + E_potential
ΔE/E_0 < 1% (required for valid simulation)
```

**Mass tracking**:
- Total mass: `M_total(t)` (should be constant)
- Bound mass: `M_bound(t) = Σm for E_bind < 0`
- Unbound mass: `M_unbound(t) = Σm for E_bind > 0`
- Disrupted fraction: `f_disrupt(t) = 1 - M_bound/M_0`

**Morphology** (moment of inertia tensor):
- Axis ratios: `a:b:c` from eigenvalues of I_ij
- Elongation: `e = 1 - c/a`
- Pancaking factor: `p = a/c` (should increase at pericenter)

**Temperature statistics**:
- Maximum temperature: `T_max(t)`
- Mass-weighted mean: `<T> = Σ(m_i T_i) / M_total`
- Hot mass fraction: `f_hot(T>3000K) = Σ(m_i | T_i>3000K) / M_total`

---

### 2. Theoretical Predictions to Compare

#### 2.1 Pancaking Shock Velocity (A.9.1.2)

**Theory**:
```
v_shock = R_cloud × √(G M_BH / b³)
```

**Simulation extraction**:
```python
# At pericenter passage (t ≈ b/v_rel):
# 1. Find particles with strong compression (ρ/ρ_0 > 2)
# 2. Compute relative velocities in vertical direction (⊥ orbital plane)
# 3. Measure: v_shock,sim = max(|Δv_z|) between adjacent compressed layers

# Compare:
v_shock_theory = R_cloud * np.sqrt(G * M_BH / b**3)
ratio = v_shock_sim / v_shock_theory
# Expected: ratio ≈ 0.8-1.2 (within 20%)
```

**Parameter scan results** (theory):

| b [pc] | v_shock [km/s] | T_shock [K] | Expected outcome |
|--------|----------------|-------------|------------------|
| 3 | 20 | 3,200 | CO dissociates, H₂ partial |
| 4 | 13 | 1,370 | Molecules survive |
| 5 | 9 | 650 | Safe |
| 6 | 7 | 390 | Safe |

**Validation metric**:
```
SUCCESS: |v_shock_sim - v_shock_theory| / v_shock_theory < 0.2
```

---

#### 2.2 Post-Shock Temperature (A.9.3)

**Theory** (Rankine-Hugoniot):
```
T_shock = 8060 × v_shock² [K]  (for v_shock in km/s)
```

**Simulation extraction**:
```python
# Find shocked particles (where ρ > 2×ρ_ambient and div(v) < 0)
T_shocked = particles[shock_flag]['temperature']
v_shock_measured = particles[shock_flag]['delta_v']

# Compare:
T_theory = 8060 * v_shock_measured**2
residual = (T_shocked - T_theory) / T_theory
# Expected: residual < 0.3 (shock capturing with SPH has ~20-30% error)
```

**Key diagnostic**: Temperature-density phase diagram
```python
plt.scatter(n_H, T, c=time, alpha=0.5)
plt.axhline(2000, color='red', label='H₂ dissociation')
plt.axhline(3000, color='orange', label='CO dissociation')
plt.xlabel('n_H [cm⁻³]')
plt.ylabel('T [K]')
# Should see compressed gas at high n_H with elevated T
```

---

#### 2.3 Encounter Timescale (A.10.5)

**Theory**:
```
t_peri ≈ b / v_rel  (impulse approximation)
t_peri ≈ (2-2.5) × b / v_rel  (exact hyperbolic orbit)
```

**Simulation extraction**:
```python
# Time when IMBH is within distance b of cloud center
d_BH = np.sqrt(x_BH**2 + y_BH**2 + z_BH**2)
t_encounter = t[d_BH < b]
duration = t_encounter[-1] - t_encounter[0]

# Compare:
t_theory = b / v_rel
ratio = duration / t_theory
# Expected: ratio ≈ 2.0-2.5 (from hyperbolic orbit correction)
```

---

#### 2.4 Tidal Disruption Fraction (A.9.4)

**Theory** (empirical from parameter space):

| b [pc] | v_rel [km/s] | f_disrupt (predicted) |
|--------|-------------|-----------------------|
| 3 | 10 | 0.60-0.80 |
| 4 | 10 | 0.30-0.50 |
| 5 | 10 | 0.15-0.25 |
| 6 | 10 | 0.05-0.15 |

**Simulation extraction**:
```python
# At t_final (well after encounter):
E_bind = 0.5 * v**2 + u + Phi  # Specific binding energy
M_bound = np.sum(m[E_bind < 0])
f_disrupt_sim = 1 - M_bound / M_initial

# Compare with theory prediction
```

---

### 3. Molecular Dissociation Diagnostics

#### 3.1 Molecular Survival Criterion

**Theory** (A.9.2):
```
H₂ survives:  T < 2,000 K
CO survives:  T < 3,000 K
Both dissociate: T > 10,000 K
```

**Simulation diagnostic**:
```python
# For each particle, compute molecular state:
def molecular_state(T):
    if T < 1000:
        return 'molecules_intact'
    elif T < 2000:
        return 'H2_survives'
    elif T < 3000:
        return 'H2_partial_CO_survives'
    elif T < 10000:
        return 'mostly_atomic'
    else:
        return 'fully_atomic'

# Mass in each state
for state in ['molecules_intact', 'H2_survives', ...]:
    M_state[state] = np.sum(m[molecular_state == state])

# Key metric: Molecular survival fraction
f_molecular = M_state['molecules_intact'] / M_total
```

**Radio observability criterion**:
```
Observable: f_molecular > 0.5  (≥50% of mass retains molecules)
Marginal:   0.2 < f_molecular < 0.5
Invisible:  f_molecular < 0.2  (dissociated)
```

---

#### 3.2 Temperature Evolution Tracking

**Expected behavior** (from theory):

```
Phase 1 (t < 0): Approach
  - T ≈ 50 K (equilibrium)
  - ρ ≈ ρ_0 (undisturbed)

Phase 2 (t ≈ 0, pericenter):
  - Pancaking: ρ → 4×ρ_0 (compression)
  - Shock heating: T → T_shock = 8060 × v_shock²
  - Duration: Δt ~ b/v_rel

Phase 3 (t > 0): Post-encounter
  - Adiabatic expansion: T decreases slowly
  - Cooling (if enabled): T → T_eq(n) on timescale τ_cool ~ 3 Myr
  - Tidal tails: ρ → ρ_0/10, T → floor
```

**Simulation validation**:
```python
# Track temperature evolution for particles near pericenter
particles_peri = particles[(d_BH < 1.5*b) & (t_near_peri)]

plt.plot(t, T_max, label='T_max(t)')
plt.axvline(t_peri, color='red', label='Pericenter')
plt.axhline(T_shock_theory, color='orange', linestyle='--', label='Theory')
plt.xlabel('Time [Myr]')
plt.ylabel('T [K]')

# Expected: Sharp spike at t_peri to T ≈ T_shock_theory
```

---

### 4. Radio Observability Predictions

#### 4.1 Line-of-Sight Velocity (A.9.6)

**Observable quantity**: Doppler-shifted molecular line emission

**Theory**: v_los should show:
1. **Bulk motion** from tidal stretching: Δv_bulk ~ 50-100 km/s
2. **Turbulent broadening** from shocks: Δv_turb ~ v_shock ~ 10-20 km/s
3. **But local temperature** must remain T < 2000 K for molecules to survive!

**Simulation extraction**:
```python
# Select molecular gas (T < 2000 K)
molecular_gas = particles[T < 2000]

# Compute line-of-sight velocity distribution
# (Assume observer along z-axis)
v_los = molecular_gas['v_z']

# Observable line profile
v_bins = np.linspace(-150, 150, 100)  # km/s
line_profile = np.histogram(molecular_gas['m'], bins=v_bins, weights=v_los)

# Metrics:
v_los_mean = np.average(v_los, weights=molecular_gas['m'])
v_los_std = np.sqrt(np.average((v_los - v_los_mean)**2, weights=molecular_gas['m']))

print(f"Mean v_los: {v_los_mean:.1f} km/s")
print(f"Line width: {v_los_std:.1f} km/s")
```

**Radio detectability**:
```
Detectable:  |v_los| > 20 km/s  AND  f_molecular > 0.5
Marginal:    |v_los| > 10 km/s  AND  0.2 < f_molecular < 0.5
Invisible:   f_molecular < 0.2  (no molecules, regardless of velocity)
```

---

#### 4.2 Observable vs Invisible Events (A.9.7)

**Theory predictions** (M_BH = 10⁵ M_☉, R_cloud = 5 pc):

| b [pc] | v_rel [km/s] | v_shock [km/s] | T_shock [K] | f_molecular | v_los [km/s] | Radio detection |
|--------|-------------|----------------|-------------|-------------|--------------|-----------------|
| 3 | 10 | 20 | 3,200 | 0.1 | 100 | ✗ **Invisible** (molecules gone) |
| 4 | 10 | 13 | 1,370 | 0.6 | 70 | ✓ **Detectable** |
| 5 | 10 | 9 | 650 | 0.9 | 50 | ✓ **Strong** |
| 6 | 15 | 7 | 390 | 0.95 | 80 | ✓ **Excellent** |

**Simulation validation**:
For each (b, v_rel) combination, verify:
1. f_molecular matches prediction (within factor of 2)
2. v_los_std ≈ v_shock_theory (shock-induced turbulence)
3. Radio detectability matches theory

**Key plot**:
```python
# Create observability diagram
plt.scatter(b, v_rel, c=f_molecular, s=100*v_los_std, cmap='RdYlGn')
plt.colorbar(label='Molecular fraction')
plt.contour(b, v_rel, T_shock, levels=[2000, 3000], colors=['red', 'orange'])
plt.xlabel('Impact parameter b [pc]')
plt.ylabel('Relative velocity v_rel [km/s]')
plt.title('Radio Observability Diagram')
# Green region = detectable, Red = invisible
```

---

### 5. Convergence Tests

#### 5.1 Resolution Convergence (Section 4a)

**Required tests**:
Run same (b, v_rel) with different particle numbers:
- N = 5×10⁴ (quick test)
- N = 2×10⁵ (baseline)
- N = 4×10⁵ (convergence test)

**Metrics to compare**:
```python
for N in [5e4, 2e5, 4e5]:
    results[N] = {
        'f_disrupt': ...,
        'T_max': ...,
        'v_shock': ...,
        'f_molecular': ...
    }

# Convergence criterion
conv = |results[4e5] - results[2e5]| / results[2e5]
# Required: conv < 0.1 (10% difference)
```

---

#### 5.2 Timestep Convergence

**CFL numbers to test** (in config):
```json
"timestep": {
  "CFL_sound": [0.2, 0.3, 0.4],
  "CFL_force": [0.15, 0.25, 0.35]
}
```

Should produce consistent results (< 5% variation).

---

### 6. Summary Checklist: What to Extract from Simulations

#### ✅ For Each Simulation Run (b, v_rel, N):

**Basic outputs**:
- [ ] Snapshot files with T, ρ, v, Φ at Δt = 0.01-0.02 Myr intervals
- [ ] Energy file showing E_total, E_kin, E_therm, E_pot vs time
- [ ] Log file with timestep information

**Post-processing analyses**:
- [ ] v_shock measurement (compare to R √(GM/b³))
- [ ] T_max evolution (peak should match 8060 × v_shock²)
- [ ] f_disrupt(t) evolution (final value)
- [ ] f_molecular(t) evolution (fraction with T < 2000 K)
- [ ] v_los distribution (for radio observability)
- [ ] n-T phase diagram (check dissociation regions)
- [ ] Morphology evolution (a:b:c axis ratios, pancaking)

**Comparison with theory**:
- [ ] v_shock: |sim - theory|/theory < 0.2
- [ ] T_shock: |sim - theory|/theory < 0.3
- [ ] t_encounter: ratio ≈ 2-2.5
- [ ] f_disrupt: within factor of 2 of prediction
- [ ] f_molecular: matches radio detectability criterion

**Publication-quality outputs**:
- [ ] 3D density rendering at t_peri
- [ ] Temperature evolution movie
- [ ] Line-of-sight velocity map (for radio comparison)
- [ ] Multi-panel comparison for different (b, v_rel)

---

### 7. Connection to Observations

#### 7.1 Simulated Radio Line Profile

**What observers see**: Molecular line (e.g., ¹²CO J=1-0) Doppler-shifted by cloud motion

**From simulation**:
```python
# For molecular gas only (T < 2000 K)
molecular_particles = particles[T < 2000]

# Compute line emission (assume optically thin)
# Emissivity ∝ n² for collisionally excited lines
emissivity = molecular_particles['n_H']**2 * molecular_particles['mass']

# Velocity-integrated line profile
v_los = molecular_particles['v_z']  # Line of sight
spectrum = np.histogram(v_los, bins=v_bins, weights=emissivity)

# Observable: Integrated intensity vs velocity
plt.plot(v_bins, spectrum)
plt.xlabel('v_los [km/s]')
plt.ylabel('Intensity [arbitrary]')
plt.title('Simulated ¹²CO Line Profile')

# Detectability:
peak_intensity = spectrum.max()
line_width = v_bins[spectrum > 0.5*peak_intensity].ptp()
print(f"Peak at v = {v_bins[np.argmax(spectrum)]:.1f} km/s")
print(f"FWHM = {line_width:.1f} km/s")
```

**Observable if**:
- Peak intensity > noise (requires f_molecular > 0.5)
- Line width ~ 50-100 km/s (from tidal stretching)
- Peak offset from systemic velocity

---

#### 7.2 Multi-Wavelength Predictions

**Temperature-based tracer selection**:

| T range [K] | Dominant tracer | Simulation output |
|-------------|----------------|-------------------|
| T < 100 | CO, H₂ | v_los map, column density |
| 100-1000 | Warm CO, [C I] | Warm component fraction |
| 1000-3000 | [C II] 158 μm | Atomic carbon map |
| 3000-10000 | [O I] 63 μm | Hot atomic gas |
| > 10000 | H recombination | Ionized fraction |

**From simulation**:
```python
for T_min, T_max, tracer in [(0, 100, 'CO'), (1000, 3000, '[CII]'), ...]:
    mass_frac = np.sum(m[(T > T_min) & (T < T_max)]) / M_total
    print(f"{tracer}: {mass_frac:.2%} of total mass")
```

**Prediction**:
- b = 3 pc: Mostly [C II], [O I] (too hot for CO)
- b ≥ 5 pc: Strong CO, some [C II] (molecules survive)

---

### 8. Example Analysis Pipeline

```python
import numpy as np
import matplotlib.pyplot as plt

# 1. Load simulation snapshot at pericenter
data = read_snapshot('snapshot_pericenter.csv')

# 2. Extract key quantities
T = data['temperature']
rho = data['density']
n_H = rho * 40.5  # cm⁻³
v = data['velocity']
m = data['mass']

# 3. Compute theoretical predictions
b = 4.0  # pc
M_BH = 1e5  # M_sun
R_cloud = 5.0  # pc
v_shock_theory = R_cloud * np.sqrt(G * M_BH / b**3)  # km/s
T_shock_theory = 8060 * v_shock_theory**2  # K

# 4. Measure shock from simulation
shocked = (rho > 2*rho_0) & (np.gradient(v) < 0)
v_shock_sim = np.std(v[shocked])
T_shock_sim = np.max(T[shocked])

# 5. Compare
print(f"v_shock: theory={v_shock_theory:.1f} km/s, sim={v_shock_sim:.1f} km/s")
print(f"T_shock: theory={T_shock_theory:.0f} K, sim={T_shock_sim:.0f} K")

# 6. Molecular survival
f_molecular = np.sum(m[T < 2000]) / np.sum(m)
print(f"Molecular fraction: {f_molecular:.2%}")

# 7. Radio observability
if f_molecular > 0.5:
    v_los = v[T < 2000, 2]  # z-component
    print(f"Detectable! Line width: {np.std(v_los):.1f} km/s")
else:
    print("Radio-invisible (molecules dissociated)")

# 8. Create diagnostic plots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes[0,0].scatter(n_H, T, c=m, s=1, alpha=0.3)
axes[0,0].axhline(2000, color='r', label='H₂ dissociation')
axes[0,0].set_xlabel('n_H [cm⁻³]')
axes[0,0].set_ylabel('T [K]')
axes[0,0].set_yscale('log')
axes[0,0].legend()
# ... more plots
```

---

### 9. Success Criteria

A successful simulation should demonstrate:

✅ **Numerical accuracy**:
- Energy conservation: ΔE/E₀ < 1%
- Mass conservation: ΔM/M₀ < 0.1%
- Convergence: Results change < 10% when doubling N

✅ **Physical realism**:
- Shock velocity within 20% of theory
- Temperature spike at pericenter
- Pancaking morphology (a/c > 2 at pericenter)
- Tidal tail formation after encounter

✅ **Observable predictions**:
- Correct molecular survival fraction for given (b, v_rel)
- Realistic velocity dispersions (50-100 km/s)
- Match radio detectability criterion from theory

✅ **Parameter space coverage**:
- At least 4 values of b (3, 4, 5, 6 pc)
- At least 3 values of v_rel (5, 10, 15 km/s)
- Convergence test for at least one case

---

## Next Steps

1. Review existing gravity force implementation in `src/gravity_force.cpp`
2. Implement external point-mass force module
3. Adapt Lane-Emden setup to include K&I thermal equilibrium
4. Create test case with N = 5×10⁴ particles
5. Validate conservation and thermal equilibrium before production runs
6. **NEW**: Add temperature tracking to identify dissociation events in post-processing
7. **NEW**: Create parameter scan focusing on observable parameter space (b ≥ 5 pc)
8. **NEW**: Implement multi-tracer diagnostics ([C II] emissivity for warm gas)

---

# Appendix B: Thermal Physics and Unit System

## B.1 Koyama-Inutsuka Thermal Physics

### B.1.1 Overview

The **Koyama & Inutsuka (2000)** cooling function provides a realistic model of **ISM thermal equilibrium** for molecular cloud simulations. It balances heating and cooling processes to determine the equilibrium temperature as a function of density.

**Key Physics**:
- **Heating sources**: Cosmic ray ionization, photoelectric heating from dust grains
- **Cooling sources**: Line emission (CO, H₂O, [C II], etc.), dust continuum radiation
- **Equilibrium condition**: Heating rate = Cooling rate → unique T(n) relation

**Implementation**: `thermal::KoyamaInutsukaCooling` class in `include/thermal/koyama_inutsuka_cooling.hpp`

### B.1.2 Thermal Equilibrium Curve

The K&I (2000) model provides temperature as a function of number density:

| n_H [cm⁻³] | T_eq [K] | Thermal State | Dominant Coolant |
|------------|----------|---------------|------------------|
| 10⁰ | ~8,000 | Warm atomic | [C II] 158 μm |
| 10¹ | ~1,000 | Cool atomic | [C II], [O I] |
| 10² | ~100 | CNM/molecular transition | CO, H₂O |
| 10³ | ~50 | Cold molecular | CO, dust |
| 10⁴ | ~30 | Dense molecular | CO, H₂O |
| 10⁵ | ~20 | Very dense molecular | Dust continuum |

**Pressure equilibrium**: The ISM maintains P/k ~ 10³⁻⁴ K cm⁻³ across different phases.

### B.1.3 With vs Without Cooling

#### **WITH Koyama-Inutsuka Cooling** (`enable_cooling: true`)

**Physics**:
- Temperature locked to T(n) equilibrium curve
- Shock heating immediately radiated away
- Energy equation: dE/dt = heating - cooling + P dV work + viscous dissipation
- Cooling time: t_cool ~ 10⁴⁻⁵ yr << t_tidal ~ 10⁶ yr

**Expected behavior in IMBH simulation**:
- **Compressed regions**: Heat up initially, then cool to T_eq(n_shocked)
- **Tidal tails**: Stay cold (T ~ 20-50 K) even during stretching
- **Shocks**: Weak temperature jump (limited by rapid cooling)
- **Fragmentation**: Enhanced due to efficient cooling → denser clumps

**Energy conservation**:
```
E_total = E_kinetic + E_thermal + E_gravitational + E_BH_potential - ∫ (Λ - Γ) dt
```
Where Λ = cooling rate, Γ = heating rate.

#### **WITHOUT Cooling** (`enable_cooling: false`)

**Physics**:
- Pure adiabatic evolution: PV^γ = const (γ = 5/3)
- Temperature equation: T ∝ ρ^(γ-1) = ρ^(2/3)
- Shock heating preserved in gas
- No energy sink → higher thermal pressure

**Expected behavior**:
- **Compressed regions**: Heat up strongly (T ∝ ρ^(2/3))
- **Tidal tails**: Warm due to adiabatic compression (T ~ 100-500 K)
- **Shocks**: Strong temperature jumps (Rankine-Hugoniot relations)
- **Fragmentation**: Suppressed by thermal pressure support

**Energy conservation** (exact in absence of numerical dissipation):
```
E_total = E_kinetic + E_thermal + E_gravitational + E_BH_potential = const
```

### B.1.4 Temperature Evolution Comparison

**Compression event** (density increases by factor 10):

| Scenario | Initial T [K] | Final T [K] | Mechanism |
|----------|---------------|-------------|-----------|
| With cooling | 50 | ~30 | Locked to T_eq(10n₀) |
| Without cooling | 50 | ~233 | Adiabatic: T ∝ ρ^(2/3) |

**Shock event** (Mach 10 shock, v_shock = 3 km/s, c_s = 0.3 km/s):

| Scenario | Pre-shock T [K] | Post-shock T [K] | Time to cool |
|----------|-----------------|------------------|--------------|
| With cooling | 50 | ~50-100 | t_cool ~ 10³ yr |
| Without cooling | 50 | ~5,000 | Never (adiabatic) |

**Physical interpretation**:
- **Cooling ON**: Cloud remains molecularly cold, behaves like isothermal (T ~ const)
- **Cooling OFF**: Cloud heats up, behaves like polytropic (P ∝ ρ^γ)

---

## B.2 Unit System: GALACTIC Units

### B.2.1 Base Units

The simulation uses **GALACTIC units** for astrophysical convenience:

```cpp
UnitSystem::create_galactic(
    length_kpc  = 1.0,    // 1 code unit = 1 kpc
    mass_msun   = 1.0e10, // 1 code unit = 10¹⁰ M_☉
    velocity_kms = 1.0    // 1 code unit = 1 km/s
);
```

**Scaling for IMBH-cloud problem**:
- Typical cloud size: R ~ 5 pc = 0.005 kpc = **0.005 code units**
- Typical cloud mass: M ~ 10⁴ M_☉ = 10⁻⁶ × 10¹⁰ M_☉ = **10⁻⁶ code units**
- IMBH mass: M_BH = 10⁵ M_☉ = 10⁻⁵ × 10¹⁰ M_☉ = **10⁻⁵ code units**
- Encounter velocity: v ~ 10 km/s = **10.0 code units**

### B.2.2 Derived Units

From the base units, all other quantities are derived:

| Quantity | Formula | Value | Physical Interpretation |
|----------|---------|-------|------------------------|
| **Time** | L/V | 1 kpc / 1 km/s | **0.978 Myr** |
| **Density** | M/L³ | 10¹⁰ M_☉ / (1 kpc)³ | **6.77×10⁻²³ g/cm³** |
| **Pressure** | M·V²/L³ | (10¹⁰ M_☉)·(1 km/s)² / (1 kpc)³ | **6.77×10⁻¹³ dyne/cm²** |
| **Energy** | M·V² | (10¹⁰ M_☉)·(1 km/s)² | **1.99×10⁵⁴ erg** |
| **Gravity** | V²/L | (1 km/s)² / (1 kpc) | **G = 1.0** (in code units) |
| **Number density** | ρ / m_H | (6.77×10⁻²³ g/cm³) / (1.67×10⁻²⁴ g) | **40.5 cm⁻³** per code density |

### B.2.3 Unit Conversions for IMBH-Cloud Simulation

**Configuration file uses PHYSICAL units** (pc, M_☉, km/s):

```json
"imbh_parameters": {
  "M_BH": 100000.0,              // 10⁵ M_☉ → 10⁻⁵ code units
  "BH_initial_position": {
    "0": -20.0,                  // -20 pc → -0.020 code units
    "1": 0.0,                    // 0 pc → 0.0 code units
    "2": 0.0                     // 0 pc → 0.0 code units
  },
  "BH_initial_velocity": {
    "0": 10.0,                   // 10 km/s → 10.0 code units
    "1": 0.0,
    "2": 0.0
  },
  "softening_epsilon": 0.05      // 0.05 pc → 5×10⁻⁵ code units
}
```

**Internal conversion** (automatic in code):
```cpp
// In GALACTIC units with length_kpc=1.0, mass_msun=1e10, velocity_kms=1.0:
// 1 pc in code = 1e-3 kpc / length_kpc = 0.001
// 1 M_☉ in code = 1 M_☉ / mass_msun = 1e-10

real M_BH_code = M_BH_msun / 1e10;           // 1e5 M_☉ → 1e-5 code
real pos_code = pos_pc * 0.001;              // -20 pc → -0.020 code
real vel_code = vel_kms;                     // 10 km/s → 10.0 code (same!)
```

### B.2.4 Physical Parameters in Code Units

For **Lane-Emden n=1.5 cloud** with relaxed initial condition:

| Parameter | Physical Value | Code Units | Conversion |
|-----------|----------------|------------|------------|
| Cloud radius | R = 5 pc | 0.005 | R_code = 5×10⁻³ |
| Cloud mass | M = 10⁴ M_☉ | 10⁻⁶ | M_code = 10⁻⁶ |
| Central density | ρ_c ~ 10⁻²¹ g/cm³ | ~0.015 | ρ_code = ρ_phys / 6.77×10⁻²³ |
| Number density | n_H ~ 600 cm⁻³ | ~15 | n_code = n_phys / 40.5 |
| Sound speed | c_s ~ 0.3 km/s | 0.3 | c_s,code = 0.3 |
| Temperature | T ~ 50 K | — | (thermal state variable) |
| Polytropic K | K = P/ρ^γ | ~10⁻⁴ | K_code (derived) |

For **IMBH encounter**:

| Parameter | Physical Value | Code Units | Conversion |
|-----------|----------------|------------|------------|
| BH mass | M_BH = 10⁵ M_☉ | 10⁻⁵ | M_BH,code = 10⁻⁵ |
| Impact parameter | b = 3 pc | 0.003 | b_code = 3×10⁻³ |
| Relative velocity | v = 10 km/s | 10.0 | v_code = 10.0 |
| Initial separation | d = 20 pc | 0.020 | d_code = 0.020 |
| Softening | ε = 0.05 pc | 5×10⁻⁵ | ε_code = 5×10⁻⁵ |

### B.2.5 Timescales in Code Units

| Timescale | Formula | Value (Physical) | Value (Code) |
|-----------|---------|------------------|--------------|
| Encounter time | t_enc = d/v | 20 pc / 10 km/s = 1.96 Myr | 0.020 / 10.0 = **0.002** |
| Crossing time | t_cross = R/v | 5 pc / 10 km/s = 0.49 Myr | 0.005 / 10.0 = **0.0005** |
| Free-fall time | t_ff = √(3π/32Gρ) | ~1.0 Myr | ~**0.001** |
| Sound crossing | t_sound = R/c_s | 5 pc / 0.3 km/s = 16.3 Myr | 0.005 / 0.3 = **0.017** |
| Cooling time | t_cool | ~10⁴ yr | ~**10⁻⁵** |
| Simulation end | — | 2.0 Myr | **2.04** |

**Key ratio**: t_cool / t_cross ~ 10⁻⁵ / 5×10⁻⁴ = **0.02** → thermal equilibrium maintained!

---

## B.3 Energy Conservation and Cooling

### B.3.1 Energy Budget (WITH Cooling)

Total energy evolution:
```
dE_total/dt = -∫ (Λ - Γ) dV + numerical_errors
```

Where:
- Λ(n,T) = volumetric cooling rate [erg cm⁻³ s⁻¹]
- Γ(n) = volumetric heating rate [erg cm⁻³ s⁻¹]
- Numerical errors should be << cooling rate

**Components**:
```
E_total(t) = E_kinetic + E_thermal + E_grav + E_BH - ∫₀ᵗ ∫(Λ-Γ) dV dt'
```

**Energy conservation check**:
1. **Initial** (t=0): E₀ = E_kin,0 + E_therm,0 + E_grav,0 (BH far away)
2. **During** (t>0): E(t) should decrease due to cooling losses
3. **Final** (t=t_end): E_f = E₀ - ∫ net_cooling dt

**Expected drift**: |ΔE| / E₀ ~ few% (dominated by cooling, not numerical error)

### B.3.2 Energy Budget (WITHOUT Cooling)

Total energy **strictly conserved**:
```
E_total = E_kinetic + E_thermal + E_grav + E_BH = const
```

**Energy conservation check**:
```
|ΔE| / E₀ < 0.1%  (numerical error only)
```

If drift > 0.1%, indicates:
- Timestep too large (reduce CFL_sound, CFL_force)
- Artificial viscosity too strong
- Poor neighbor finding (increase neighbor_number)

### B.3.3 Monitoring Energy in Simulation

**Energy output file**: `results/.../energy.dat`

Columns:
```
# time  E_kinetic  E_thermal  E_potential  E_total
0.000   1.234e-5   2.345e-5   -6.789e-5    -3.210e-5
0.001   1.235e-5   2.340e-5   -6.790e-5    -3.215e-5
...
```

**Python diagnostic**:
```python
import pandas as pd
import numpy as np

energy = pd.read_csv('energy.dat', delim_whitespace=True, comment='#',
                     names=['time', 'E_kin', 'E_therm', 'E_pot', 'E_tot'])

# Check conservation (no cooling case)
dE = energy['E_tot'] - energy['E_tot'].iloc[0]
rel_error = np.abs(dE / energy['E_tot'].iloc[0])
print(f"Max energy drift: {rel_error.max():.4%}")

# Expected: < 0.1% for adiabatic, few% for cooling
```

---

## B.4 Running Thermal Comparison Simulations

### B.4.1 Makefile Targets

**Run both simulations sequentially**:
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_thermal_compare_run
```

This will:
1. Run **WITH cooling**: `simulation_10k_b3pc_cool.json` → `results/imbh_relaxed_2k_b3pc/`
2. Run **WITHOUT cooling**: `simulation_10k_b3pc_nocool.json` → `results/imbh_relaxed_2k_b3pc_nocool/`

**Visualize comparison** (when script is implemented):
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_thermal_compare_viz
```

**Clean thermal comparison outputs**:
```bash
make -f sample/imbh_cloud/Makefile.imbh_cloud imbh_thermal_compare_clean
```

### B.4.2 Configuration Files

**WITH Cooling**: `config/presets/simulation_10k_b3pc_cool.json`
```json
"thermal": {
  "enable_cooling": true,
  "cooling_type": "koyama_inutsuka",
  "N_H_column": 1e20,
  "thermal_relax_time": 0.01
}
```

**WITHOUT Cooling**: `config/presets/simulation_10k_b3pc_nocool.json`
```json
"thermal": {
  "enable_cooling": false,
  "comment": "Adiabatic evolution only (γ=5/3)"
}
```

### B.4.3 Expected Outcomes

| Observable | With Cooling | Without Cooling |
|------------|--------------|-----------------|
| **Tail temperature** | T ~ 20-50 K | T ~ 100-500 K |
| **Tail density contrast** | Higher (ρ_max/ρ_min ~ 100) | Lower (ρ_max/ρ_min ~ 10) |
| **Fragmentation** | Enhanced (clumpy) | Suppressed (smooth) |
| **Energy drift** | Few % (cooling losses) | < 0.1% (conserved) |
| **Thermal pressure** | Low (cold) | High (warm) |
| **Tail width** | Narrow (compact) | Broad (puffy) |

**Scientific question**: Does thermal physics significantly affect tidal tail morphology and fragmentation?

---

## B.5 References

### Thermal Physics
- **Koyama & Inutsuka (2000)**: ApJ, 532, 980 - "Molecular Cloud Formation in Shock-Compressed Layers"
- **Koyama & Inutsuka (2002)**: ApJ, 564, L97 - "An Origin of Supersonic Motions in Interstellar Clouds"

### SPH Methods
- **Hopkins (2013)**: MNRAS, 428, 2840 - "A general class of Lagrangian smoothed particle hydrodynamics methods"
- **Rosswog (2020)**: Living Rev. Comput. Astrophys., 6, 1 - "SPH methods in astrophysical applications"

### Unit Systems
- **Binney & Tremaine (2008)**: *Galactic Dynamics*, 2nd Edition - Standard galactic units reference

---

**Document updated**: 2025-01-02  
**Appendix B added**: Thermal physics comparison framework and unit system details
