# Hydrostatic Self-Gravity Test for IMBH-Cloud Simulations

## Overview

Before running IMBH tidal disruption simulations, we must verify that the self-gravity implementation is working correctly. This test ensures that a relaxed Lane-Emden polytrope maintains hydrostatic equilibrium over timescales much longer than the tidal disruption and scattering dynamics.

## Purpose

**Why this test is critical:**
1. **Verify gravity-pressure balance:** Self-gravity must correctly balance pressure forces
2. **Detect spurious forces:** Identify any artificial forces from discretization errors
3. **Test numerical stability:** Ensure long-term integration preserves equilibrium
4. **Validate before science runs:** Catch implementation bugs before expensive production runs

## Physics

### Lane-Emden Polytrope (n=1.5, γ=5/3)

Hydrostatic equilibrium condition:
```
dP/dr = -ρ(r) ∇Φ(r)
```

For polytropic equation of state:
```
P = K ρ^γ    where γ = 5/3 (ideal gas)
```

Analytical Lane-Emden solution:
```
ρ(r) = ρ_c θ(ξ)^n    where ξ = r/α, α = R/3.654, n = 1.5
```

### Timescales

All times in code units (where R = 1.0, M = 1.0):

| Timescale | Value | Definition | Physical Meaning |
|-----------|-------|------------|------------------|
| **Dynamical time** | t_dyn = 1.0 | sqrt(R³/GM) | Free-fall collapse time (code unit) |
| **Crossing time** | t_cross ~ 1.0 | R/c_s | Sound crossing time |
| **Tidal time (IMBH)** | t_tidal ~ 3.2 | sqrt(R³/GM_BH) | For M_BH/M_cloud = 10 |
| **Scattering time** | t_scatter ~ 0.3-3 | b/v_rel | Depends on impact parameter & velocity |
| **Test duration** | t_test = 10.0 | 10 × t_cross | >> all physical timescales |

**Test requirement:** t_test >> t_tidal >> t_scatter ensures we test stability over timescales relevant to IMBH interactions.

## Test Configuration

### 2k Particle Test (Quick validation)

```json
{
  "N": 22,                          // ~2,130 particles
  "SPHType": "disph",               // Conservative formulation
  "kernel": "wendland",             // C² continuous kernel
  "useGravity": true,               // ★ Self-gravity enabled
  "neighborNumber": 50,
  
  "endTime": 10.0,                  // 10 crossing times
  "outputTime": 0.2,                // 50 snapshots
  
  "resumeFromSnapshot": "<path>"    // Relaxed initial conditions
}
```

**Runtime:** ~5-10 minutes  
**Use case:** Quick verification, development testing

### 200k Particle Test (Production validation)

```json
{
  "N": 100,                         // 200,000 particles
  "SPHType": "disph",
  "kernel": "wendland",
  "useGravity": true,
  
  "endTime": 10.0,
  "outputTime": 0.2,
  
  "resumeFromSnapshot": "<path>"
}
```

**Runtime:** ~8-16 hours  
**Use case:** Production quality verification before expensive IMBH runs

## Expected Behavior

### ✅ Pass Criteria

1. **Density profile stability**
   - RMS deviation from analytic Lane-Emden: < 5%
   - No systematic drift or oscillations
   - Profile shape unchanged over 10 t_cross

2. **Velocity suppression**
   - Maximum velocity: < 1% of sound speed
   - Median velocity: < 0.1% of sound speed
   - No coherent bulk flows

3. **Energy conservation**
   - Total energy drift: < 1% over entire test
   - Smooth evolution (no jumps or spikes)
   - Virial equilibrium maintained

4. **Structural stability**
   - No expansion or contraction
   - Spherical symmetry preserved
   - No artificial fragmentation

### ❌ Failure Modes

| Symptom | Possible Cause | Action |
|---------|----------------|--------|
| Density profile drifts | Gravity-pressure imbalance | Check gravity force calculation |
| Large velocities (> 5% c_s) | Spurious forces | Check kernel gradient, tree opening angle |
| Energy drift > 5% | Numerical instability | Reduce timestep, check integrator |
| Systematic expansion | Gravity too weak | Check G constant, softening length |
| Systematic collapse | Gravity too strong | Check force normalization |
| Oscillations | Poor relaxation | Increase relaxation steps |

## Usage

### Quick Test (2k particles)

```bash
# Step 1: Run relaxation (if not done)
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_relax_2k

# Step 2: Run hydrostatic test
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_2k

# Step 3: View results
open simulations/astrophysics/imbh_cloud/results/lane_emden_2k_hydrostatic/hydrostatic_animation.gif
```

### Production Test (200k particles)

```bash
# Step 1: Run high-resolution relaxation (4-8 hours)
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_relax_200k

# Step 2: Run production hydrostatic test (8-16 hours)
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_200k

# Step 3: Verify results
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_200k_viz
```

### Visualization Only

```bash
# Re-generate plots from existing data
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_hydrostatic_2k_viz
```

## Output Files

### Results Directory Structure

```
results/lane_emden_2k_hydrostatic/
├── snapshot_0000.csv            # Initial state (t=0)
├── snapshot_0001.csv            # t = 0.2
├── ...
├── snapshot_0049.csv            # Final state (t=10.0)
├── energy.dat                   # Time-series: [time, E_tot, E_kin, E_thermal, E_pot]
├── run.log                      # Full simulation output
├── hydrostatic_animation.gif    # 6-panel evolution movie
└── hydrostatic_summary.png      # Final state analysis
```

### Visualization Panels

**hydrostatic_animation.gif** (6 panels):
1. **3D Density Distribution:** XY projection showing particle distribution
2. **XY Slice:** Equatorial cut (|z| < 0.1) with density colormap
3. **Radial Density Profile:** SPH binned vs Lane-Emden analytic (RED LINE)
4. **Velocity Histogram:** Distribution of |v|/c_s (should be << 1%)
5. **Energy Evolution:** Fractional drift from initial energy
6. **Quality Metrics:** Time evolution of density error and max velocity

**hydrostatic_summary.png** (6 panels):
1. Final density profile comparison (SPH vs analytic)
2. Density error vs radius
3. Velocity distribution vs radius
4. RMS density error time series
5. Velocity evolution (max & median)
6. Test summary statistics

## Interpretation Guide

### Reading the Animation

**Panel 3 - Radial Density Profile** (most important):
- **Red line:** Analytical Lane-Emden solution ρ(r) = ρ_c θ^1.5
- **Blue points:** SPH binned density
- **Pass:** Blue points track red line throughout animation
- **Fail:** Blue points drift away from red line

**Panel 4 - Velocity Histogram**:
- **X-axis:** Velocity normalized by sound speed (|v|/c_s)
- **Pass:** Peak near zero, max < 0.01
- **Fail:** Broad distribution, max > 0.05

**Panel 5 - Energy Evolution**:
- **Y-axis:** (E - E₀)/|E₀| in percent
- **Pass:** Stays within ±1% (red dashed lines)
- **Fail:** Monotonic drift or oscillations > 5%

**Panel 6 - Quality Metrics**:
- **Blue line:** RMS density error (%)
- **Red line:** Max velocity (% of c_s)
- **Pass:** Both lines stay below 5% threshold
- **Fail:** Trending upward or spikes

### Summary Statistics

Located in Panel 6 of `hydrostatic_summary.png`:

```
HYDROSTATIC TEST SUMMARY
════════════════════════════

Test Duration:
  • Total time: 10.00 code units
  • Crossing times: 10.2 t_cross

Final State Quality:
  • RMS density error: 2.34%      ✓ PASS (< 5%)
  • Max velocity: 0.45% c_s       ✓ PASS (< 1%)
  • Median velocity: 0.012% c_s   ✓ PASS (< 0.1%)

Pass Criteria:
  ✓ Density error < 5%
  ✓ Max velocity < 1% c_s
  ✓ Median vel < 0.1% c_s
```

## Physical Validation

### Why DISPH Method?

**DISPH (Density-Independent SPH)** advantages for this test:
1. **Conservative formulation:** Better energy conservation
2. **Pressure-entropy formulation:** Reduces spurious density oscillations
3. **No artificial viscosity required:** Cleaner test of gravity-pressure balance
4. **Hydrostatic equilibrium:** Naturally maintains equilibrium better than SSPH

### Connection to IMBH Physics

**Tidal disruption timescale:**
```
t_tidal = sqrt(R³ / (G M_BH))
```

For M_BH = 10⁵ M_☉, M_cloud = 10⁴ M_☉:
```
M_BH / M_cloud = 10
t_tidal ≈ sqrt(10) × t_dyn ≈ 3.2 code units
```

**Test requirement:**
```
t_test = 10.0 >> t_tidal ≈ 3.2
```

If the sphere maintains equilibrium for t = 10.0, then self-gravity is accurate enough for IMBH simulations where tidal effects occur on timescale ~ 3 code units.

## Troubleshooting

### "ERROR: Relaxed initial conditions not found"

**Cause:** Haven't run relaxation yet.

**Solution:**
```bash
make -f simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud imbh_relax_2k
```

### Test fails: Density error > 10%

**Diagnosis:**
1. Check relaxation quality: `imbh_relax_2k_viz`
2. Increase relaxation steps if velocities not converged
3. Check gravity implementation (tree opening angle θ)

**Common causes:**
- Insufficient relaxation (particles not at rest)
- Tree opening angle too large (inaccurate forces)
- Softening length too large

### Test fails: Velocities > 5% c_s

**Diagnosis:**
1. Check for spurious forces from discretization
2. Verify kernel gradient normalization
3. Check timestep stability (CFL conditions)

**Actions:**
- Reduce CFL numbers (cflSound, cflForce)
- Try different kernel (cubic spline vs Wendland)
- Increase neighbor number (50 → 100)

### Energy drift > 5%

**Diagnosis:**
- Numerical integration error accumulation
- Timestep too large
- Tree force accuracy insufficient

**Actions:**
```json
{
  "cflSound": 0.2,        // Reduce from 0.3
  "cflForce": 0.1,        // Reduce from 0.125
  "maxTreeLevel": 25,     // Increase from 20
  "treeOpeningAngle": 0.4 // Reduce from 0.5
}
```

## Success Confirmation

### Before Proceeding to IMBH Simulations

**Required checks:**
- [x] Relaxation completed successfully
- [x] Hydrostatic test run for t = 10 t_cross
- [x] Density RMS error < 5%
- [x] Max velocity < 1% c_s
- [x] Energy drift < 1%
- [x] Visual inspection of animation shows stable profile

**If all checks pass:**
```
✅ Self-gravity implementation verified!
✅ Ready for IMBH tidal disruption simulations
```

**If any check fails:**
```
❌ Self-gravity issues detected
⚠️  Do NOT proceed to IMBH runs
⚠️  Debug gravity implementation first
```

## References

**Lane-Emden Solution:**
- Chandrasekhar, S. (1939), "An Introduction to the Study of Stellar Structure"

**Hydrostatic Equilibrium in SPH:**
- Springel, V. (2010), "Smoothed Particle Hydrodynamics in Astrophysics", ARAA 48, 391

**DISPH Method:**
- Hopkins, P. F. (2013), "A New Class of Accurate, Mesh-Free Hydrodynamic Simulation Methods", MNRAS 428, 2840

**IMBH Tidal Disruption:**
- Stone, N., Sari, R., & Loeb, A. (2013), "Consequences of strong compression in tidal disruption events", MNRAS 435, 1809

---

**Last updated:** 2025-12-01  
**Contact:** Guo  
**Repository:** sph-code-simulation
