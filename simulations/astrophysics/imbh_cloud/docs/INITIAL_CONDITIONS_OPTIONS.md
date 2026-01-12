# Initial Conditions for IMBH-Cloud Simulation (Oka et al. 2017)

## The Fundamental Problem

You're trying to simulate the HVCC CO-0.40-0.22 (Oka et al. 2017) using SPH. The core issue is:

**The observed cloud (M = 40 M☉, R ~ 0.2 pc, T = 60 K) is NOT in equilibrium.**

The virial ratio M_vir / M_actual ~ 100 proves the cloud is **tidally dominated**, not pressure-supported.

### Why Bonnor-Ebert + K&I 2000 Fails

| What you want | What BE physics gives | Problem |
|---------------|----------------------|---------|
| M = 40 M☉, R = 0.14 pc, T = 60K | Cannot satisfy all three | BE equations over-constrained |
| BE with T = 60K, R = 0.14 pc | M ~ 2-3 M☉ | Mass too small |
| BE with T = 60K, M = 40 M☉ | R ~ 0.8 pc | Radius too large |
| K&I 2000 at n ~ 10⁴ cm⁻³ | T_eq ~ 7-10 K | Too cold, then Jeans unstable |

## Literature Approaches

### 1. Ballone et al. (2018) - Simple Uniform Sphere (RECOMMENDED FOR FIRST RUN)

```json
{
    "sample": "uniform_sphere",
    "M_cloud": 40.0,
    "R_cloud": 0.14,
    "T_cloud": 60.0,
    "n_particles": 100000,
    "self_gravity": false,
    "M_BH": 57000,
    "impact_parameter": 0.0
}
```

**Why this works:**
- Matches observed parameters directly
- Self-gravity OFF (tidal forces dominate, so cloud self-gravity negligible)
- No need for hydrostatic equilibrium (cloud is NOT in equilibrium!)
- Reproduces observed velocity gradient

**Physics rationale:** The observed cloud has M_vir/M ~ 100, meaning kinetic energy >> gravitational energy. Self-gravity is irrelevant.

### 2. Oka et al. (2017) - Gaussian N-body Style

```json
{
    "sample": "gaussian_cloud",
    "sigma_r": 0.2,
    "sigma_v": 1.43,
    "N_particles": 1000,
    "M_BH": 100000,
    "initial_position": [9.8, -0.65, 0.0],
    "initial_velocity": [-8.19, 0.4, 0.0]
}
```

**Why this works:**
- Initially virialized (σ_v set for virial equilibrium)
- Gaussian density profile shows tidal stretching naturally
- Orbital parameters fitted to reproduce PV diagram

### 3. Pre-Encounter Cloud with Cooling (ADVANCED)

If you want to simulate the FULL evolution (pre-encounter → encounter → observed state):

```json
{
    "sample": "isothermal_bonnor_ebert",
    "T_cloud": 15.0,
    "n_center": 300,
    "xi_s": 3.0,
    "self_gravity": true,
    "cooling": "ki2000",
    "M_BH": 100000,
    "encounter_start_time": 1.0
}
```

**This gives:**
- M ~ 200-400 M☉ (larger pre-encounter cloud)
- R ~ 3 pc (pre-compression)
- T = 15 K (K&I 2000 equilibrium)
- Stable against Jeans instability (see stability table below)
- Cloud compresses and heats during encounter

**Expected evolution:**
1. Pre-encounter: T = 15 K, R ~ 3 pc, hydrostatic
2. Pericenter passage: Tidal compression, shock heating
3. Post-encounter: T ~ 60 K, R ~ 0.2 pc (matches observations!)

## Stability Table for Pre-Encounter Clouds

| n_center (cm⁻³) | T_eq (K) | M_BE (M☉) | M_J(T_eq) | M/M_J | Status |
|-----------------|----------|-----------|-----------|-------|--------|
| 100 | 80 | 1276 | 2192 | 0.58 | ✓ STABLE |
| **200** | **45** | **453** | **654** | **0.69** | **✓ BEST** |
| 300 | 28 | 224 | 262 | 0.86 | ⚠ RISKY |
| 500 | 17 | 112 | 96 | 1.16 | ✗ UNSTABLE |

**Maximum safe density: n ~ 200 cm⁻³** for K&I 2000 cooling with self-gravity.

## Recommended Path Forward

### Option A: Quick Start (Match Observations Directly)

1. Use **uniform sphere** with Ballone parameters
2. **Turn OFF self-gravity** (tidal dominated regime)
3. **Turn OFF cooling** (or use isothermal EOS at T = 60 K)
4. Focus on tidal dynamics and velocity structure

```bash
# Create simple config
cat > config/presets/ballone_simple.json << 'EOF'
{
    "sample": {
        "name": "uniform_cloud",
        "M_cloud": 40.0,
        "R_cloud": 0.14,
        "T_cloud": 60.0,
        "N": 50000,
        "mu": 2.33
    },
    "gravity": {
        "enabled": false
    },
    "external_force": {
        "type": "point_mass",
        "M_BH": 57000,
        "position": [0.5, 0, 0],
        "softening": 0.01
    },
    "physics": {
        "gamma": 1.0001,
        "isothermal": true
    }
}
EOF
```

### Option B: Full Physical Simulation (COMPLEX)

1. Start with **low-density BE sphere** (n ~ 200 cm⁻³)
2. Enable **K&I 2000 cooling** + **self-gravity**
3. **Verify Jeans stability** before adding IMBH
4. Introduce IMBH at t = 1 Myr (after thermal equilibration)

This approach is physically richer but requires careful setup to avoid the Jeans collapse problem you encountered.

## Why Your Current Approach Failed

Your phase 2.5 simulation failed because:

1. **n_center = 3000 cm⁻³** → T_eq = 7 K (K&I 2000)
2. **M_cloud = 22.5 M☉** → exceeds M_J(7K) = 10 M☉
3. **M/M_J = 2.74** → gravitational collapse before IMBH encounter

The K&I 2000 cooling drove the cloud to thermal equilibrium, which destabilized it gravitationally.

## The Key Insight

**The Oka/Ballone simulations work because they ignore the thermal equilibrium problem entirely:**

- Ballone: T = 60 K **imposed**, not from cooling equilibrium
- Self-gravity OFF (tidal dominated)
- No need to worry about Jeans stability

**Your approach** (K&I 2000 + BE sphere + self-gravity) is more physically complete but hits the fundamental conflict: K&I 2000 equilibrium is too cold for high-density self-gravitating clouds.

## Summary: Choose Your Level of Physical Realism

| Approach | Self-Gravity | Cooling | Jeans Stable? | Matches Obs? |
|----------|--------------|---------|---------------|--------------|
| Ballone (simple) | OFF | OFF | N/A | ✓ Yes |
| Oka N-body | OFF | OFF | N/A | ✓ Yes |
| Low-density BE | ON | K&I 2000 | ✓ if n<300 | Partially |
| High-density BE | ON | K&I 2000 | ✗ Collapses | No |

**Recommendation:** Start with Option A (Ballone-style uniform sphere, no self-gravity). Once you reproduce the basic dynamics, you can add physics incrementally.

---
*Analysis based on Oka et al. (2017) Nature Astronomy and Ballone et al. (2018) MNRAS*
