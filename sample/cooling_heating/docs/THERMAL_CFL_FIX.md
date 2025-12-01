# Thermal CFL Condition - Analysis and Fix

## Problem Identified

The timestep (`dt`) was NOT properly constrained by the **thermal/cooling timescale**, leading to potential instabilities in thermal evolution.

## Root Cause Analysis

### 1. Current Timestep Calculation (`src/timestep.cpp`)

The code only considers:
```cpp
dt_sound = CFL_sound * (h / c_s)     // Sound crossing time
dt_force = CFL_force * sqrt(h / |a|)  // Acceleration time
dt = min(dt_sound, dt_force)
```

**MISSING:** Thermal timescale constraint!

### 2. Physical Cooling Timescale

For CNM at n = 10 cm⁻³ (from K&I Figure 1d):
- **Cooling time**: τ_cool ~ 10⁵ years
- **In code units** (1 code time = 1 Myr): τ_cool ~ 0.1 code units

### 3. Thermal CFL Condition

For explicit thermal evolution with cooling/heating:
```
dt_thermal < CFL_thermal * τ_thermal
```

where:
- `τ_thermal = thermalRelaxationTime = 0.1` (code units)
- `CFL_thermal ~ 0.1` (typical for explicit methods)

Therefore:
```
dt < 0.1 * 0.1 = 0.01 code units
```

### 4. Configuration Problem

**Old config:**
```json
"thermalRelaxationTime": 0.1,   // Correct physical timescale
"dtMax": 1.0,                    // ❌ TOO LARGE! (10x thermal timescale!)
```

This violates the thermal CFL condition by factor of 100!

## The Fix

**Updated config:**
```json
"thermalRelaxationTime": 0.1,   // Physical cooling time
"dtMax": 0.01,                  // ✅ Respects thermal CFL (0.1 * τ_thermal)
```

Now:
```
dt ≤ 0.01 = 0.1 * τ_thermal  ✅ Stable thermal evolution
```

## Why This Matters

### Without Thermal CFL

With `dt = 1.0` (old config):
```
Δu = (u_eq - u_current) * dt / τ_thermal
    = (u_eq - u_current) * 1.0 / 0.1
    = 10 * (u_eq - u_current)
```

**Result:** Energy changes by 10× the deficit in ONE timestep!
- Massive overshooting
- Numerical oscillations
- Potential instability

### With Proper Thermal CFL

With `dt = 0.01` (new config):
```
Δu = (u_eq - u_current) * 0.01 / 0.1
    = 0.1 * (u_eq - u_current)
```

**Result:** Energy changes by only 10% of deficit per timestep
- Smooth relaxation
- Numerically stable
- Accurate thermal evolution

## Physical Interpretation

Think of thermal relaxation like a damped spring:
- `τ_thermal` = damping timescale
- Must resolve the exponential decay: `T(t) = T_eq + (T_0 - T_eq) * exp(-t/τ)`
- Need ~10 timesteps per e-folding time
- Therefore: `dt ~ 0.1 * τ_thermal`

## Comparison with Other CFL Conditions

| Constraint | Formula | Typical CFL | This Test |
|------------|---------|-------------|-----------|
| Sound waves | dt < CFL_s * h/c_s | 0.3 | 0.3 |
| Acceleration | dt < CFL_a * sqrt(h/a) | 0.25 | 0.25 |
| **Thermal** | **dt < CFL_th * τ_th** | **0.1** | **0.1** |

The **most restrictive** condition sets the timestep.

For this CNM test:
- Sound crossing: ~0.001 (very restrictive due to small h)
- Force time: N/A (no gravity, weak pressure forces)
- **Thermal time: 0.01** ← **This should limit dt**

## Implementation Note

**Current code limitation:** `src/timestep.cpp` does NOT check thermal timescale.

**Workaround:** Set `dtMax` manually to respect thermal CFL.

**Better solution (future):** Add thermal constraint to timestep calculation:
```cpp
void TimeStep::calculation(Simulation* sim) {
    // ... existing sound and force constraints ...
    
    // Add thermal constraint if cooling enabled
    if (sim->has_cooling()) {
        real dt_thermal = c_thermal * sim->get_thermal_timescale();
        dt = min(dt, dt_thermal);
    }
    
    sim->set_dt(dt);
}
```

## Verification

After fix, check that simulation respects thermal timescale:

```bash
# Run simulation
make -f sample/cooling_heating/Makefile.cooling cnm_run

# Check actual timestep used
tail -20 sample/cooling_heating/results/cnm_relaxation/energy.dat

# Should see dt ~ 0.001-0.01, NOT dt ~ 1.0
```

Expected output:
```
loop: 100, time: 0.95, dt: 0.0095, num: 1000
```

## Summary

✅ **Fixed:** `dtMax = 1.0` → `dtMax = 0.01`

✅ **Reason:** Respect thermal CFL condition: `dt < 0.1 * τ_thermal`

✅ **Result:** Stable, accurate thermal relaxation to K&I equilibrium

The simulation will now properly resolve the thermal evolution timescale!
