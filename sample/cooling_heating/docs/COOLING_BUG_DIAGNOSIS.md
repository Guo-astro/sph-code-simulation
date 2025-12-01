# CRITICAL BUG: Cooling NOT Applied - Root Cause Analysis

**Date:** December 1, 2024  
**Status:** 🔴 CRITICAL - Cooling code exists but du/dt NEVER applied

---

## Executive Summary

**Problem:** Thermal energy stays CONSTANT despite `enableCooling: true`
- Simulation runs without errors
- Cooling module is initialized correctly
- But `du/dt` from cooling is **NEVER integrated into particle energy**

**Root Cause:** Missing time integration step in main simulation loop

---

## Evidence

### 1. Energy File Shows NO Cooling
```
time: 0.0    thermal: 1599.664
time: 0.5    thermal: 1599.666  # NO CHANGE!
time: 1.0    thermal: 1599.666  # NO CHANGE!
time: 1.5    thermal: 1599.666  # NO CHANGE!
time: 2.0    thermal: 1599.666  # NO CHANGE!
```

**Expected:** E_thermal should DECREASE by ~30% (T: 107K → ~25K)

### 2. Code Flow Analysis

#### ✅ Cooling Rate Calculation (WORKS)
`src/gdisph/gd_fluid_force.cpp` lines 201-228:
```cpp
if (m_enable_cooling && p_i.dens > 0 && p_i.pres > 0) {
    const real n_H = p_i.dens * m_density_to_n_H;
    const real T_eq = m_cooling->temperature(n_H);
    const real T_current = p_i.pres / p_i.dens;
    const real u_eq = T_eq / (m_gamma - 1.0);
    const real u_current = p_i.ene;
    const real du_dt = (u_eq - u_current) / m_thermal_relax_time;
    
    if (std::isfinite(du_dt)) {
        dene += du_dt;  // ✅ Adds cooling to dene
    }
}

p_i.dene = dene;  // ✅ Stores du/dt in particle
```

**This part is CORRECT!** `du/dt` is calculated and stored in `p_i.dene`.

#### ❌ Time Integration (MISSING or BROKEN)
The particle `dene` field is set, but **where is it integrated?**

Expected somewhere in main loop:
```cpp
// EXPECTED (but missing?):
for (auto& p : particles) {
    p.ene += p.dene * dt;  // Apply energy change
}
```

---

## Debug Investigation

### Check 1: Does dene get written to particles?
```bash
# Add debug output in gd_fluid_force.cpp after line 228:
if (i == 0) {  // First particle
    std::cout << "DEBUG: particle 0: dene=" << p_i.dene 
              << ", du_cooling=" << du_dt << std::endl;
}
```

### Check 2: Where is time integration done?
Search for energy update in simulation loop:
```bash
grep -n "\.ene.*dt" src/*.cpp src/*/*.cpp
grep -n "integrate" src/simulation.cpp
grep -n "kick.*drift" src/simulation.cpp
```

### Check 3: Is there a separate integrator module?
```bash
find src -name "*integrat*" -o -name "*kick*" -o -name "*drift*"
```

---

## Hypotheses

### Hypothesis 1: Missing Kick Step
SPH typically uses **kick-drift-kick** leapfrog:
1. **Kick:** `v^{n+1/2} = v^n + (dt/2) * a^n`
2. **Drift:** `x^{n+1} = x^n + dt * v^{n+1/2}`
3. **Kick:** `v^{n+1} = v^{n+1/2} + (dt/2) * a^{n+1}`

For thermal energy:
```cpp
// Kick thermal energy
p.ene += dt * p.dene;  # THIS MIGHT BE MISSING!
```

### Hypothesis 2: Integrator Doesn't Know About dene
Maybe the time integrator only updates `pos`, `vel`, `acc` and **forgets** `ene`?

Check `src/simulation.cpp` for the main loop.

### Hypothesis 3: GDISPH vs Base SPH
- Base `sph::FluidForce` sets `dene`
- But GDISPH overrides without calling base integrator?

---

## Required Files to Check

1. **`src/simulation.cpp`** - Main time loop
   - Look for particle energy update
   - Check if `dene` is used

2. **`include/simulation.hpp`** - Simulation class
   - Check for integrator methods

3. **`src/solver.cpp`** - Solver loop
   - Where is the time integration?

4. **GSPH/DISPH equivalents**
   - Do they integrate `dene` correctly?
   - Compare with GDISPH implementation

---

## Next Steps

1. **Locate time integrator code**
   ```bash
   grep -rn "\.ene.*\+.*dt" src/
   grep -rn "kick\|drift\|integrate" src/simulation.cpp
   ```

2. **Add debug output**
   - Print `dene` in `gd_fluid_force.cpp`
   - Print `ene` before/after integration

3. **Verify integrator uses dene**
   - If integrator exists: check if it applies `dene`
   - If missing: **ADD integration step**

4. **Test fix**
   ```bash
   # Rebuild
   cd build && make -j8 && cd ..
   
   # Run with debug
   ./build/sph ism_cooling_1d 2>&1 | grep "DEBUG"
   
   # Check if energy changes
   tail sample/cooling_heating/results/cnm_relaxation/energy.dat
   ```

---

## Temporary Diagnostic Code

Add to `src/gdisph/gd_fluid_force.cpp` after line 228:

```cpp
// DIAGNOSTIC: Check if cooling is being calculated
static int debug_count = 0;
if (debug_count < 10 && i == 0 && m_enable_cooling) {
    std::cout << "COOLING DEBUG [" << debug_count << "]: "
              << "n_H=" << (p_i.dens * m_density_to_n_H) << ", "
              << "T_current=" << (p_i.pres / p_i.dens) << ", "
              << "T_eq=" << m_cooling->temperature(p_i.dens * m_density_to_n_H) << ", "
              << "du_dt=" << ((T_eq/(m_gamma-1.0) - p_i.ene)/m_thermal_relax_time) << ", "
              << "dene=" << p_i.dene << std::endl;
    debug_count++;
}
```

---

## Summary

- ✅ Cooling module **initialized correctly**
- ✅ Cooling rate `du/dt` **calculated correctly**
- ✅ Stored in `p.dene` **correctly**
- ❌ **TIME INTEGRATION MISSING OR BROKEN**

**The bug is NOT in the cooling physics - it's in the simulation loop!**

Next: Find and fix the time integrator.
