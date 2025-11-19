# SR-GSPH Variable Naming and Field Usage

## Overview

SR-GSPH (Special Relativistic Godunov SPH) uses **different meanings** for the standard SPHParticle fields compared to classical SPH. This document clarifies the mapping and explains the design decisions.

## Key Principle

**SR-GSPH integrates CONSERVED variables (S, e), not PRIMITIVE variables (v, u)**

This is fundamentally different from standard SPH, which directly integrates primitive variables like velocity and internal energy.

## Field Mapping Table

| Field Name | Standard SPH | SR-GSPH | Notes |
|------------|--------------|---------|-------|
| **`S`** | *unused* | **Canonical momentum** S = γHv (per baryon) | **PRIMARY: time-integrated** |
| **`e`** | *unused* | **Canonical energy** e = γH - P/(Nc²) (per baryon) | **PRIMARY: time-integrated** |
| **`dS`** | *unused* | **Time derivative** dS/dt | Computed by force calculation |
| **`de`** | *unused* | **Time derivative** de/dt | Computed by force calculation |
| **`vel`** | Velocity (integrated) | **Primitive velocity** v | Recovered from S, NOT integrated |
| **`ene`** | Internal energy (integrated) | **Primitive internal energy** u = P/[(γ-1)n] | Recovered, NOT integrated |
| **`acc`** | Acceleration dv/dt | **ALIAS for dS** | ⚠️ MISLEADING NAME! |
| **`dene`** | Energy derivative du/dt | **ALIAS for de** | Copy for compatibility |
| **`mass`** | Particle mass | **Baryon number** ν | Constant per particle |
| **`dens`** | Mass density ρ | **Rest-frame density** n | For output only |
| **`N`** | *unused* | **Lab-frame density** N = γn | Used in primitive recovery |
| **`gamma_lor`** | *unused* | **Lorentz factor** γ = 1/√(1-v²/c²) | Derived quantity |
| **`enthalpy`** | *unused* | **Specific enthalpy** H | Derived quantity |
| **`nu`** | *unused* | **Baryon number** ν | Constant (same as mass) |

## Why This Design?

### Problem: Name Collisions

The SPHParticle structure was designed for **classical SPH**, where:
- Velocity `vel` is time-integrated
- Acceleration `acc` is the force per mass (F/m)
- Internal energy `ene` is time-integrated

In **SR-GSPH**, we need:
- Conserved variables S and e (time-integrated)
- Primitive variables v and u (recovered, not integrated)
- Time derivatives dS/dt and de/dt

### Solution: Dual Usage with Clear Documentation

**Strategy 1: Dedicated SR-GSPH fields**
- Added `S`, `e`, `dS`, `de`, `N`, `gamma_lor`, `enthalpy`, `nu` to SPHParticle
- These are the **PRIMARY** storage for SR-GSPH

**Strategy 2: Reuse standard fields for output**
- `vel` stores **recovered primitive velocity** (not time-integrated!)
- `ene` stores **recovered primitive internal energy** u
- `acc` stores **copy of dS** (for CSV/VTK output compatibility)
- `dene` stores **copy of de** (for output)

**Strategy 3: Extensive documentation**
- Comments in `particle.hpp` explain dual usage
- Comments in integration code clarify what's being integrated
- This document provides the complete mapping

## Critical Warning: Don't Mix Paradigms!

**WRONG** ❌:
```cpp
// In SR-GSPH code, DO NOT:
p[i].vel += p[i].acc * dt;  // This would integrate v, not S!
```

**CORRECT** ✅:
```cpp
// SR-GSPH integrates conserved variables:
p[i].S += p[i].dS * dt;   // Integrate S
p[i].e += p[i].de * dt;   // Integrate e

// Then recover primitives:
auto prim = PrimitiveRecovery::conserved_to_primitive(p[i].S, p[i].e, p[i].N, ...);
p[i].vel = prim.vel;  // Store recovered velocity for output
```

## Data Flow in SR-GSPH

```
1. INITIAL CONDITIONS
   └─> Set (S, e, N) from (v, P, n)
       Store primitives in (vel, ene, pres) for initial output

2. PRE-INTERACTION (each timestep)
   ├─> Compute N from kernel sum
   ├─> Read conserved (S, e) from dedicated fields
   ├─> Recover primitives (v, P, n) using quartic solver
   └─> Store primitives in (vel, pres, dens) for Riemann solver

3. FORCE CALCULATION
   ├─> Read primitives (vel, pres) from standard fields
   ├─> Solve Riemann problem at each particle pair
   ├─> Compute dS/dt and de/dt
   ├─> Store in dedicated fields (dS, de)
   └─> Copy to (acc, dene) for output compatibility

4. TIME INTEGRATION
   ├─> Integrate ONLY conserved variables: S += dS*dt, e += de*dt
   ├─> Recover primitives at half-step for position update
   ├─> Update positions: pos += v_half*dt
   ├─> Recover primitives at full step
   └─> Store in (vel, ene, pres) for output and next step

5. OUTPUT (CSV/VTK)
   ├─> Writes vel → primitive velocity v ✓
   ├─> Writes ene → primitive internal energy u ✓
   ├─> Writes acc → dS/dt (MISLEADING: not acceleration!) ⚠️
   └─> Writes dene → de/dt ✓
```

## Output File Interpretation

When reading CSV/VTK output files:

- **`vel_x, vel_y, vel_z`**: Primitive velocity v (valid for plotting)
- **`acc_x, acc_y, acc_z`**: ⚠️ NOT acceleration! Contains dS/dt (canonical momentum derivative)
- **`ene`**: Primitive internal energy u = P/[(γ-1)n] (valid for plotting)
- **`dens`**: Rest-frame density n (NOT lab-frame N!)
- **`pres`**: Pressure P (valid)

If you need the **conserved variables**:
- Canonical momentum S is NOT in CSV output (use S field in code)
- Canonical energy e is NOT in CSV output (use e field in code)

## Why Not Rename `acc` to `dS_dt_output`?

**Reason 1: Backward compatibility**
- Existing analysis scripts expect `acc` field in CSV/VTK
- Standard SPH still uses `acc` as actual acceleration
- Changing field names would break all output readers

**Reason 2: Output format consistency**
- CSV and VTK formats have fixed column headers
- Many visualization tools expect specific field names
- Renaming would require updating all parsers

**Reason 3: Field aliasing is clear in code**
- The code explicitly documents: `p_i.acc = p_i.dS;  // Copy for output`
- Comments in particle.hpp clarify the dual usage
- Integration code never uses `acc` directly in SR-GSPH path

## Recommendations for Future Work

1. **For new SR-GSPH users**: Always use the dedicated fields (S, e, dS, de) in your code
2. **For output analysis**: Remember that `acc` contains dS/dt in SR-GSPH outputs, not acceleration
3. **For extending SR-GSPH**: Follow the pattern: integrate conserved, recover primitives, store in standard fields
4. **Consider**: Adding SR-GSPH-specific output format that writes S and e directly

## References

- Kitajima, Inutsuka, Seno (2025): "Special relativistic smoothed particle hydrodynamics with improved calculation of particle volume", arXiv:2510.18251v1
- Pons, Martí, Müller (1999): "The exact solution of the Riemann problem with non-zero tangential velocities in relativistic hydrodynamics", J. Fluid Mech. 422:125-139

---
**Last Updated**: 2025-11-18  
**Author**: SR-GSPH Implementation Team
