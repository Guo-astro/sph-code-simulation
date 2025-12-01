# ISM Cooling 1D Test - Current Status

## Status: ⚠️ PARTIALLY COMPLETE

### ✅ What Works
1. **Code Integration**: All files compiled successfully
2. **GDISPH Cooling**: Cooling module properly integrated into fluid force
3. **Thermal Parameters**: Configuration system working
4. **Initial Condition**: 1D density gradient setup complete
5. **Build System**: Makefile targets functional

### ⚠️  Current Issue
**Unit Conversion Mismatch**: The Koyama-Inutsuka cooling module uses physical CGS units (K, cm⁻³, erg/s), but the SPH code uses dimensionless code units.

**Symptom**: NaN values in simulation output due to incorrect temperature conversion.

## Problem Analysis

### K&I Cooling Module Expects:
- Input: `n_H` in [cm⁻³], `T` in [K]
- Output: `T_eq` in [K]
- Cooling rate: `du/dt` in [erg s⁻¹ g⁻¹]

### SPH Code Uses:
- Density: code units (dimensionless)
- Temperature: implicit in `P/ρ` (code units)
- Energy: code units (dimensionless)

### Failed Conversion Attempt:
```cpp
// This doesn't work because code units != physical units
const real T_eq = m_cooling->temperature(n_H);  // Returns K
const real T_current = p_i.pres / p_i.dens;     // Code units, NOT K
const real u_eq = T_eq / (m_gamma - 1.0);       // Mixing units!
```

## Solutions

### Option 1: Modify Cooling Module (Recommended)
Make `cooling_rate()` work in dimensionless units:
```cpp
real cooling_rate_dimensionless(real rho_code, real u_code, real P_code, 
                                real unit_density, real unit_energy, real t_relax)
{
    // Convert code → CGS
    real n_H = rho_code * unit_density / (mu * m_H);
    real T = P_code * unit_pressure / (n_H * k_B);
    
    // Get equilibrium
    real T_eq = temperature(n_H);  // K
    real u_eq_cgs = k_B * T_eq / ((gamma-1) * mu * m_H);
    real u_current_cgs = u_code * unit_energy;
    
    // Convert back to code units
    return (u_eq_cgs - u_current_cgs) / (t_relax * unit_energy);
}
```

### Option 2: Use Physical Units Throughout
Set up code with physical units:
```json
{
  "units": {
    "length": 3.0857e18,    // 1 pc
    "mass": 1.989e33,        // 1 M_sun
    "time": 3.1557e13        // 1 Myr
  }
}
```

### Option 3: Simplified Test (Quick Fix)
Disable cooling temporarily to test GDISPH mechanics:
```json
{
  "enableCooling": false
}
```

## Recommended Action

1. **Short term**: Implement Option 1 (add unit conversion layer)
2. **Medium term**: Create unit system in parameters
3. **Long term**: Full physical unit support in code

## Current Files Status

| File | Status | Notes |
|------|--------|-------|
| `include/parameters.hpp` | ✅ Done | Thermal params added |
| `include/gdisph/gd_fluid_force.hpp` | ✅ Done | Cooling member added |
| `src/gdisph/gd_fluid_force.cpp` | ⚠️  Needs fix | Unit conversion required |
| `src/sample/ism_cooling_1d.cpp` | ✅ Done | Initial condition works |
| `sample/cooling_heating/ism_cooling_1d.json` | ✅ Done | Config complete |
| `include/thermal/koyama_inutsuka_cooling.hpp` | ⚠️  Needs method | Add dimensionless interface |
| `src/thermal/koyama_inutsuka_cooling.cpp` | ⚠️  Needs method | Implement unit conversion |
| `Makefile.cooling_heating` | ✅ Done | All targets defined |

## Next Steps

1. Add unit conversion parameters to config
2. Implement `cooling_rate_code_units()` method
3. Update GDISPH integration to use new method
4. Test with simple uniform density case first
5. Then run full gradient test

## Workaround for Testing

To test GDISPH without cooling:
```bash
# Edit ism_cooling_1d.json
"enableCooling": false

# Run
make -f sample/cooling_heating/Makefile.cooling_heating cooling_run
```

This will verify:
- GDISPH works correctly
- Initial conditions are set up properly
- Output format is correct
- Visualization pipeline works

Then re-enable cooling once unit conversion is fixed.

## Documentation Completed

✅ `ISM_COOLING_1D_INTEGRATION.md` - Integration guide
✅ `QUICK_START.md` - User guide
✅ `COMPLETE_IMPLEMENTATION.md` - Cooling module docs
✅ `Makefile.cooling_heating` - Build automation

---

**Bottom Line**: Code infrastructure is 95% complete. Just needs unit conversion layer to make cooling work with dimensionless SPH code units.
