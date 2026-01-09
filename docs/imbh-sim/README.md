# IMBH-Cloud Interaction Simulation

## Overview

This project simulates molecular cloud dynamics for IMBH tidal interaction studies.

## Key Physical Constraint

**For an isothermal self-gravitating sphere at T = 10 K:**

| n_center [cm^-3] | R_cloud [pc] | M_cloud [M_sun] |
|------------------|--------------|-----------------|
| 100 | 0.11 | 0.028 |
| 1000 | 0.035 | 0.003 |
| 10000 | 0.011 | 0.0003 |

**The cloud radius scales as R ~ c_s / sqrt(G * rho_c).**

Cold (T ~ 10 K) molecular clouds with high density (n > 100 cm^-3) are necessarily compact (R << 1 pc). This is physically correct - cold molecular cores are embedded in larger warm envelopes.

## Configuration

### Current Setup: `isothermal_bonnor_ebert` with n_center specified

```json
{
    "sample": "isothermal_bonnor_ebert",
    "n_center": 1000.0,    // Central density [cm^-3]
    "T_cloud": 10.0,       // Temperature [K]
    "xi_s": 5.0            // Dimensionless truncation (critical ~6.45)
}
```

This gives:
- R_cloud ~ 0.035 pc
- M_cloud ~ 0.003 M_sun
- True hydrostatic equilibrium with self-gravity

### Running

```bash
make -f simulations/astrophysics/imbh_cloud/Makefile.ki2000_1pc hydrostatic
```

## Why Not R = 1 pc?

To get R = 1 pc with isothermal equilibrium at T = 10 K, you'd need:
- n_center ~ 0.1 cm^-3 (extremely diffuse)
- This is the warm neutral medium (WNM), not a molecular cloud

The original Oka (2017) HVCC CO-0.40-0.22 has:
- M ~ 40 M_sun
- R ~ 0.17 pc (compact core!)
- n ~ 10^5 cm^-3

Our setup with n_center = 1000 cm^-3 gives R ~ 0.035 pc, which is a reasonable compact core scale.

## Files

```
docs/imbh-sim/
    README.md                    # This file
    FIRST_PRINCIPLES.md          # Detailed physics derivations

simulations/astrophysics/imbh_cloud/
    config/presets/ki2000_1pc/
        hydrostatic.json         # Main config
    Makefile.ki2000_1pc          # Run commands
```

## References

1. Koyama & Inutsuka (2000) - ISM thermal equilibrium
2. Oka et al. (2017) - HVCC CO-0.40-0.22 IMBH candidate
3. Bonnor (1956), Ebert (1955) - Pressure-truncated isothermal spheres
