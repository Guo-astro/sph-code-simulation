# Unit System Presets

This directory contains unit preset files for common astrophysical scenarios
in special relativistic SPH simulations.

## Available Presets

### 1. `dimensionless_sr.json`
**Purpose**: Code validation and test problems  
**Scales**: All quantities dimensionless, c=1  
**Use case**: Sod shock tube, blast wave tests, Kelvin-Helmholtz tests

### 2. `neutron_star_merger.json`
**Purpose**: Binary neutron star simulations  
**Length scale**: 10 km (typical NS radius)  
**Density scale**: 10¹⁴ g/cm³ (nuclear density)  
**Time scale**: ~0.03 ms (light-crossing time)  
**Use case**: BNS mergers, post-merger disk, jet formation

### 3. `relativistic_jet.json`
**Purpose**: AGN and astrophysical jets  
**Length scale**: 1 pc  
**Density scale**: 10⁻²⁴ g/cm³ (ISM proton density)  
**Time scale**: ~3.3 years  
**Use case**: AGN jets, microquasar jets, jet-ISM interaction

### 4. `grb_afterglow.json`
**Purpose**: Gamma-ray burst afterglow dynamics  
**Length scale**: 10¹⁷ cm  
**Density scale**: 10⁻²⁴ g/cm³  
**Time scale**: ~40 days  
**Use case**: GRB external shocks, afterglow light curves

### 5. `supernova_remnant.json`
**Purpose**: Relativistic supernova remnants  
**Length scale**: 1 pc  
**Density scale**: 10⁻²⁴ g/cm³ (ISM)  
**Time scale**: ~3 kyr  
**Use case**: Young SNR, relativistic ejecta

## Usage

### In Simulation Config (JSON)

Include the `units` section in your simulation config:

```json
{
  "units": {
    "type": "RELATIVISTIC",
    "preset": "neutron_star",
    "length_km": 15.0,
    "density_g_cm3": 2.8e14
  }
}
```

Or specify custom scales directly:

```json
{
  "units": {
    "type": "RELATIVISTIC",
    "length_cm": 1.0e17,
    "density_gcm3": 1.67e-24,
    "labels": {
      "length": "10^17 cm",
      "time": "days",
      "density": "n_0"
    }
  }
}
```

### In Python Visualization

```python
from units.unit_system import RelativisticUnits

# From config file
with open('config.json') as f:
    config = json.load(f)
units = RelativisticUnits.from_config(config)

# Or create directly
units = RelativisticUnits.create_neutron_star(length_km=12, density_scale=2.8e14)

# Convert to physical units
physical_time = units.to_physical_time(code_time)  # Returns seconds
physical_length = units.to_physical_length(code_length)  # Returns cm
```

## Physical Constants

The unit system uses the following CGS constants:

| Constant | Value | Units |
|----------|-------|-------|
| Speed of light (c) | 2.998 × 10¹⁰ | cm/s |
| Gravitational constant (G) | 6.674 × 10⁻⁸ | cm³ g⁻¹ s⁻² |
| Solar mass (M☉) | 1.989 × 10³³ | g |
| Proton mass (mₚ) | 1.673 × 10⁻²⁴ | g |
| Parsec (pc) | 3.086 × 10¹⁸ | cm |
| Kiloparsec (kpc) | 3.086 × 10²¹ | cm |

## Adding New Presets

1. Create a new JSON file in this directory
2. Include `_comment` and `_description` fields for documentation
3. Specify `units.type` as "RELATIVISTIC" for SR simulations
4. Set appropriate `length_cm` and `density_gcm3` scales
5. Add `_physical_context` with typical values for reference
