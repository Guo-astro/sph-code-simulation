# ISM Cooling Functions Implementation

This directory contains two ISM cooling function implementations suitable for different simulation needs.

## Available Implementations

### 1. Inoue & Inutsuka (2008) - Simplified Analytic Cooling
**File:** `inoue_inutsuka_cooling.{hpp,cpp}`  
**Reference:** Inoue, T. & Inutsuka, S. (2008), ApJ

#### When to Use
- **Diffuse ISM simulations** (n = 10⁻² - 10² cm⁻³)
- **Fast evaluation** required (no table lookups)
- **Simple two-phase medium** (WNM/CNM)
- **IMBH-cloud scattering** where τ_AD >> τ_dyn

#### Features
- Analytic cooling coefficient: Λ/Γ = 10⁷ exp(-114800/(T+1000)) + 0.014 √T exp(-92/T)
- Constant heating rate: Γ = 2×10⁻²⁶ erg s⁻¹
- Newton-Raphson equilibrium solver
- Thermal stability checker (Balbus criterion)
- **Valid range:** 10⁻² - 10² cm⁻³

#### Limitations
- ❌ No column density dependence
- ❌ No molecular cooling (H₂, CO)
- ❌ No chemistry tracking
- ❌ **Invalid above 10⁴ cm⁻³** (breaks down at high density)

---

### 2. Koyama & Inutsuka (2000) - Full ISM Cooling with Chemistry
**File:** `koyama_inutsuka_cooling.{hpp,cpp}`  
**Reference:** Koyama, H. & Inutsuka, S. (2000), ApJ 532, 980

#### When to Use
- **Molecular cloud simulations** (n = 10² - 10⁶ cm⁻³)
- **High-density regions** (star formation, dense cores)
- **Chemistry tracking** needed (H₂, CO, e⁻)
- **Column density dependence** matters

#### Features
- **Full cooling processes:**
  - [CII] 158 μm fine-structure
  - [OI] 63/145 μm fine-structure
  - Lyman-α
  - H₂ rovibrational cooling
  - CO rotational cooling
  - Gas-grain collisions
  
- **Full heating processes:**
  - Photoelectric heating from dust/PAHs
  - Cosmic ray ionization
  - X-ray heating
  - H₂ formation heating

- **Chemistry equilibrium:**
  - H/H₂ fraction from tabulated data
  - CO fraction
  - Electron fraction
  
- **Tabulated data:** Extracted from K&I 2000 figures
  - N_H = 10¹⁹ cm⁻² (65 density points)
  - N_H = 10²⁰ cm⁻² (65 density points)
  - Log-log interpolation in (n, N_H) space

- **Valid range:** 10⁻² - 10⁶ cm⁻³

#### Requirements
- Data files in `data/ki2000_extracted/`:
  - `f1a_temperature_N{19,20}.txt`
  - `f1a_pressure_N{19,20}.txt`
  - `f1b_x_H2_N{19,20}.txt`
  - `f1b_x_CO_N{19,20}.txt`
  - `f1b_x_e_N{19,20}.txt`

#### Fallback Mode
If tabulated data is unavailable, uses simplified approximations:
- T_eq ~ 6000 K for n < 1 cm⁻³ (WNM)
- T_eq ~ 100 (n/30)⁻⁰·³ K for n > 1 cm⁻³ (CNM)
- Basic H₂ and e⁻ estimates

---

## Comparison Table

| Feature | Inoue & Inutsuka 2008 | Koyama & Inutsuka 2000 |
|---------|----------------------|------------------------|
| **Density range** | 10⁻² - 10² cm⁻³ | 10⁻² - 10⁶ cm⁻³ |
| **Molecular cooling** | ❌ No | ✅ Yes (H₂, CO) |
| **Chemistry** | ❌ No | ✅ Yes (H₂, CO, e⁻) |
| **Column dependence** | ❌ No | ✅ Yes (N_H) |
| **Data files needed** | ❌ None | ✅ Required |
| **Evaluation speed** | ⚡ Very fast (analytic) | 🐢 Moderate (interpolation) |
| **Memory usage** | 💾 Minimal | 💾 ~10 KB (tables) |
| **Best for** | Diffuse ISM, fast sims | Molecular clouds, chemistry |

---

## Usage Examples

### C++ Code

```cpp
#include "thermal/inoue_inutsuka_cooling.hpp"
#include "thermal/koyama_inutsuka_cooling.hpp"

// Simplified cooling (diffuse ISM)
sph::thermal::InoueInutsukaCooling cooling_simple(5.0/3.0);

real n_H = 1.0;  // cm^-3
real T_eq = cooling_simple.equilibrium_temperature(n_H);
// T_eq ~ 4900 K (WNM)

// Full cooling (molecular clouds)
sph::thermal::KoyamaInutsukaCooling cooling_full(
    5.0/3.0,                        // gamma
    "data/ki2000_extracted",        // data directory
    1.7,                            // G0 (FUV field)
    1.0e-17                         // zeta_CR
);

real n_H = 1000.0;    // cm^-3
real N_H = 1e19;      // cm^-2 (column density)
real T_eq = cooling_full.equilibrium_temperature(n_H, N_H);
real x_H2 = cooling_full.h2_fraction(n_H, N_H);
// T_eq ~ 20 K, x_H2 ~ 0.01 (molecular cloud)
```

### SPH Integration

Both classes provide `cooling_rate_sph()` for direct SPH integration:

```cpp
real dudt = cooling.cooling_rate_sph(
    rho,         // density [code units]
    u,           // specific internal energy [code units]
    N_H,         // column density [cm^-2] (K&I 2000 only)
    dt,          // timestep [code units]
    n_to_cgs,    // conversion factor to cm^-3
    u_to_cgs,    // conversion factor to erg/g
    t_to_cgs     // conversion factor to seconds
);
```

---

## Visualization

### Plot Inoue & Inutsuka (2008)
```bash
python scripts/plot_inoue_inutsuka_cooling.py --output-dir results/cooling_analysis
```

Generates:
- `cooling_coefficients.png` - Λ/Γ and Λ(T)
- `net_cooling_rate.png` - Cooling rate vs temperature
- `equilibrium_curve.png` - T_eq and P_eq vs density
- `phase_diagram.png` - Thermal stability map
- `cooling_timescale.png` - Cooling time in (n,T) space
- `cooling_summary.png` - 4-panel summary

### Plot Koyama & Inutsuka (2000)
```bash
python scripts/plot_koyama_inutsuka_cooling.py --output-dir results/ki2000_analysis
```

Generates:
- `ki2000_equilibrium_curves.png` - T_eq and P_eq for different N_H
- `ki2000_chemistry.png` - H₂, CO, e⁻ fractions
- `ki2000_vs_inoue2008.png` - Comparison with simplified formula
- `ki2000_summary.png` - 4-panel comprehensive summary

---

## Physical Parameters

### Inoue & Inutsuka (2008)
```
Γ = 2×10⁻²⁶ erg s⁻¹          (constant heating)
Λ/Γ = 10⁷ exp(-114800/(T+1000)) + 0.014 √T exp(-92/T)
m_n = 1.27 m_p                (mean particle mass)
γ = 5/3                       (adiabatic index)
```

**Equilibrium phases:**
- WNM: n ~ 0.57 cm⁻³, T ~ 6000 K, P/k ~ 3500 K cm⁻³
- CNM: n ~ 30 cm⁻³, T ~ 100 K, P/k ~ 3000 K cm⁻³

### Koyama & Inutsuka (2000)
```
G0 = 1.7                      (FUV field in Habing units)
ζ_CR = 1×10⁻¹⁷ s⁻¹           (cosmic ray ionization rate)
N_H = 10¹⁹ or 10²⁰ cm⁻²      (column density)
```

**Cooling mechanisms:**
- [CII] 158 μm: Dominant in WNM (T ~ 6000 K)
- [OI] 63/145 μm: Important in CNM (T ~ 100 K)
- H₂ cooling: Dominates at n > 1000 cm⁻³
- CO cooling: Important at n > 10⁴ cm⁻³

---

## Selection Guide

### Use Inoue & Inutsuka (2008) if:
✅ Simulating **diffuse ISM** (n < 100 cm⁻³)  
✅ Need **fast computation** (real-time, large N)  
✅ Two-phase medium is sufficient  
✅ Chemistry is not tracked  
✅ Column density effects negligible

### Use Koyama & Inutsuka (2000) if:
✅ Simulating **molecular clouds** (n > 100 cm⁻³)  
✅ Need **chemistry** (H₂, CO, ionization)  
✅ **Column density** matters (shielding)  
✅ Willing to pay interpolation cost  
✅ Have extracted tabulated data

### For High Density (n > 10⁶ cm⁻³):
⚠️ **Neither formula is valid!**  
Need specialized cooling:
- Glover & Mac Low (2007) - chemistry network
- GRACKLE library - full non-equilibrium chemistry
- Opacity-limited cooling for optically thick regions

---

## Implementation Notes

### Thread Safety
Both implementations are **thread-safe** for read-only operations:
- ✅ `equilibrium_temperature()`, `cooling_rate_sph()` safe
- ⚠️ Constructor modifies internal state (call before parallel region)

### Performance
- **Inoue & Inutsuka 2008:** ~100 ns per evaluation (analytic)
- **Koyama & Inutsuka 2000:** ~1 μs per evaluation (table lookup + interpolation)

### Memory
- **Inoue & Inutsuka 2008:** Negligible (~100 bytes)
- **Koyama & Inutsuka 2000:** ~10 KB (2 tables × 65 points × 5 arrays)

---

## References

1. **Inoue, T. & Inutsuka, S. (2008)**  
   "Two-Fluid MHD Simulations of Converging HI Flows in the Interstellar Medium"  
   ApJ

2. **Koyama, H. & Inutsuka, S. (2000)**  
   "Two-dimensional Simulations of Thermal Instability in the Interstellar Medium"  
   ApJ 532, 980

3. **Wolfire, M. G. et al. (1995)**  
   "The Neutral Atomic Phases of the Interstellar Medium"  
   ApJ 443, 152

4. **Bakes, E. L. O. & Tielens, A. G. G. M. (1994)**  
   "The Photoelectric Heating Mechanism for Very Small Graphitic Grains and PAHs"  
   ApJ 427, 822
