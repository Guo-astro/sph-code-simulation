# ISM Cooling and Heating Chemistry Network

Complete implementation of the thermal and chemical network described in **Koyama & Inutsuka (2000, ApJ 533, 793)** for reproducing Figure 1.

## Overview

This module implements the full chemistry network from the paper's Appendix, including:

### Chemistry
- **H ionization/recombination**: UV, cosmic ray, X-ray, and collisional ionization
- **H₂ formation/dissociation**: Grain catalysis, associative detachment, photodissociation with self-shielding, collisional dissociation
- **CO formation/destruction**: Simplified C⁺ → CO chemistry with photodissociation

### Heating Processes
- Photoelectric heating from small grains and PAHs (Bakes & Tielens 1994)
- Cosmic ray ionization heating
- Soft X-ray heating (with column density attenuation)
- H₂ formation heating
- H₂ photodissociation heating

### Cooling Processes
- Fine-structure lines: C II (158 μm), O I (63 μm), Fe II (26 μm), Si II (35 μm)
- Hydrogen Lyman-α cooling
- H₂ rovibrational line cooling (Hollenbach & McKee 1979, updated with Galli & Palla 1998)
- CO rotational and vibrational cooling
- Gas-grain collisional cooling

## Files

```
simulations/astrophysics/cooling_heating/
├── scripts/
│   ├── chemistry_network.py       # Full chemistry network implementation
│   ├── thermal_equilibrium.py     # Equilibrium solver
│   ├── reproduce_figure1.py       # Main script to generate Figure 1
│   └── run_reproduction.sh        # Convenience script
├── results/
│   ├── figure1_reproduction.png   # Combined 4-panel figure
│   ├── f1a_reproduction.png       # Panel (a): T and P vs n
│   ├── f1b_reproduction.png       # Panel (b): Chemical fractions
│   ├── f1c_reproduction.png       # Panel (c): Heating/cooling rates
│   └── f1d_reproduction.png       # Panel (d): Timescales
└── README.md                       # This file
```

## Figure 1 Description

Figure 1 from Koyama & Inutsuka (2000) shows the equilibrium thermal and chemical state of the ISM as a function of density:

- **(a) Temperature and Pressure**: Shows equilibrium T and P for two column densities (10¹⁹ and 10²⁰ cm⁻²)
- **(b) Chemical Fractions**: Electron fraction (xₑ), H₂ fraction (x₂), CO fraction (xCO), and HI fraction
- **(c) Heating and Cooling Rates**: Individual heating processes (PE, CR, XR, H₂) and cooling processes (CII, OI, Lyman-α, H₂, CO, grain collisions)
- **(d) Physical Timescales**: Cooling time, recombination time, free-fall time, and H₂ formation time

## Dependencies

```bash
pip install numpy scipy matplotlib
```

## Usage

### Quick Run

```bash
cd simulations/astrophysics/cooling_heating/scripts
./run_reproduction.sh
```

### Manual Run

```bash
cd simulations/astrophysics/cooling_heating/scripts
python3 reproduce_figure1.py
```

### Output

The script generates:
- `../results/figure1_reproduction.png` - Complete 4-panel figure
- `../results/f1{a,b,c,d}_reproduction.png` - Individual panels for detailed comparison

## Implementation Notes

### Equilibrium Solver

The solver finds thermal and chemical equilibrium by:

1. **Chemical Equilibrium**: At fixed (n, T), solve the system:
   - H ionization equilibrium: ζ·xHI = α·xₑ²·n
   - H₂ formation/destruction balance
   - CO formation/destruction balance

2. **Thermal Equilibrium**: At fixed n, find T where:
   - Net heating = Total heating - Total cooling = 0

3. **Iteration**: Compute equilibrium curves over density range 0.1-10⁶ cm⁻³

### Key Physics

- **Self-shielding**: H₂ photodissociation includes self-shielding via column density
- **UV attenuation**: Photoionization attenuated by HI column density
- **Two-phase medium**: Natural emergence of cold (CNM) and warm (WNM) phases
- **Molecular cloud conditions**: At n > 10³ cm⁻³, H₂ and CO become important coolants

### Trace Element Abundances

Following Wolfire et al. (1995):
- He/H = 0.10
- O/H = 4.6×10⁻⁴
- C/H = 3.0×10⁻⁴
- Si/H = 3.55×10⁻⁶
- Fe/H = 7.08×10⁻⁷

### Radiation Field

- UV field: G₀ = 1.7 (Habing units)
- Cosmic ray ionization rate: ζCR = 1.8×10⁻¹⁷ s⁻¹
- Column densities: NH = 10¹⁹ cm⁻² (solid) and 10²⁰ cm⁻² (dashed)

## Comparison with Original Figure

The original Figure 1 (f1a.ps, f1b.ps, f1c.ps, f1d.ps) shows:

1. **Thermal bistability**: Two stable branches at low densities (WNM ~8000K, CNM ~100K)
2. **Molecular transition**: Rapid H₂ formation at n ~ 100-1000 cm⁻³
3. **CO formation**: Begins at n ~ 10³-10⁴ cm⁻³
4. **Dominant coolants**: 
   - Low n: CII, OI fine-structure
   - High n: H₂ rovibrational, CO rotational

## References

- Koyama, H. & Inutsuka, S. 2000, ApJ, 533, 793
- Bakes, E. L. O. & Tielens, A. G. G. M. 1994, ApJ, 427, 822
- Wolfire, M. G. et al. 1995, ApJ, 443, 152
- Hollenbach, D. & McKee, C. F. 1979, ApJS, 41, 555
- Hollenbach, D. & McKee, C. F. 1989, ApJ, 342, 306
- Tielens, A. G. G. M. & Hollenbach, D. 1985, ApJ, 291, 722
- Galli, D. & Palla, F. 1998, A&A, 335, 403

## Contact

For issues or questions about this implementation, see the main repository README.
