# Chemistry Network Test Case
## ISM Cooling and Heating (Koyama & Inutsuka 2000)

**Status**: ✅ Complete  
**Purpose**: Full reproduction of ISM thermal and chemical equilibrium (Figure 1)  
**Directory**: `sample/cooling_heating/`

---

## Quick Start

```bash
# View available targets
make -f sample/cooling_heating/Makefile.cooling_heating cooling_help

# Reproduce Figure 1 from the paper
make -f sample/cooling_heating/Makefile.cooling_heating cooling_reproduce

# Compare with original
make -f sample/cooling_heating/Makefile.cooling_heating cooling_compare

# Clean up
make -f sample/cooling_heating/Makefile.cooling_heating cooling_clean
```

## What This Test Case Does

This test case implements the **complete ISM chemistry network** from Koyama & Inutsuka (2000, ApJ 533, 793) and reproduces their Figure 1, which shows:

1. **Panel (a)**: Equilibrium temperature and pressure vs density
2. **Panel (b)**: Chemical fractions (e⁻, H₂, CO, HI) vs density
3. **Panel (c)**: Individual heating and cooling rates vs density
4. **Panel (d)**: Physical timescales vs density

## Physics Included

### Chemistry (6 species):
- H, H⁺, e⁻ (ionization/recombination)
- H₂ (formation/dissociation with self-shielding)
- CO (formation/destruction)
- Trace elements: C⁺, O I

### Heating (5 processes):
- Photoelectric (grains/PAHs)
- Cosmic rays
- X-rays
- H₂ formation
- H₂ photodissociation

### Cooling (6 major processes):
- Fine-structure lines (CII, OI, FeII, SiII)
- Lyman-α
- H₂ rovibrational
- CO rotational
- CO vibrational
- Grain collisions

## Output

All figures saved to `sample/cooling_heating/results/`:

```
figure1_reproduction.png    # 4-panel combined figure
f1a_reproduction.png        # Temperature & Pressure
f1b_reproduction.png        # Chemical fractions
f1c_reproduction.png        # Heating/Cooling rates
f1d_reproduction.png        # Timescales
```

## Implementation Notes

- **Total Code**: ~1,500 lines of Python
- **Runtime**: 2-5 minutes
- **Dependencies**: numpy, scipy, matplotlib
- **All formulas**: Transcribed directly from paper appendix

## Documentation

- `README.md` - User guide and physics overview
- `IMPLEMENTATION_SUMMARY.md` - Detailed technical summary
- `COMPLETION_REPORT.md` - Full implementation verification

## References

Koyama, H. & Inutsuka, S. 2000, ApJ, 533, 793  
(Plus 8 supporting references for individual processes)

---

**Note**: This is a standalone Python implementation for thermal equilibrium calculations. It is separate from the main SPH code but could be integrated for time-dependent ISM simulations.
