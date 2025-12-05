# Grad-H Correction Study: Essential for Hydrostatic Equilibrium

This test case demonstrates that the grad-h correction term (Ω) is **essential** for preventing artificial core collapse in hydrostatic equilibrium simulations using SPH methods with variable smoothing length.

## Theory

### The Grad-H Correction Factor

When using SPH with adaptive smoothing length $h$, the density summation:

$$\rho_i = \sum_j m_j W(|\mathbf{r}_i - \mathbf{r}_j|, h_i)$$

becomes inconsistent with the pressure gradient if $h$ varies spatially. The grad-h correction factor addresses this:

$$\Omega_i = \frac{1}{1 + \frac{h_i}{D \rho_i} \frac{\partial \rho_i}{\partial h_i}}$$

where $D$ is the number of dimensions and $\frac{\partial \rho_i}{\partial h_i}$ is computed from the kernel derivative with respect to $h$.

### Physical Interpretation

- **With grad-h correction (Ω ≠ 1)**: The pressure force properly accounts for spatial variations in the smoothing length, maintaining hydrostatic balance.

- **Without grad-h correction (Ω = 1)**: Spurious forces arise from the inconsistency between density estimation and pressure gradients. In self-gravitating systems or any density gradient, this leads to artificial core collapse.

## Directory Structure

```
sample/gradh_study/
├── README.md                    # This file
├── gradh_study.json             # Active configuration (overwritten by presets)
├── Makefile.gradh_study         # Build targets
├── config/
│   └── presets/
│       ├── hydrostatic_gsph_with_gradh.json    # Periodic box, with Ω
│       ├── hydrostatic_gsph_no_gradh.json      # Periodic box, without Ω
│       ├── selfgrav_gsph_with_gradh.json       # Self-gravitating, with Ω
│       └── selfgrav_gsph_no_gradh.json         # Self-gravitating, without Ω
├── scripts/
│   ├── compare_gradh.py                  # Comparison visualization
│   ├── compare_gradh_selfgrav.py         # Self-gravitating comparison
│   ├── animate_gradh_comparison.py       # Animation generator
│   └── generate_paper_figures.py         # Publication-ready figures
├── results/                     # Simulation output
├── figures/                     # Generated plots
└── animations/                  # Generated animations
```

## Usage

### Prerequisites

1. **Build configuration**: 
   - For periodic box tests: DIM=2
   - For self-gravitating tests: DIM=3

2. **For self-gravitating tests**: Requires relaxed initial conditions from IMBH cloud relaxation:
   ```bash
   make -f sample/imbh_cloud/Makefile.relaxation relax_oneshot SIZE=61k METHOD=gsph
   ```

### Running the Study

#### Periodic Box Hydrostatic Test (2D)

```bash
# Run both with and without grad-h simulations
make -f sample/gradh_study/Makefile.gradh_study gradh_hydro_run

# Generate comparison plots
make -f sample/gradh_study/Makefile.gradh_study gradh_hydro_viz

# Generate animation
make -f sample/gradh_study/Makefile.gradh_study gradh_hydro_animate

# Run all (simulation + visualization + animation)
make -f sample/gradh_study/Makefile.gradh_study gradh_hydro_all
```

#### Self-Gravitating Cloud Test (3D)

```bash
# Run both cases (requires 3D build and relaxed IC)
make -f sample/gradh_study/Makefile.gradh_study gradh_selfgrav_run

# Generate comparison plots
make -f sample/gradh_study/Makefile.gradh_study gradh_selfgrav_viz

# Run all
make -f sample/gradh_study/Makefile.gradh_study gradh_selfgrav_all
```

#### Paper Figures

```bash
# Generate publication-ready figures
make -f sample/gradh_study/Makefile.gradh_study gradh_paper_figs
```

### Configuration Parameter

The `useGradH` parameter in JSON config files controls the grad-h correction:

```json
{
    "SPHType": "gsph",
    "useGradH": true,    // true (default) or false
    ...
}
```

## Expected Results

### With Grad-H Correction (Stable)

- Central density remains constant (ρc/ρc,0 ≈ 1.0)
- Particle distribution stays uniform
- Energy is well-conserved
- Cloud radius remains stable in self-gravitating tests

### Without Grad-H Correction (Collapse)

- Central density increases rapidly
- Particles drift toward center
- Energy conservation degrades
- Artificial core collapse occurs

## Quantitative Metrics

The study produces several diagnostic plots:

1. **Central Density Evolution**: Shows ρc(t)/ρc(0) vs time
2. **Particle Distribution**: Side-by-side snapshots at multiple times
3. **Density Profile**: Radial density distribution ρ(r)
4. **Energy Conservation**: (E - E₀)/|E₀| vs time
5. **Grad-H Factor Distribution**: Histogram of Ω values

## References

- Springel, V. & Hernquist, L. (2002). "Cosmological smoothed particle hydrodynamics simulations: the entropy equation." MNRAS 333, 649.
- Hopkins, P. F. (2013). "A general class of Lagrangian smoothed particle hydrodynamics methods." MNRAS 428, 2840.
- Price, D. J. (2012). "Smoothed particle hydrodynamics and magnetohydrodynamics." J. Comput. Phys. 231, 759.

## Paper Outline

This study supports a paper demonstrating:

1. **Introduction**: Variable-h SPH requires careful treatment of kernel gradients
2. **Theory**: Derivation of the grad-h correction factor Ω
3. **Numerical Tests**: 
   - Periodic box hydrostatic equilibrium
   - Self-gravitating polytropic cloud
4. **Results**: 
   - Without Ω: artificial core collapse
   - With Ω: stable equilibrium maintained
5. **Discussion**: Physical interpretation and implications
6. **Conclusion**: Grad-h correction is essential for accurate hydrostatic simulations

## Troubleshooting

### DIM mismatch error
```bash
# For 2D tests:
sed -i '' 's/#define DIM [0-9]/#define DIM 2/' include/defines.hpp
cd build && make -j8 && cd ..

# For 3D tests:
sed -i '' 's/#define DIM [0-9]/#define DIM 3/' include/defines.hpp
cd build && make -j8 && cd ..
```

### Missing relaxed IC
For self-gravitating tests, first run the relaxation:
```bash
make -f sample/imbh_cloud/Makefile.relaxation relax_oneshot SIZE=61k METHOD=gsph
```

### Python dependencies
```bash
pip install numpy matplotlib
```
