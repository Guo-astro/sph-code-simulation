# Grad-H Correction Study# Grad-H Correction Study# Grad-H Correction Study# Grad-H Correction Study: Essential for Hydrostatic Equilibrium



## Overview



This study demonstrates that the **grad-h correction term (Ω)** is essential for maintaining hydrostatic equilibrium in SPH simulations with variable smoothing length. Without this correction, self-gravitating clouds experience artificial core collapse even when they should remain stable.## Overview



## Physics



### The Grad-H Correction FactorThis study demonstrates that the **grad-h correction term (Ω)** is essential for maintaining hydrostatic equilibrium in SPH simulations with variable smoothing length. Without this correction, self-gravitating clouds experience artificial core collapse even when they should remain stable.## OverviewThis test case demonstrates that the grad-h correction term (Ω) is **essential** for preventing artificial core collapse in hydrostatic equilibrium simulations using SPH methods with variable smoothing length.



When using SPH with adaptive smoothing length *h*, the density summation:



```## Physics

ρ_i = Σ_j m_j W(|r_i - r_j|, h_i)

```



creates an implicit dependence `ρ_i ↔ h_i` because `h_i = η (m_i/ρ_i)^(1/D)`.### The Grad-H Correction FactorThis study demonstrates that the **grad-h correction term** is essential for maintaining hydrostatic equilibrium in SPH simulations. Without this correction, self-gravitating clouds experience artificial core collapse even when they should remain in stable equilibrium.## Theory



The grad-h correction factor accounts for this:



```When using SPH with adaptive smoothing length $h$, the density summation:

Ω_i = [1 + (h_i / D ρ_i) Σ_j m_j ∂W_ij/∂h_i]^(-1)

```



where *D* is the number of dimensions.$$\rho_i = \sum_j m_j W(|\mathbf{r}_i - \mathbf{r}_j|, h_i)$$## Theory### The Grad-H Correction Factor



### Physical Interpretation



| Condition | Behavior |creates an implicit dependence $\rho_i \leftrightarrow h_i$ because $h_i = \eta (m_i/\rho_i)^{1/D}$.

|-----------|----------|

| **With Ω ≠ 1** | Pressure gradient properly accounts for spatial variations in smoothing length. Hydrostatic balance maintained. |

| **Without Ω (Ω = 1)** | Spurious forces arise from inconsistency between density and pressure gradients. Leads to artificial core collapse. |

The grad-h correction factor accounts for this:The grad-h correction factor accounts for variable smoothing length in the SPH density summation:When using SPH with adaptive smoothing length $h$, the density summation:

### Analytic Solution: Lane-Emden Equation



The equilibrium state is described by the Lane-Emden equation for polytropic index *n*:

$$\Omega_i = \left[1 + \frac{h_i}{D \rho_i} \sum_j m_j \frac{\partial W_{ij}}{\partial h_i}\right]^{-1}$$

```

(1/ξ²) d/dξ (ξ² dθ/dξ) = -θ^n

```

where $D$ is the number of dimensions.```$$\rho_i = \sum_j m_j W(|\mathbf{r}_i - \mathbf{r}_j|, h_i)$$

For an ideal gas with γ = 5/3, we have n = 1.5 and:

- Surface radius: ξ₁ ≈ 3.6538

- Physical density profile: ρ(r) = ρ_c θ(ξ)^(3/2)

### Physical InterpretationΩ_i = 1 / (1 + (h_i / D ρ_i) × dρ_i/dh_i)

The interactive visualization overlays this analytic solution to verify equilibrium maintenance.



## Quick Start

| Condition | Behavior |```becomes inconsistent with the pressure gradient if $h$ varies spatially. The grad-h correction factor addresses this:

```bash

# Run 2-case comparison (Ω ON vs OFF) + generate interactive HTML|-----------|----------|

make -f sample/gradh_study/Makefile.gradh_study gradh_all

| **With Ω ≠ 1** | Pressure gradient properly accounts for spatial variations in smoothing length. Hydrostatic balance maintained. |

# Open interactive viewer

open sample/gradh_study/figures_v2/interactive_viewer.html| **Without Ω (Ω = 1)** | Spurious forces arise from inconsistency between density and pressure gradients. Leads to artificial core collapse. |

```

Where:$$\Omega_i = \frac{1}{1 + \frac{h_i}{D \rho_i} \frac{\partial \rho_i}{\partial h_i}}$$

## Directory Structure

### Analytic Solution: Lane-Emden Equation

```

sample/gradh_study/- `h_i` is the smoothing length of particle i

├── Makefile.gradh_study         # Main build targets

├── README.md                    # This fileThe equilibrium state is described by the Lane-Emden equation for polytropic index $n$:

├── config/

│   └── presets/                 # All config files (SSOT)- `ρ_i` is the densitywhere $D$ is the number of dimensions and $\frac{\partial \rho_i}{\partial h_i}$ is computed from the kernel derivative with respect to $h$.

│       ├── selfgrav_gradh_hk.json           # Ω ON + Hernquist-Katz

│       ├── selfgrav_no_gradh_hk.json        # Ω OFF + Hernquist-Katz$$\frac{1}{\xi^2}\frac{d}{d\xi}\left(\xi^2\frac{d\theta}{d\xi}\right) = -\theta^n$$

│       ├── selfgrav_gradh_wendland.json     # Ω ON + Wendland C4

│       ├── selfgrav_no_gradh_wendland.json  # Ω OFF + Wendland C4- `D` is the number of spatial dimensions

│       └── *_long.json                      # Long-run variants (t=150)

├── scripts/For an ideal gas with $\gamma = 5/3$, we have $n = 1.5$ and:

│   ├── gen_html_viewer.py               # Interactive HTML generator (MAIN)

│   ├── compare_gradh_selfgrav_v2.py     # Legacy matplotlib comparison- Surface radius: $\xi_1 \approx 3.6538$- `dρ/dh` is computed from the density summation derivative### Physical Interpretation

│   └── kernel_gravity_comparison.py     # 4-case energy analysis

├── results/- Physical density profile: $\rho(r) = \rho_c \theta(\xi)^{3/2}$

│   ├── case1_gradh_hk/          # Results: Ω ON + H-K (stable)

│   ├── case2_no_gradh_hk/       # Results: Ω OFF + H-K (collapse)

│   ├── case3_gradh_wendland/    # Results: Ω ON + Wendland

│   ├── case4_no_gradh_wendland/ # Results: Ω OFF + WendlandThe interactive visualization overlays this analytic solution to verify equilibrium maintenance.

│   ├── kernel_test/             # 4-case kernel comparison

│   └── long_test/               # Long-duration collapse runs**Without this correction (Ω = 1):**- **With grad-h correction (Ω ≠ 1)**: The pressure force properly accounts for spatial variations in the smoothing length, maintaining hydrostatic balance.

└── figures_v2/

    └── interactive_viewer.html  # Main visualization output## Quick Start

```

- The pressure gradient becomes inconsistent with the density estimate

## Test Matrix

```bash

| Case | Grad-H (Ω) | Grav Softening | Expected Result |

|------|------------|----------------|-----------------|# Run 2-case comparison (Ω ON vs OFF) + generate interactive HTML- Spurious forces arise that push particles toward density peaks- **Without grad-h correction (Ω = 1)**: Spurious forces arise from the inconsistency between density estimation and pressure gradients. In self-gravitating systems or any density gradient, this leads to artificial core collapse.

| 1    | ON         | Hernquist-Katz | ✅ Stable       |

| 2    | OFF        | Hernquist-Katz | ❌ Collapse     |make -f sample/gradh_study/Makefile.gradh_study gradh_all

| 3    | ON         | Wendland C4    | ✅ Stable       |

| 4    | OFF        | Wendland C4    | ❓ KEY TEST     |- Results in artificial core collapse in equilibrium configurations



**HYPOTHESIS**: True kernel-convolved gravity (Wendland C4) where gravitational softening matches the SPH kernel may eliminate the need for Ω correction.# Open interactive viewer



## Available Make Targetsopen sample/gradh_study/figures_v2/interactive_viewer.html## Directory Structure



```bash```

# Quick 2-case study (recommended)

make -f sample/gradh_study/Makefile.gradh_study gradh_run     # Run case1 + case2**With the correction:**

make -f sample/gradh_study/Makefile.gradh_study viz_html      # Generate HTML viewer

make -f sample/gradh_study/Makefile.gradh_study gradh_all     # Run + visualize## Directory Structure



# 4-case kernel comparison- Density and pressure gradient calculations are self-consistent```

make -f sample/gradh_study/Makefile.gradh_study kernel_all    # All 4 cases + viz

```

# Individual cases

make -f sample/gradh_study/Makefile.gradh_study case1         # Ω ON + H-Ksample/gradh_study/- Equilibrium configurations remain stablesample/gradh_study/

make -f sample/gradh_study/Makefile.gradh_study case2         # Ω OFF + H-K

make -f sample/gradh_study/Makefile.gradh_study case3         # Ω ON + Wendland├── Makefile.gradh_study         # Main build targets

make -f sample/gradh_study/Makefile.gradh_study case4         # Ω OFF + Wendland

├── README.md                    # This file- Essential for any simulation involving variable smoothing lengths├── README.md                    # This file

# Long runs (t=150)

make -f sample/gradh_study/Makefile.gradh_study long_all├── gradh_study.json             # Active config (overwritten by presets)



# Cleanup├── config/├── gradh_study.json             # Active configuration (overwritten by presets)

make -f sample/gradh_study/Makefile.gradh_study clean

│   └── presets/

# Help

make -f sample/gradh_study/Makefile.gradh_study help│       ├── selfgrav_gradh_hk.json           # Ω ON + Hernquist-Katz## Quick Start├── Makefile.gradh_study         # Build targets

```

│       ├── selfgrav_no_gradh_hk.json        # Ω OFF + Hernquist-Katz

## Interactive Viewer Features

│       ├── selfgrav_gradh_wendland.json     # Ω ON + Wendland C4├── config/

The HTML viewer (`figures_v2/interactive_viewer.html`) provides:

│       ├── selfgrav_no_gradh_wendland.json  # Ω OFF + Wendland C4

1. **Animated particle scatter plots** - Color-coded by density

2. **Radial density profile** - With Lane-Emden n=1.5 analytic overlay (cyan line)│       └── *_long.json                      # Long-run variants (t=150)```bash│   └── presets/

3. **Ω distribution vs radius** - Shows deviation from Ω=1

4. **Energy conservation** - E/|E₀| over time├── scripts/

5. **First-principles derivation** - LaTeX-rendered theory explanation

│   ├── gen_html_viewer.py               # Interactive HTML generator (MAIN)# Run the complete study (simulations + visualization + interactive viewer)│       ├── hydrostatic_gsph_with_gradh.json    # Periodic box, with Ω

## Key Results

│   ├── compare_gradh_selfgrav_v2.py     # Legacy matplotlib comparison

- **With Ω correction**: Cloud maintains hydrostatic equilibrium, density profile matches Lane-Emden solution

- **Without Ω correction**: Central density increases ~20%+ over time, artificial collapse occurs│   └── kernel_gravity_comparison.py     # 4-case energy analysismake -f sample/gradh_study/Makefile.gradh_study gradh_all│       ├── hydrostatic_gsph_no_gradh.json      # Periodic box, without Ω

- **Energy conservation**: Proper Ω correction ensures E_total remains constant

├── results/

## Initial Conditions

│   ├── case1_gradh_hk/          # Results: Ω ON + H-K (stable)│       ├── selfgrav_gsph_with_gradh.json       # Self-gravitating, with Ω

- **Relaxed polytrope**: 10k particles from `simulations/astrophysics/imbh_cloud/results/relaxation/10k/GSPH/snapshot_0020.csv`

- **Polytropic index**: n = 1.5 (γ = 5/3 ideal gas)│   ├── case2_no_gradh_hk/       # Results: Ω OFF + H-K (collapse)

- **Gravitational softening**: Hernquist-Katz with ε = h/2

│   ├── case3_gradh_wendland/    # Results: Ω ON + Wendland# Open the interactive comparison viewer│       └── selfgrav_gsph_no_gradh.json         # Self-gravitating, without Ω

## References

│   ├── case4_no_gradh_wendland/ # Results: Ω OFF + Wendland

- Springel & Hernquist (2002): "Cosmological SPH simulations"

- Hopkins (2013): "A general class of Lagrangian SPH methods"│   ├── kernel_test/             # 4-case kernel comparisonopen sample/gradh_study/figures/interactive_viewer.html├── scripts/

- Hernquist & Katz (1989): "TREESPH"

│   └── long_test/               # Long-duration collapse runs

└── figures_v2/```│   ├── compare_gradh.py                  # Comparison visualization

    └── interactive_viewer.html  # Main visualization output

```│   ├── compare_gradh_selfgrav.py         # Self-gravitating comparison



## Test Matrix## Directory Structure│   ├── animate_gradh_comparison.py       # Animation generator



| Case | Grad-H (Ω) | Grav Softening | Expected Result |│   └── generate_paper_figures.py         # Publication-ready figures

|------|------------|----------------|-----------------|

| 1    | ON         | Hernquist-Katz | ✅ Stable       |```├── results/                     # Simulation output

| 2    | OFF        | Hernquist-Katz | ❌ Collapse     |

| 3    | ON         | Wendland C4    | ✅ Stable       |sample/gradh_study/├── figures/                     # Generated plots

| 4    | OFF        | Wendland C4    | ❓ KEY TEST     |

├── Makefile.gradh_study           # Build and run commands└── animations/                  # Generated animations

**HYPOTHESIS**: True kernel-convolved gravity (Wendland C4) where gravitational softening matches the SPH kernel may eliminate the need for Ω correction.

├── README.md                       # This file```

## Available Make Targets

├── gradh_study.json               # Active config (copied from presets)

```bash

# Quick 2-case study (recommended)├── config/## Usage

make -f sample/gradh_study/Makefile.gradh_study gradh_run     # Run case1 + case2

make -f sample/gradh_study/Makefile.gradh_study viz_html      # Generate HTML viewer│   └── presets/

make -f sample/gradh_study/Makefile.gradh_study gradh_all     # Run + visualize

│       ├── selfgrav_gsph_with_gradh.json    # WITH grad-h (stable)### Prerequisites

# 4-case kernel comparison

make -f sample/gradh_study/Makefile.gradh_study kernel_all    # All 4 cases + viz│       └── selfgrav_gsph_no_gradh.json      # WITHOUT grad-h (collapse)



# Individual cases├── scripts/1. **Build configuration**: 

make -f sample/gradh_study/Makefile.gradh_study case1         # Ω ON + H-K

make -f sample/gradh_study/Makefile.gradh_study case2         # Ω OFF + H-K│   └── compare_gradh_selfgrav.py  # Visualization script   - For periodic box tests: DIM=2

make -f sample/gradh_study/Makefile.gradh_study case3         # Ω ON + Wendland

make -f sample/gradh_study/Makefile.gradh_study case4         # Ω OFF + Wendland├── results/                       # Simulation output   - For self-gravitating tests: DIM=3



# Long runs (t=150)│   ├── selfgrav_with_gradh/       # Results with grad-h

make -f sample/gradh_study/Makefile.gradh_study long_all

│   └── selfgrav_no_gradh/         # Results without grad-h2. **For self-gravitating tests**: Requires relaxed initial conditions from IMBH cloud relaxation:

# Cleanup

make -f sample/gradh_study/Makefile.gradh_study clean└── figures/                       # Generated visualizations   ```bash



# Help    ├── selfgrav_gradh_central_density.png   make -f simulations/astrophysics/imbh_cloud/Makefile.relaxation relax_oneshot SIZE=61k METHOD=gsph

make -f sample/gradh_study/Makefile.gradh_study help

```    ├── selfgrav_gradh_radius.png   ```



## Interactive Viewer Features    ├── gradh_comparison.gif



The HTML viewer (`figures_v2/interactive_viewer.html`) provides:    ├── interactive_viewer.html    # Interactive HTML viewer### Running the Study



1. **Animated particle scatter plots** - Color-coded by density    └── frames/                    # Individual animation frames

2. **Radial density profile** - With Lane-Emden n=1.5 analytic overlay (cyan line)

3. **Ω distribution vs radius** - Shows deviation from Ω=1```#### Periodic Box Hydrostatic Test (2D)

4. **Energy conservation** - E/|E₀| over time

5. **First-principles derivation** - LaTeX-rendered theory explanation



## Key Results## Make Targets```bash



- **With Ω correction**: Cloud maintains hydrostatic equilibrium, density profile matches Lane-Emden solution# Run both with and without grad-h simulations

- **Without Ω correction**: Central density increases ~20%+ over time, artificial collapse occurs

- **Energy conservation**: Proper Ω correction ensures E_total remains constant| Target | Description |make -f sample/gradh_study/Makefile.gradh_study gradh_hydro_run



## Initial Conditions|--------|-------------|



- **Relaxed polytrope**: 10k particles from `simulations/astrophysics/imbh_cloud/results/relaxation/10k/GSPH/snapshot_0020.csv`| `gradh_help` | Show help message |# Generate comparison plots

- **Polytropic index**: n = 1.5 (γ = 5/3 ideal gas)

- **Gravitational softening**: Hernquist-Katz with ε = h/2| `gradh_run` | Run both WITH and WITHOUT grad-h simulations |make -f sample/gradh_study/Makefile.gradh_study gradh_hydro_viz



## References| `gradh_run_with` | Run only WITH grad-h simulation |



- Springel & Hernquist (2002): "Cosmological SPH simulations: A hybrid multiphase model for star formation"| `gradh_run_no` | Run only WITHOUT grad-h simulation |# Generate animation

- Hopkins (2013): "A general class of Lagrangian smoothed particle hydrodynamics methods"

- Hernquist & Katz (1989): "TREESPH: A unification of SPH with the hierarchical tree method"| `gradh_viz` | Generate comparison plots and interactive viewer |make -f sample/gradh_study/Makefile.gradh_study gradh_hydro_animate


| `gradh_all` | Run everything (simulations + visualizations) |

| `gradh_clean` | Clean all results |# Run all (simulation + visualization + animation)

make -f sample/gradh_study/Makefile.gradh_study gradh_hydro_all

## Output Files```



### Static Plots#### Self-Gravitating Cloud Test (3D)

- `selfgrav_gradh_central_density.png` - Central density evolution over time

- `selfgrav_gradh_radius.png` - Cloud radius (90% mass) evolution```bash

# Run both cases (requires 3D build and relaxed IC)

### Animationmake -f sample/gradh_study/Makefile.gradh_study gradh_selfgrav_run

- `gradh_comparison.gif` - Side-by-side animation showing both cases

# Generate comparison plots

### Interactive Viewermake -f sample/gradh_study/Makefile.gradh_study gradh_selfgrav_viz

- `interactive_viewer.html` - HTML page with:

  - Frame slider navigation# Run all

  - Play/Pause controlsmake -f sample/gradh_study/Makefile.gradh_study gradh_selfgrav_all

  - Speed adjustment (0.5x to 8x)```

  - Keyboard shortcuts (Space, Arrow keys, Home, End)

  - Theory explanation panel#### Paper Figures



## Expected Results```bash

# Generate publication-ready figures

### With Grad-H Correction (Blue)make -f sample/gradh_study/Makefile.gradh_study gradh_paper_figs

- Central density remains constant at initial value```

- Cloud radius stays at ~1.0 (normalized)

- Radial density profile maintains Lane-Emden shape### Configuration Parameter

- System stays in hydrostatic equilibrium

The `useGradH` parameter in JSON config files controls the grad-h correction:

### Without Grad-H Correction (Orange)

- Central density increases dramatically (artificial collapse)```json

- Cloud radius contracts{

- Density profile develops a sharp central cusp    "SPHType": "gsph",

- System collapses even though it should be in equilibrium    "useGradH": true,    // true (default) or false

    ...

## Prerequisites}

```

1. **DIM=3 build:**

   ```bash## Expected Results

   cd build && cmake -DSPH_DIM=3 .. && make -j8

   ```### With Grad-H Correction (Stable)



2. **Relaxed initial conditions:**- Central density remains constant (ρc/ρc,0 ≈ 1.0)

   ```bash- Particle distribution stays uniform

   make -f simulations/astrophysics/imbh_cloud/Makefile.relaxation relax_oneshot SIZE=10k METHOD=gsph- Energy is well-conserved

   ```- Cloud radius remains stable in self-gravitating tests



3. **Python packages:**### Without Grad-H Correction (Collapse)

   ```bash

   pip install numpy matplotlib pillow- Central density increases rapidly

   ```- Particles drift toward center

- Energy conservation degrades

## Code Implementation- Artificial core collapse occurs



The grad-h correction is controlled by the `useGradH` parameter in the JSON config:## Quantitative Metrics



```jsonThe study produces several diagnostic plots:

{

    "useGradH": true,   // Enable grad-h correction (default)1. **Central Density Evolution**: Shows ρc(t)/ρc(0) vs time

    // or2. **Particle Distribution**: Side-by-side snapshots at multiple times

    "useGradH": false   // Disable for demonstration3. **Density Profile**: Radial density distribution ρ(r)

}4. **Energy Conservation**: (E - E₀)/|E₀| vs time

```5. **Grad-H Factor Distribution**: Histogram of Ω values



The correction is applied in `src/gsph/g_pre_interaction.cpp`:## References



```cpp- Springel, V. & Hernquist, L. (2002). "Cosmological smoothed particle hydrodynamics simulations: the entropy equation." MNRAS 333, 649.

// With grad-h: compute Ω from density derivative- Hopkins, P. F. (2013). "A general class of Lagrangian smoothed particle hydrodynamics methods." MNRAS 428, 2840.

if (m_use_gradh) {- Price, D. J. (2012). "Smoothed particle hydrodynamics and magnetohydrodynamics." J. Comput. Phys. 231, 759.

    real omega_inv = 1.0 + (p[i].sml / (DIM * p[i].dens)) * drhodh;

    p[i].omega = 1.0 / omega_inv;## Paper Outline

} else {

    // Without grad-h: Ω = 1 (no correction)This study supports a paper demonstrating:

    p[i].omega = 1.0;

}1. **Introduction**: Variable-h SPH requires careful treatment of kernel gradients

```2. **Theory**: Derivation of the grad-h correction factor Ω

3. **Numerical Tests**: 

## References   - Periodic box hydrostatic equilibrium

   - Self-gravitating polytropic cloud

- Hopkins, P. F. (2013). "A general class of Lagrangian smoothed particle hydrodynamics methods..."4. **Results**: 

- Springel, V. & Hernquist, L. (2002). "Cosmological smoothed particle hydrodynamics simulations..."   - Without Ω: artificial core collapse

- Price, D. J. (2012). "Smoothed particle hydrodynamics and magnetohydrodynamics" (JCPH)   - With Ω: stable equilibrium maintained

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
make -f simulations/astrophysics/imbh_cloud/Makefile.relaxation relax_oneshot SIZE=61k METHOD=gsph
```

### Python dependencies
```bash
pip install numpy matplotlib
```
