# BNS Cocoon Shock Breakout and Fast Ejecta

Special relativistic Godunov SPH simulations of cocoon shock breakout and fast ejecta 
from binary neutron star (BNS) mergers, with application to GRB 170817A and predictions 
for future events.

## Quick Start

The C++ code includes built-in 2D sample initialization for BNS cocoon simulations:

```bash
# Build for 2D
cd build
cmake .. -DDIM=2
make -j8

# Run quick test
cd ..
make -f simulations/astrophysics/bns_cocoon/Makefile.bns_cocoon cocoon_2d_test

# Visualize results  
make -f simulations/astrophysics/bns_cocoon/Makefile.bns_cocoon cocoon_2d_viz

# Or full simulation
make -f simulations/astrophysics/bns_cocoon/Makefile.bns_cocoon cocoon_2d_run
```

## Scientific Background

This project models outflows from binary neutron star mergers in 1D and 2D using SR-GSPH, 
focusing on two related scenarios:

1. **Jet-driven cocoon**: A relativistic jet pushing through merger ejecta, producing 
   cocoon shock breakout with gamma/X-ray emission
2. **Fast dynamical ejecta**: Fast tail of merger ejecta interacting with the external 
   medium, driving a relativistic blast wave

The simulations use special relativity (no dynamical GR), with SR-GSPH capturing strong 
shocks, steep density gradients, and mildly to ultra-relativistic flows.

## Initial Condition Construction

### Method 1: Direct Use of GR Simulation Data

Based on **Radice et al. (2018, ApJ 869:130)** public dataset on Zenodo:
"Binary Neutron Star Mergers: Mass Ejection, Electromagnetic Counterparts, and Nucleosynthesis"

**Available data files per binary model:**
- `outflow.txt`: Angle-integrated outflow rate and cumulative ejecta mass vs time
- `profile.txt`: Time-integrated ejecta profiles as function of polar angle
- `hist_vinf.dat`: Histogram of ejecta mass vs asymptotic velocity
- `hist_vinf_theta.h5`: Ejecta mass as function of velocity and polar angle
- `hist_ye_theta.h5`: Electron fraction distribution

**Reconstruction procedure:**

1. Select BNS model (masses, EOS) representative of GW170817-like systems
2. Extract differential mass distribution dM/(dΩ dv∞)(θ, v) from `hist_vinf_theta.h5`
3. Assume homologous expansion at reference time t₀: r ≈ v·t₀
4. Compute shell density: ρ(θⱼ, rₖ) ≈ ΔM(θⱼ, vₖ) / V_shell
5. Sample SRGSPH particles to reproduce ρ(r,θ) and v(r,θ)
6. For 1D: Average over angle using `profile.txt`

### Method 2: Analytic Fits to GR Simulations

**Hotokezaka et al. (2013, 2015):**
- Kinetic energy distribution: E(≥β) ∝ β⁻⁰·⁵ up to β ≈ 0.4
- Typical average velocity: β ≈ 0.2 for 1.4-1.4 M☉ with APR4 EOS

**Bauswein, Goriely & Janka (2013):**
- Angular distribution concentrated within ~60° from orbital plane
- Parametric form: ρ(θ) ∝ 1 + A·cos^n(θ)

**Radice et al. (2016, 2018, 2020), Foucart et al. (2024):**
- Detailed angular distributions and EOS/mass-ratio dependence

## Observables Extraction

### 1. Shock Dynamics
```
Γ_sh(t) = 1 / √(1 - v_sh²/c²)
Γ_bo ≈ Γ_sh(t_bo)  [breakout Lorentz factor]
```

### 2. Internal Energy
```
E_int(t) = Σ_shell mᵢ uᵢ
E_int,bo = E_int(t_bo)  [at breakout]
```

### 3. Optical Depth and Breakout
```
τ(t) = ∫_{r_sh}^∞ κ ρ(r) dr
Breakout: τ(t_bo) ≈ 1, r_sh(t_bo) = R_bo
```

### 4. Radiation Temperature
```
e_rad ≈ E_int,bo / V_shell
e_rad = a T_bo⁴
T_bo ≈ (E_int,bo / (a V_shell))^(1/4)
T_obs ≈ Γ_bo T_bo / (1 + z)
```

### 5. Flash Properties
```
E_rad ≈ f_rad E_int,bo
E_γ,iso ≈ (4π / ΔΩ) E_rad
t_curv ≈ R_bo / (2c Γ_bo²)
t_flash ≈ max(t_curv, ΔR_shell/c)
L_peak ≈ E_rad / t_flash
```

## Directory Structure

```
simulations/astrophysics/bns_cocoon/
├── README.md                          # This file
├── Makefile.bns_cocoon               # Build and run targets
├── config/
│   └── presets/
│       ├── cocoon_1d_fiducial.json   # 1D cocoon shock breakout
│       ├── cocoon_2d_axisym.json     # 2D axisymmetric cocoon
│       ├── ejecta_1d_spherical.json  # 1D fast ejecta
│       ├── ejecta_2d_axisym.json     # 2D fast ejecta-ISM
│       └── grb170817a_like.json      # GRB 170817A-matched parameters
├── data/
│   └── radice2018/                   # Radice+2018 Zenodo data (user-supplied)
│       ├── README.md                 # Instructions for downloading
│       └── [model directories]
├── scripts/
│   ├── ejecta_models/
│   │   ├── __init__.py
│   │   ├── radice_loader.py          # Load Radice+2018 HDF5 data
│   │   ├── homologous_expansion.py   # Convert v→r at reference time
│   │   ├── analytic_profiles.py      # Hotokezaka, Bauswein fits
│   │   ├── particle_sampler.py       # Generate SRGSPH particles
│   │   └── generate_ic.py            # Main IC generation script
│   ├── analysis/
│   │   ├── shock_tracker.py          # Track shock position and Lorentz factor
│   │   ├── optical_depth.py          # Compute optical depth τ(t)
│   │   ├── breakout_detector.py      # Detect shock breakout (τ=1)
│   │   ├── observables.py            # Compute E_rad, T_obs, L_peak
│   │   └── grb_comparison.py         # Compare to GRB 170817A
│   └── visualization/
│       ├── ejecta_structure.py       # Plot ρ(r,θ), v(r,θ)
│       ├── shock_evolution.py        # Shock position, Γ vs time
│       ├── breakout_flash.py         # Breakout observables
│       └── animate_cocoon.py         # 2D animation
├── results/                          # Simulation output
└── docs/
    ├── PHYSICS_NOTES.md              # Detailed physics derivations
    ├── RADICE_DATA_FORMAT.md         # Radice+2018 data description
    └── REFERENCES.md                 # Bibliography
```

## Quick Start

### 1. Download Radice+2018 Data

```bash
# See data/radice2018/README.md for Zenodo download instructions
```

### 2. Generate Initial Conditions

```bash
# From GR simulation data
python scripts/ejecta_models/generate_ic.py \
    --source radice2018 \
    --model SFHo_M140140 \
    --t0 1.0 \
    --output config/presets/ejecta_radice.json

# From analytic fit
python scripts/ejecta_models/generate_ic.py \
    --source analytic \
    --profile hotokezaka2013 \
    --total-mass 0.01 \
    --avg-velocity 0.2 \
    --output config/presets/ejecta_analytic.json
```

### 3. Run Simulation

```bash
# Build for 1D or 2D
make -j4 DIM=1  # or DIM=2

# Run cocoon simulation
make -f simulations/astrophysics/bns_cocoon/Makefile.bns_cocoon cocoon_1d_run

# Or with specific config
./build/sph simulations/astrophysics/bns_cocoon/config/presets/cocoon_1d_fiducial.json
```

### 4. Analyze Results

```bash
# Track shock and detect breakout
python scripts/analysis/shock_tracker.py results/cocoon_1d/

# Compute observables
python scripts/analysis/observables.py results/cocoon_1d/ \
    --output results/cocoon_1d/observables.json

# Compare to GRB 170817A
python scripts/analysis/grb_comparison.py results/cocoon_1d/observables.json
```

## Physical Parameters

### GW170817/GRB 170817A Reference Values

| Parameter | Value | Reference |
|-----------|-------|-----------|
| Total ejecta mass | 0.01-0.05 M☉ | Cowperthwaite+2017, Villar+2017 |
| Average ejecta velocity | 0.1-0.3 c | Mooley+2018 |
| Fast tail velocity | 0.4-0.8 c | Hotokezaka+2018 |
| Jet opening angle | ~5° | Mooley+2018 |
| Cocoon opening angle | ~30° | Nakar & Piran 2017 |
| Isotropic gamma-ray energy | ~10⁴⁶ erg | Abbott+2017 |
| Gamma-ray duration | ~2 s | Abbott+2017 |

### Simulation Parameters

| Parameter | Fiducial | Range |
|-----------|----------|-------|
| Ejecta mass M_ej | 0.01 M☉ | 0.001-0.1 M☉ |
| Reference time t₀ | 1 s | 0.1-10 s |
| Cocoon energy E_coc | 10⁵⁰ erg | 10⁴⁸-10⁵² erg |
| Jet Lorentz factor Γ_jet | 100 | 10-1000 |
| ISM density n_ISM | 10⁻³ cm⁻³ | 10⁻⁵-1 cm⁻³ |
| Opacity κ | 0.2 cm²/g | 0.1-10 cm²/g |

## References

### Primary Sources
- Radice et al. (2018), ApJ 869:130 - BNS mass ejection dataset
- Hotokezaka et al. (2013), PRD 87:024001 - Dynamical ejecta structure
- Bauswein, Goriely & Janka (2013), ApJ 773:78 - Ejecta systematics

### Cocoon and Jet Physics
- Nakar & Piran (2017), ApJ 834:28 - Cocoon emission
- Gottlieb et al. (2018), MNRAS 479:588 - Jet-cocoon structure
- Mooley et al. (2018), Nature 561:355 - GRB 170817A structure

### Kilonova and Afterglow
- Abbott et al. (2017), ApJL 848:L12 - GW170817 multi-messenger
- Villar et al. (2017), ApJL 851:L21 - Kilonova modeling
- Hajela et al. (2022), ApJL 927:L17 - Late-time X-ray excess

## License

This research code is part of the SPH simulation package. See main repository LICENSE.
