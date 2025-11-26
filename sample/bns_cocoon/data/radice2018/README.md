# Radice et al. (2018) Ejecta Data

## Data Source

**Paper**: Radice, D., Perego, A., Hotokezaka, K., et al. (2018)  
"Binary Neutron Star Mergers: Mass Ejection, Electromagnetic Counterparts, and Nucleosynthesis"  
ApJ 869:130, [arXiv:1809.11161](https://arxiv.org/abs/1809.11161)

**Dataset**: [Zenodo DOI: 10.5281/zenodo.1298474](https://zenodo.org/record/1298474)

## Download Instructions

```bash
# Option 1: Download via Zenodo web interface
# Visit: https://zenodo.org/record/1298474
# Download the full dataset archive

# Option 2: Use zenodo_get (Python tool)
pip install zenodo_get
zenodo_get 10.5281/zenodo.1298474 -o .

# Option 3: Direct wget (if individual files are available)
# Check Zenodo page for direct links
```

## Available Models

The dataset includes multiple BNS models with different:
- Total masses (e.g., 2.7-2.8 M☉)
- Mass ratios (q = M₂/M₁)
- Equations of state (EOS): SFHo, DD2, LS220, BHBlp, etc.

### Recommended Models for GW170817

For GW170817-like simulations, consider models with:
- Total mass ≈ 2.74 M☉ (chirp mass constraint)
- Mass ratio q ≈ 0.7-1.0
- Stiff EOS (SFHo, DD2) for typical NS radii

## File Format

Each model directory contains:

### `outflow.txt`
```
# Column 1: time [ms]
# Column 2: outflow rate [M_sun/s]
# Column 3: cumulative ejecta mass [M_sun]
```

### `profile.txt`
```
# Column 1: cos(theta)
# Column 2: time-integrated mass in angular bin [M_sun]
# Column 3: average velocity in angular bin [c]
# Column 4: average Y_e in angular bin
```

### `hist_vinf.dat`
```
# Column 1: v_inf [c] (asymptotic velocity)
# Column 2: dM/dv [M_sun/c] (differential mass)
```

### `hist_vinf_theta.h5` (HDF5 format)
```python
# Structure:
# 'v_inf': velocity bins [c]
# 'cos_theta': angular bins
# 'dM_dv_dOmega': differential mass dM/(dv dΩ) [M_sun/c/sr]
```

### `hist_ye_theta.h5` (HDF5 format)
```python
# Structure:
# 'Y_e': electron fraction bins
# 'cos_theta': angular bins
# 'dM_dYe_dOmega': differential mass dM/(dY_e dΩ) [M_sun/sr]
```

## Usage Example

```python
import h5py
import numpy as np

# Load velocity-angle distribution
with h5py.File('hist_vinf_theta.h5', 'r') as f:
    v_inf = f['v_inf'][:]           # velocity bins
    cos_theta = f['cos_theta'][:]   # angular bins
    dM_dv_dOmega = f['dM_dv_dOmega'][:, :]  # 2D mass distribution

# Convert to theta
theta = np.arccos(cos_theta)

# Integrate over angle for 1D profile
dv = np.diff(v_inf, prepend=0)
dOmega = 2 * np.pi * np.abs(np.diff(cos_theta, prepend=1))
M_v = np.sum(dM_dv_dOmega * dOmega[:, np.newaxis], axis=0) * dv
```

## Physical Units

- **Mass**: Solar masses (M☉)
- **Velocity**: Speed of light (c)
- **Time**: Milliseconds (ms)
- **Angles**: Radians or cos(θ)

## Coordinate System

- **θ = 0**: Polar axis (perpendicular to orbital plane)
- **θ = π/2**: Equatorial plane (orbital plane)
- Dynamical ejecta concentrated near equator (θ ≈ π/2)

## Notes

1. Data represents **dynamical ejecta** only (not disk wind ejecta)
2. Homologous expansion assumed at late times (t > 10 ms)
3. Asymptotic velocities valid for r >> initial separation
4. Y_e important for nucleosynthesis but not needed for hydrodynamics

## Citation

If using this data, please cite:
```bibtex
@article{Radice2018,
  author = {Radice, David and Perego, Albino and Hotokezaka, Kenta and others},
  title = {Binary Neutron Star Mergers: Mass Ejection, Electromagnetic Counterparts, 
           and Nucleosynthesis},
  journal = {ApJ},
  volume = {869},
  pages = {130},
  year = {2018},
  doi = {10.3847/1538-4357/aaf054}
}
```
