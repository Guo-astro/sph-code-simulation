# IMBH-Cloud Interaction Simulation

SPH simulation of a molecular cloud interacting with an Intermediate-Mass Black Hole (IMBH), based on the HVCC CO-0.40-0.22 observations from Oka et al. (2017).

## Two-Phase Workflow

The simulation uses a **two-phase approach**:

### Phase 1: Relaxation with External Pressure Confinement

Create a modified isothermal sphere in hydrostatic equilibrium:

```bash
make relax
```

This phase:
1. Creates a modified isothermal sphere: `rho = rho_c / (1 + (r/r_c)^2)`
2. Adds ghost envelope particles for external pressure confinement
3. Runs GLASS pre-relaxation to uniformize particle spacing
4. Runs main relaxation (subtracts equilibrium pressure gradient)
5. Outputs relaxed configuration (ghost particles removed)

### Phase 2: IMBH Tidal Interaction

Resume from relaxed snapshot with IMBH:

```bash
make interact
```

This phase:
1. Loads the relaxed BE sphere (without ghost envelope)
2. Adds IMBH at specified distance
3. Runs tidal interaction simulation
4. Cloud freely evolves under IMBH tidal forces

### Full Workflow

Run both phases:

```bash
make run           # Full workflow
make run-quick     # Quick test with fewer relaxation steps
```

## Physical Setup

### Option A: Atomic CNM Critical BE Sphere (Recommended)

This configuration achieves a 40 M_sun **critical** Bonnor-Ebert sphere on the K&I 2000 thermal equilibrium curve using **atomic gas** (mu = 1.27).

| Parameter | Value | Derivation |
|-----------|-------|------------|
| mu | 1.27 | Atomic HI + 10% He (CNM) |
| M_cloud | 40 M_sun | Critical BE mass |
| T_cloud | 9 K | K&I equilibrium at n_edge |
| n_center | 980 cm^-3 | 14.04 x n_edge (BE profile) |
| n_edge | 70 cm^-3 | K&I curve for M_crit = 40 |
| P_ext | 630 K cm^-3 | n_edge x T |
| R_cloud | 1.2 pc | xi_s x r_0 |
| xi_s | 6.45 | Critical truncation |
| Stability | M/M_crit = 1 | Marginally stable |

**Why atomic gas?** The critical BE mass scales as M_crit ~ mu^-2. Using mu = 2.33 (molecular) would only give M_crit ~ 13 M_sun at these conditions. Atomic gas allows 40 M_sun while staying on the K&I equilibrium curve.

**Physical picture:** The cloud starts as atomic CNM. During IMBH compression, increased density and UV shielding trigger H2 formation, converting it to molecular gas (observed as CO emission).

See `docs/imbh-sim/observational_constraints.tex` Section "Critical Bonnor-Ebert Mass" for full derivation.

### Legacy: Modified Isothermal Sphere

| Parameter | Value | Notes |
|-----------|-------|-------|
| M_cloud | 40 M_sun | Observed mass |
| T_cloud | 15 K | KI2000 equilibrium |
| n_center | ~200 cm^-3 | Centrally concentrated |
| P_ext | 3000 K cm^-3 | External pressure (K&I 2000) |
| Profile | rho = rho_c / (1 + (r/r_c)^2) | Modified isothermal |

### IMBH Configuration

| Parameter | Value | Notes |
|-----------|-------|-------|
| M_BH | 10^5 M_sun | From Oka et al. (2017) |
| d_BH | 10 pc | Initial distance |
| Position | (+10, 0, 0) pc | Along +x axis |

## Why Relaxation is Necessary

The initial particle distribution doesn't exactly match the equilibrium density profile. Relaxation:

1. **External Pressure**: Ghost envelope provides confining pressure at the edge
2. **Equilibrium Convergence**: Particles settle to positions where SPH forces match analytical equilibrium
3. **Envelope Penetration Prevention**: Ghost particles at correct pressure stay outside
4. **Clean Starting Point**: Removes transient oscillations before IMBH interaction

### Relaxation Process

1. **GLASS Pre-Relaxation**: Uniformizes particle spacing using repulsive forces
2. **Main Relaxation**:
   - Computes SPH pressure gradient forces
   - Subtracts analytical equilibrium forces
   - Net force drives particles toward equilibrium positions
   - Velocities zeroed each step (quasi-static)
3. **Ghost Removal**: After relaxation, ghost particles are removed

## Expected Evolution

| Phase | Time | R [pc] | n_c [cm^-3] | T [K] | sigma_v [km/s] |
|-------|------|--------|-------------|-------|----------------|
| Relaxed | 0 | ~1-2 | ~200 | 15 | 0 |
| Tidal shaping | 1 Myr | 1-1.5 | ~1000 | 15 | 0.5-1 |
| Pre-pericenter | 2-3 Myr | 0.3-0.5 | ~10^4 | 15 | 1-2 |
| Pericenter | +0.1 Myr | stretching | 10^5-10^6 | 10^4 (shock) | 50-100 |
| Post-encounter | +0.2 Myr | 0.3 | 10^6.5 | **60** | **22** |

## Configuration Files

```
config/
├── presets/
│   ├── option_a_atomic_cnm.json        # Recommended: 40 Msun atomic CNM BE sphere
│   ├── be_sphere_relaxation.json       # Legacy: Phase 1 Relaxation
│   ├── be_sphere_imbh_interaction.json # Legacy: Phase 2 IMBH interaction
│   └── supercritical_be_*.json         # Super-critical configurations
└── active_config.json                  # Symlink to current config
```

### Key Relaxation Parameters

```json
{
    "relaxation": {
        "useRelaxation": true,
        "relaxationSteps": 5000,
        "relaxationOnly": true,
        "useGlassRelaxation": true,
        "glassRelaxationSteps": 2000
    },
    "sample": {
        "useEnvelope": true,
        "envelopeLayers": 6
    }
}
```

## References

1. Oka, T. et al. (2017). "Millimetre-wave Emission from an Intermediate-Mass Black Hole Candidate in the Milky Way." Nature Astronomy, 1, 709-712.

2. Ballone, A. et al. (2018). "Tidal disruption of molecular clouds by intermediate-mass black holes." MNRAS, 480, 4684-4692.

3. Koyama, H. & Inutsuka, S. (2000). "Molecular Cloud Formation in Shock-compressed Layers." ApJ, 532, 980-993.

## Directory Structure

```
imbh_cloud/
├── README.md           # This file
├── Makefile            # Build and run automation
├── config/
│   └── presets/
│       ├── be_sphere_relaxation.json
│       └── be_sphere_imbh_interaction.json
├── results/
│   ├── be_sphere_relaxation/  # Relaxation snapshots
│   └── be_sphere_imbh/        # IMBH interaction snapshots
└── scripts/
    ├── visualize_bonnor_ebert_relaxation.py
    ├── hvcc_visualize.py
    └── quick_viz.py
```
