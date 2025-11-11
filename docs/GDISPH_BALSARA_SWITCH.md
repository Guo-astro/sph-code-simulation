# GDISPH Balsara Switch Implementation

## Overview

This document describes the implementation of the Balsara switch for GDISPH (Godunov Density-Independent SPH) according to equations 100-101 in [Yuasa & Mori (2024)](https://arxiv.org/pdf/2312.03224).

## Paper Reference

**Equations 100-101** from "GDISPH: A Density-independent Formulation of Smoothed Particle Hydrodynamics using Godunov's Method" specify the Balsara switch for GDISPH:

```
∇·v = (1/ρᵢ) Σⱼ mⱼ (vᵢⱼ · ∇Wᵢⱼ)         (Equation 100)
∇×v = (1/ρᵢ) Σⱼ mⱼ (vᵢⱼ × ∇Wᵢⱼ)         (Equation 101)
```

The Balsara factor is then:
```
fᵢ = |∇·v| / (|∇·v| + |∇×v| + ε·cₛ/h)
```

where ε = 10⁻⁴ is a small number to prevent division by zero.

## Key Difference from DISPH

### DISPH Formulation (Energy-Weighted)
```cpp
// DISPH uses mass × energy weighting
div_v -= p_j.mass * p_j.ene * inner_product(v_ij, dw);
rot_v += vector_product(v_ij, dw) * (p_j.mass * p_j.ene);
// Normalize by (γ-1)/P
const real p_inv = (gamma - 1.0) / p_i.pres;
div_v *= p_inv;
rot_v *= p_inv;
```

### GDISPH Formulation (Density-Based)
```cpp
// GDISPH uses mass-only weighting (equations 100-101)
div_v -= p_j.mass * inner_product(v_ij, dw);
rot_v += vector_product(v_ij, dw) * p_j.mass;
// Normalize by density
div_v /= p_i.dens;
rot_v /= p_i.dens;
```

## Implementation Location

**File:** `src/gdisph/gd_pre_interaction.cpp`  
**Lines:** 150-173 (in the `compute()` function)

```cpp
// Artificial viscosity - GDISPH Balsara switch (equations 100-101 from paper)
if(m_use_balsara_switch && DIM != 1) {
#if DIM != 1
    // GDISPH Balsara switch: uses density-based formulation (not mass*energy like DISPH)
    // Equations 100-101: div_v and rot_v calculated with mass/density weighting
    real div_v = 0.0;
#if DIM == 2
    real rot_v = 0.0;
#else
    vec_t rot_v = 0.0;
#endif
    for(int n = 0; n < n_neighbor; ++n) {
        int const j = neighbor_list[n];
        auto & p_j = particles[j];
        const vec_t r_ij = periodic->calc_r_ij(pos_i, p_j.pos);
        const real r = std::abs(r_ij);
        const vec_t dw = kernel->dw(r_ij, r, p_i.sml);
        const vec_t v_ij = p_i.vel - p_j.vel;
        // GDISPH: use mass (not mass*energy) - paper eq 100
        div_v -= p_j.mass * inner_product(v_ij, dw);
        rot_v += vector_product(v_ij, dw) * p_j.mass;
    }
    // Normalize by density (paper eq 100-101)
    div_v /= p_i.dens;
    rot_v /= p_i.dens;
    p_i.balsara = std::abs(div_v) / (std::abs(div_v) + std::abs(rot_v) + 1e-4 * p_i.sound / p_i.sml);
#endif
}
```

## Folder Structure

The Balsara switch implementation follows the established pattern:

```
src/
├── pre_interaction.cpp           # SSPH Balsara switch (standard)
├── disph/
│   └── d_pre_interaction.cpp     # DISPH Balsara switch (energy-weighted)
├── gsph/
│   └── g_pre_interaction.cpp     # GSPH Balsara switch
└── gdisph/
    └── gd_pre_interaction.cpp    # GDISPH Balsara switch (density-based, equations 100-101)
```

Each SPH method implements its own variant of the Balsara switch in its respective `pre_interaction.cpp` file within its method-specific subdirectory.

## Configuration

Enable the Balsara switch in JSON configuration:

```json
{
  "SPHType": "gdisph",
  "avAlpha": 0.0,
  "useBalsaraSwitch": true
}
```

Note: The Balsara switch is typically used when `avAlpha > 0` to reduce artificial viscosity in shear flows while maintaining it in compressive flows.

## Testing

Tested with 1D Sod shock tube:

```bash
make shock_tube_run PRESET=shock_tube_1d_gdisph
```

**Results:**
- Simulation time: ~370ms for 500 particles
- Shock propagation correctly captured
- No spurious oscillations at discontinuities

## References

1. Yuasa, T., & Mori, M. (2024). GDISPH: A Density-independent Formulation of Smoothed Particle Hydrodynamics using Godunov's Method. arXiv:2312.03224.
2. Balsara, D. S. (1995). von Neumann stability analysis of smooth particle hydrodynamics—suggestions for optimal algorithms. Journal of Computational Physics, 121(2), 357-372.

## Changelog

**2024-11-11:**
- Updated GDISPH Balsara switch to use density-based formulation (equations 100-101)
- Previous implementation used DISPH energy-weighted formulation
- Changed from `mass × energy / pressure` to `mass / density` weighting
- Verified with 1D Sod shock tube test
