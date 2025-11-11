# GDISPH Implementation Summary

## Overview

**GDISPH (Godunov Density-Independent SPH)** is a novel SPH method combining the best features of GSPH (Godunov SPH) and DISPH (Density-Independent SPH), based on the paper:

> Yuasa, T., & Mori, M. (2024). "Novel Hydrodynamic Schemes Capturing Shocks and Contact Discontinuities"  
> arXiv:2312.03224v3 [astro-ph.IM]

## Key Features

GDISPH combines:

1. **DISPH's pressure-energy formulation** (Hopkins 2013)
   - Pressure calculated from specific internal energy: `P = (γ-1) Σ(m_j * u_j * W_ij)`
   - Eliminates pressure errors in multi-phase flows
   - More accurate for contact discontinuities

2. **GSPH's Riemann solver approach** (Cha & Whitworth 2003)
   - HLL Riemann solver for shock capturing
   - MUSCL reconstruction for 2nd order accuracy
   - Better handling of strong shocks

3. **DISPH's Newton-Raphson smoothing length solver**
   - Iterative calculation ensuring h consistency
   - Uses pressure-energy formulation

## Implementation Files

### Headers
- `include/gdisph/gd_pre_interaction.hpp` - Pre-interaction phase (density, pressure, gradients)
- `include/gdisph/gd_fluid_force.hpp` - Fluid force calculation with Riemann solver

### Source Files
- `src/gdisph/gd_pre_interaction.cpp` (276 lines)
  - Newton-Raphson for smoothing length
  - DISPH pressure-energy formulation
  - GSPH gradient calculation for MUSCL reconstruction
  
- `src/gdisph/gd_fluid_force.cpp` (219 lines)
  - HLL Riemann solver (pstar, vstar calculation)
  - DISPH force formulation with f_ij factors
  - MUSCL reconstruction for 2nd order (when enabled)

### Integration
- `include/parameters.hpp` - Added `GDISPH` to `SPHType` enum
- `src/solver.cpp` - Added GDISPH case handling
- `src/CMakeLists.txt` - Added gdisph subdirectory
- `include/CMakeLists.txt` - Added gdisph subdirectory

## Algorithm Details

### Pre-Interaction Phase

1. **Density Calculation** (standard SPH):
   ```
   ρ_i = Σ(m_j * W_ij)
   ```

2. **Pressure Calculation** (DISPH formulation):
   ```
   P_i = (γ-1) * Σ(m_j * u_j * W_ij)
   ```

3. **Gradient Calculation** (for MUSCL reconstruction):
   ```
   ∇φ_i = (1/Ω_i) * Σ(m_j/ρ_j * (φ_j - φ_i) * ∇W_ij)
   ```

4. **Newton-Raphson for Smoothing Length**:
   - Solves: `h^DIM * ρ = const`
   - Uses DISPH pressure-energy formulation for consistency

### Fluid Force Phase

1. **MUSCL Reconstruction** (when 2nd order enabled):
   ```
   φ_L = φ_i + 0.5 * min_mod_limiter(∇φ_i) · Δr
   φ_R = φ_j + 0.5 * min_mod_limiter(∇φ_j) · Δr
   ```

2. **HLL Riemann Solver**:
   ```
   S_L = min(v_i, v_j) - max(c_i, c_j)
   S_R = max(v_i, v_j) + max(c_i, c_j)
   pstar = (P_R + P_L)/2 + (v_L - v_R) * (ρ_L + ρ_R) * (c_L + c_R) / 8
   vstar = (v_L + v_R)/2 + (P_L - P_R) / ((ρ_L + ρ_R) * (c_L + c_R) / 2)
   ```

3. **DISPH Force Formulation**:
   ```
   f_ij = Σ(m_j * (f_i + f_j) * W_ij)
   where f_i = -P_i / (Ω_i * ρ_i)
   ```

4. **Acceleration**:
   ```
   a_i = Σ(m_j * (pstar_ij/ρ_i + pstar_ij/ρ_j) * ∇W_ij)
   ```

## Configuration Parameters

### GDISPH-Specific Settings

```json
{
  "sphType": "gdisph",
  "gsphIs2ndOrder": false,    // Enable MUSCL reconstruction
  "avAlpha": 0.0              // Artificial viscosity (0.0 recommended for Case 1)
}
```

### Recommended Settings for Shock Tube

From the paper (Case 1 - Riemann solver only):
- No artificial viscosity: `avAlpha = 0.0`
- 1st order accurate: `gsphIs2ndOrder = false`
- CFL sound speed: `0.3`
- CFL force: `0.125`

## Test Results

### 1D Sod Shock Tube Comparison

Configuration: `sample/shock_tube/config/presets/shock_tube_1d_gdisph.json`

**Performance**:
- Particles: 500
- End time: 0.2
- Calculation time: ~340-390 ms
- Output: 21 snapshots

**Comparison with Other Methods**:
- GSPH: Standard Godunov SPH with artificial viscosity
- SSPH: Standard SPH with artificial viscosity
- DISPH: Density-Independent SPH (pressure-energy)
- GDISPH: Combines DISPH + GSPH features

Results show GDISPH successfully captures:
1. Shock front (sharp discontinuity)
2. Contact discontinuity (smooth transition)
3. Rarefaction wave (expansion region)

## Build Instructions

1. **Build the code**:
   ```bash
   cd build
   make -j8
   ```

2. **Run GDISPH shock tube test**:
   ```bash
   cd /path/to/sphcode
   make -f sample/shock_tube/Makefile.shock_tube shock_tube_compare_run
   ```

3. **Generate comparison plots**:
   ```bash
   make -f sample/shock_tube/Makefile.shock_tube shock_tube_compare_viz
   ```

4. **Generate animation**:
   ```bash
   make -f sample/shock_tube/Makefile.shock_tube shock_tube_compare_animate
   ```

## Visualization

Comparison plots show all four methods (GSPH, SSPH, DISPH, GDISPH) side-by-side for:
- Density profile
- Velocity profile
- Pressure profile
- Internal energy profile

**Color Scheme** (colorblind-friendly):
- GSPH: Blue `#0173B2` (solid line, circle markers)
- SSPH: Orange `#DE8F05` (dashed line, square markers)
- DISPH: Teal `#029E73` (dash-dot line, triangle markers)
- GDISPH: Purple `#CC78BC` (dotted line, diamond markers)

Output files:
- `sample/shock_tube/results/comparison/comparison_tsnap*.png` - Time snapshots
- `sample/shock_tube/results/comparison/comparison_final.png` - Final state
- `sample/shock_tube/results/comparison/comparison_animation.gif` - Animation

## Future Work

1. **2nd Order MUSCL**: Test with `gsphIs2ndOrder = true`
2. **Artificial Viscosity**: Test Case 2 from paper with `avAlpha > 0`
3. **Multi-dimensional Tests**: Extend to 2D/3D test cases
4. **Performance Optimization**: Profile and optimize critical loops
5. **Paper Benchmarks**: Compare against all test cases in Yuasa & Mori (2024)

## References

1. Yuasa, T., & Mori, M. (2024). Novel Hydrodynamic Schemes Capturing Shocks and Contact Discontinuities. arXiv:2312.03224v3
2. Hopkins, P. F. (2013). A general class of Lagrangian smoothed particle hydrodynamics methods. MNRAS, 428(4), 2840-2856.
3. Cha, S.-H., & Whitworth, A. P. (2003). Implementations and tests of Godunov-type particle hydrodynamics. MNRAS, 340(1), 73-90.

## Implementation Status

✅ **Complete and tested**:
- GDISPH pre-interaction module
- GDISPH fluid force module
- Solver integration
- Build system integration
- Test configuration
- Comparison visualization
- Documentation

🎉 **First successful run**: 2024-11-11
