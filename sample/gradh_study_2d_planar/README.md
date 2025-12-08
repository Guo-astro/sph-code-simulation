# 2D Planar Lane-Emden Slab Test Suite

This test suite studies GSPH/SSPH energy conservation and grad-h effects
for a 2D planar Lane-Emden slab (γ = 5/3, polytropic index n = 1.5).

## Geometry

- **Configuration**: Infinite slab in 2D
- **Density profile**: Planar Lane-Emden in y-direction: ρ(y) = ρ_c θ^n
- **Uniformity**: Constant in x-direction
- **Gravity**: 1D planar gravity (y-direction only)

The planar Lane-Emden equation:
```
d²θ/dξ² = -θⁿ   with θ(0) = 1, θ'(0) = 0
```

## Physical Setup

- **γ = 5/3** (adiabatic index)
- **n = 1.5** (polytropic index, n = 1/(γ-1))
- **Wendland C4 kernel** for better energy conservation
- **Kernel-convolved 1D gravity** in y-direction

## Test Cases

| Config | SPH Type | Grad-h | Expected Result |
|--------|----------|--------|-----------------|
| `gsph_gradh.json` | GSPH | Yes | Stable equilibrium |
| `gsph_nogradh.json` | GSPH | No | **Core collapse** (diffusive instability) |
| `ssph_gradh.json` | SSPH | Yes | Stable equilibrium |
| `ssph_nogradh.json` | SSPH | No | **Core collapse** (diffusive instability) |

## Build & Run

### Build for 2D
```bash
cd build
cmake -DSPH_DIM=2 ..
make -j4
```

### Run simulations
```bash
# From project root
./build/sph sample/gradh_study_2d_planar/config/presets/gsph_gradh.json
./build/sph sample/gradh_study_2d_planar/config/presets/gsph_nogradh.json
./build/sph sample/gradh_study_2d_planar/config/presets/ssph_gradh.json
./build/sph sample/gradh_study_2d_planar/config/presets/ssph_nogradh.json
```

## Key Parameters

| Parameter | Value |
|-----------|-------|
| N | 50 (50² = 2500 particles) |
| ρ_center | 1.0 |
| K | 1.0 |
| L_x | 2.0 |
| γ | 5/3 |
| kernel | Wendland C4 |
| neighborNumber | 32 |
| useKernelGravityPlanar2D | true |

## Expected Physics

Without grad-h correction, the GSPH/SSPH methods exhibit a diffusive
instability that leads to artificial core collapse. The growth rate
is approximately:

```
Γ ≈ ε × c_s × h / L²
```

where ε ~ 0.4 is the force error coefficient, c_s is the sound speed,
h is the smoothing length, and L is the characteristic scale.

With grad-h correction, this instability is suppressed and the
equilibrium remains stable.
