# 3D Cylindrical Lane-Emden Test Suite

This test suite studies GSPH/SSPH energy conservation and grad-h effects
for a 3D cylindrical Lane-Emden configuration (γ = 5/3, polytropic index n = 1.5).

## Geometry

- **Configuration**: Infinite cylinder along z-axis in 3D
- **Density profile**: Cylindrical Lane-Emden radially in xy-plane: ρ(r) = ρ_c θ^n
- **Uniformity**: Constant in z-direction
- **Gravity**: 2D radial gravity (xy-plane only, like an infinite line mass)

The cylindrical Lane-Emden equation:
```
(1/ξ) d/dξ (ξ dθ/dξ) = -θⁿ   with θ(0) = 1, θ'(0) = 0
```

This is equivalent to:
```
d²θ/dξ² + (1/ξ) dθ/dξ = -θⁿ
```

## Physical Setup

- **γ = 5/3** (adiabatic index)
- **n = 1.5** (polytropic index, n = 1/(γ-1))
- **Wendland C4 kernel** for better energy conservation
- **Kernel-convolved 2D radial gravity** in xy-plane

## Test Cases

| Config | SPH Type | Grad-h | Expected Result |
|--------|----------|--------|-----------------|
| `gsph_gradh.json` | GSPH | Yes | Stable equilibrium |
| `gsph_nogradh.json` | GSPH | No | **Core collapse** (diffusive instability) |
| `ssph_gradh.json` | SSPH | Yes | Stable equilibrium |
| `ssph_nogradh.json` | SSPH | No | **Core collapse** (diffusive instability) |

## Build & Run

### Build for 3D
```bash
cd build
cmake -DSPH_DIM=3 ..
make -j4
```

### Run simulations
```bash
# From project root
./build/sph simulations/stability/gradh_study_3d_cylinder/config/presets/gsph_gradh.json
./build/sph simulations/stability/gradh_study_3d_cylinder/config/presets/gsph_nogradh.json
./build/sph simulations/stability/gradh_study_3d_cylinder/config/presets/ssph_gradh.json
./build/sph simulations/stability/gradh_study_3d_cylinder/config/presets/ssph_nogradh.json
```

## Key Parameters

| Parameter | Value |
|-----------|-------|
| N | 25 (25³ = 15625 particles) |
| R | 1.0 |
| L_z | 1.0 |
| M_total | 1.0 |
| γ | 5/3 |
| kernel | Wendland C4 |
| neighborNumber | 50 |
| useKernelGravityCylinder3D | true |

## Expected Physics

The cylindrical configuration tests a different geometry from the
1D planar slab and 2D disk cases. The cylindrical Lane-Emden has:

- Different surface location ξ₁
- Different mass distribution (∝ ξ θⁿ instead of θⁿ alone)
- 2D logarithmic gravity potential

Without grad-h correction, the same diffusive instability occurs,
leading to artificial core collapse. With grad-h correction, the
equilibrium remains stable.

## Notes

- Particles are placed using random sampling with radii mapped
  through the inverse cumulative mass function
- For better initial conditions, run relaxation first
- The z-coordinate is uniform random within [-L_z/2, +L_z/2]
