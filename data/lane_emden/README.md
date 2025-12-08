# Lane-Emden Data Files

This directory contains pre-computed solutions to the Lane-Emden equation for
polytropic indices relevant to astrophysical simulations.

## Lane-Emden Equation

For a polytrope of index n, the Lane-Emden equation in 3D is:

```
(1/ξ²) d/dξ (ξ² dθ/dξ) = -θⁿ
```

with boundary conditions θ(0) = 1, dθ/dξ(0) = 0.

## File Naming Convention

Files are named as `n{index}_{dim}d.dat` where:
- `{index}` is the polytropic index (e.g., 1.5)
- `{dim}` is the dimension (2 or 3)

## Canonical Files (USE THESE)

| File | Index n | Dimension | First Zero ξ₁ | |dθ/dξ|₁ |
|------|---------|-----------|---------------|----------|
| `n1.5_3d.dat` | 1.5 | 3D | 3.6537540101 | 0.20331244769 |
| `n1.5_2d.dat` | 1.5 | 2D | (varies) | (varies) |

## Key Physical Relations

For n=1.5 (γ=5/3 polytrope):

- **Density**: ρ(ξ) = ρ_c θ^(3/2)
- **Pressure**: P(ξ) = K ρ_c^γ θ^(n+1)
- **Scaling**: r = α ξ, where α = R / ξ₁
- **Mass**: M = 4π α³ ρ_c ξ₁² |dθ/dξ|₁

## File Format

```
# Lane-Emden solution for 3D geometry, n=1.5
# xi_1 = 3.6537540101e+00
# dtheta_1 = -2.0331244769e-01
# n_points = 366
# Columns: xi  theta  dtheta/dxi
1.0e-10  1.0  0.0
...
```

## IMPORTANT NOTES

1. Always validate against known values: ξ₁ ≈ 3.6538 for n=1.5 in 3D
2. Use the shared Python module (`scripts/shared/lane_emden.py`) for SSOT

## Python Usage

Use the shared module for SSOT:

```python
from scripts.shared.lane_emden import load_lane_emden_solution, get_density_profile

# Load solution
solution = load_lane_emden_solution(n=1.5, dim=3)

# Get physical density profile
r, rho = get_density_profile(solution, rho_c=1.0, R=1.0)
```

## C++ Usage

The C++ code loads these files via `LaneEmdenData` class in:
- `include/relaxation/lane_emden_data.hpp`
- `src/relaxation/lane_emden_data.cpp`

## References

- Chandrasekhar, S. (1939). "An Introduction to the Study of Stellar Structure"
- Horedt, G. P. (2004). "Polytropes: Applications in Astrophysics and Related Fields"
