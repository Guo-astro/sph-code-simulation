# SR-GSPH Tangential Velocity Tests

Special relativistic Riemann problem tests with non-zero tangential velocity component.

## References

1. **Pons et al. (2000)** - "The exact solution of the Riemann problem with non-zero tangential velocities in relativistic hydrodynamics"
   - DOI: [10.1017/S0022112000001439](https://doi.org/10.1017/S0022112000001439)
   - Table 1: Test cases with tangential velocities v^t_L, v^t_R ∈ {0, 0.9, 0.99}

2. **Rezzolla et al. (2003)** - "An improved exact Riemann solver for relativistic hydrodynamics"
   - DOI: [10.1017/S0022112002003506](https://doi.org/10.1017/S0022112002003506)
   - Table 1: Test cases with various wave patterns (SR, 2S, 2R)

## Test Cases

### Pons2000 Tests (Strong Pressure Jump)

Initial conditions:
- Left:  (P, ρ, v^x) = (1000, 1, 0)
- Right: (P, ρ, v^x) = (0.01, 1, 0)
- γ = 5/3

| Test | v^t_L | v^t_R | P* | v^x* | Wave Pattern |
|------|-------|-------|-----|------|--------------|
| P00  | 0.00  | 0.00  | 18.6 | 0.960 | R-S |
| P01  | 0.00  | 0.90  | 42.8 | 0.913 | R-S |
| P02  | 0.00  | 0.99  | 127  | 0.767 | R-S |
| P10  | 0.90  | 0.00  | 0.189| 0.328 | R-S |
| P11  | 0.90  | 0.90  | 0.904| 0.319 | R-S |
| P12  | 0.90  | 0.99  | 8.48 | 0.292 | R-S |
| P20  | 0.99  | 0.00  | 0.0316| 0.099 | R-S |
| P21  | 0.99  | 0.90  | 0.0927| 0.098 | R-S |
| P22  | 0.99  | 0.99  | 0.706| 0.095 | R-S |

### Rezzolla2003 Tests (Moderate Pressure Jump)

Initial conditions:
- Left:  (P, ρ) = (1.0, 1.0)
- Right: (P, ρ) = (0.1, 0.125)
- γ = 5/3

| Test | v^x_L | v^x_R | v^t_L | v^t_R | P* | v^x* | Pattern |
|------|-------|-------|-------|-------|-----|------|---------|
| R01  | 0.5   | 0.0   | 0.0   | 0.0   | 0.597 | 0.640 | SR |
| R02  | 0.5   | 0.0   | 0.0   | 0.3   | 0.621 | 0.631 | SR |
| R03  | 0.5   | 0.0   | 0.0   | 0.5   | 0.673 | 0.611 | SR |
| R04  | 0.5   | 0.0   | 0.0   | 0.7   | 0.787 | 0.570 | SR |
| R05  | 0.5   | 0.0   | 0.0   | 0.9   | 1.150 | 0.455 | 2S |
| R06  | 0.5   | 0.0   | 0.0   | 0.99  | 2.199 | 0.212 | 2S |
| R07  | 0.5   | 0.0   | 0.0   | 0.999 | 3.011 | 0.078 | 2S |
| R08  | 0.0   | 0.5   | 0.0   | 0.0   | 0.154 | 0.620 | SR |
| R09  | 0.0   | 0.5   | 0.3   | 0.0   | 0.139 | 0.594 | SR |
| R10  | 0.0   | 0.5   | 0.5   | 0.0   | 0.115 | 0.542 | SR |
| R11  | 0.0   | 0.5   | 0.7   | 0.0   | 0.085 | 0.450 | 2R |
| R12  | 0.0   | 0.5   | 0.9   | 0.0   | 0.051 | 0.280 | 2R |
| R13  | 0.0   | 0.5   | 0.99  | 0.0   | 0.031 | 0.095 | 2R |
| R14  | 0.0   | 0.5   | 0.999 | 0.0   | 0.026 | 0.031 | 2R |

## Usage

```bash
# Run all tangent velocity tests
make sr_tangent_test_all

# Run individual test suites
make sr_tangent_test_pons2000
make sr_tangent_test_rezzolla2003

# Clean results
make sr_tangent_clean
```

## Physical Background

In special relativistic hydrodynamics, the tangential velocity component `v^t` introduces additional coupling through the Lorentz factor:

```
W = 1 / √(1 - v^x² - v^t²)
```

The conserved quantity across rarefaction waves is:
```
A = h · W · v^t = const
```

This means the tangential velocity `v^t` is NOT simply advected but transforms as:
```
v^t_b = A · √((1 - v^x_b²) / (h_b² + A²))
```

## Implementation Notes

The exact Riemann solver uses the arctanh transformation for rarefaction waves:
```
B = arctanh(v^x_a) + sign · ∫[P_a to P_b] f(P) dP
v^x_b = tanh(B)
```

where the integrand is:
```
f(P) = √(h² + A²(1 - c_s²)) / ((h² + A²) · ρ · c_s)
```

This matches the Python reference implementation in `docs/papers/sg-gsph/srrp/`.
