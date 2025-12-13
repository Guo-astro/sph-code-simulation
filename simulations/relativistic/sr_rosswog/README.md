# Rosswog (2010) SR-SPH Benchmark Tests

This directory contains benchmark tests from:

**Rosswog, S. (2010)** "Conservative, special-relativistic smoothed particle hydrodynamics"  
arXiv:0907.4890v3, Journal of Computational Physics

## Test Cases

### 1D Tests

| Test | Name | Description | γ_eos | Key Challenge |
|------|------|-------------|-------|---------------|
| 1 | Standard SR Shock Tube | Marti & Mueller benchmark | 5/3 | Mildly relativistic (γ_max≈1.4) |
| 2 | Strong Blast | Norman blast wave | 5/3 | Thin dense shell (γ_shell=3.6, γ_shock=6.0) |
| 3 | Perturbed Shock Tube | Sinusoidal density perturbation | 5/3 | Smooth structure transport |
| 4 | Wall Shock | Ultra-relativistic wall reflection | 4/3 | Extreme Lorentz factor (γ=50,000) |
| 5 | Sine Wave Advection | Smooth advection test | 4/3 | Conservation (γ≈12.92) |
| 6 | Square Wave Advection | Steep gradient advection | 4/3 | Sharp feature preservation |
| 7 | Simple Wave | Relativistic simple wave evolution | 4/3 | Shock formation & dissipation |

### 2D Tests

| Test | Name | Description |
|------|------|-------------|
| 8 | 2D Shock Tube | Multi-dimensional version of Test 1 |

## Initial Conditions

### Test 1: Standard Relativistic Shock Tube
```
Left:  (N, v, P) = (10, 0, 40/3)
Right: (N, v, P) = (1, 0, 10⁻⁶)
Domain: [-0.5, 0.5], t_end = 0.35
```

### Test 2: Strong Blast
```
Left:  (N, v, P) = (1, 0, 1000)
Right: (N, v, P) = (1, 0, 0.01)
Domain: [-0.5, 0.5], t_end = 0.16
```

### Test 3: Perturbed Shock Tube
```
Left:  (N, v, P) = (5, 0, 50)
Right: (N, v, P) = (2 + 0.3*sin(50x), 0, 5)
Domain: [0, 1], t_end = 0.35
```

### Test 4: Ultra-Relativistic Wall Shock
```
Initial: (N, v, u) = (1, 0.9999999998, 10⁻⁵)
Wall at x = 1 (reflecting boundary)
γ = 50,000, Domain: [0, 1], t_end = 1.0
Post-shock: N_p = 4*N_i, P_p = γ*Γ*N_i, u_p = γ
```

### Test 5: Sine Wave Advection
```
N(x) = 1 + 0.5*sin(2πx), v = 0.997 (γ≈12.92)
Constant pressure P_0 = (Γ-1)*n_0*u_0
Periodic box [0, 1], 100 box crossings
```

### Test 6: Square Wave Advection
```
N = N_low + (N_high - N_low)*Fermi_step
N_low = 1.0, N_high = 1.1, v = 0.997
Periodic box, 10 box crossings
```

### Test 7: Simple Wave
```
v_max = 0.7, c_s,0 = 0.3, n_0 = 1
Radiation-dominated: P = k(s)*ρ^(4/3)
Pulse width l_0 = 100 (characteristic length)
```

## Usage

```bash
# From repository root:
make -f simulations/relativistic/sr_rosswog/Makefile.sr_rosswog sr_rosswog_help

# Run individual tests
make -f simulations/relativistic/sr_rosswog/Makefile.sr_rosswog sr_rosswog_test1_run
make -f simulations/relativistic/sr_rosswog/Makefile.sr_rosswog sr_rosswog_test2_run

# Run all 1D tests
make -f simulations/relativistic/sr_rosswog/Makefile.sr_rosswog sr_rosswog_1d_all
```

## References

1. Rosswog, S. (2010), JCP 229, 8591-8612, arXiv:0907.4890
2. Marti, J.M. & Müller, E. (2003), Living Rev. Relativ. 6, 7
3. Norman, M.L. & Winkler, K.-H.A. (1986), NATO ASI 188, 187
4. Anile, A.M. (1989), Relativistic Fluids and Magneto-fluids
