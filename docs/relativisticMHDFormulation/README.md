# Special Relativistic Magnetohydrodynamics in Godunov SPH

This directory contains the mathematical formulation and implementation guide for Special Relativistic Magnetohydrodynamics (SRMHD) using the Godunov Smoothed Particle Hydrodynamics (GSPH) framework.

## Overview

This formulation combines two key methodologies:

1. **SR-GSPH** (Kitajima et al. 2025): Special relativistic hydrodynamics with Riemann solver
2. **GSPMHD** (Iwasaki & Inutsuka 2011): MHD extension using Method of Characteristics

## Documents

| File | Description |
|------|-------------|
| `SRMHD_GSPH_Formulation.tex` | Complete LaTeX document with step-by-step derivations |
| `README.md` | This overview and implementation guide |

## Quick Reference

### Fundamental Variables

| Symbol | Name | Definition |
|--------|------|------------|
| `n` | Rest-frame density | Baryon number per proper volume |
| `N` | Lab-frame density | N = gamma * n |
| `gamma` | Lorentz factor | 1/sqrt(1 - v^2/c^2) |
| `H` | Specific enthalpy | 1 + u/c^2 + P/(n*c^2) |
| `S` | Canonical momentum | S = gamma * H * v |
| `e` | Canonical energy | e = gamma*H - P/(N*c^2) |
| `B` | Magnetic field | Lab-frame magnetic field |

### Key Equations

#### 1. Momentum Equation (SPH Form)

```
dS_i/dt = -sum_j m_j * F_ij * [
    -(P_t^* - P_B_parallel) * n_ij   // Compressive (Riemann)
    + B_parallel^* * B_perp^*         // Alfvenic (MOC)
]
```

#### 2. Energy Equation (SPH Form)

```
de_i/dt = -sum_j m_j * F_ij * [
    -(P_t^* - P_B_parallel) * v_parallel^*   // Compressive flux
    + B_parallel^* * (B_perp^* . v_perp^*)   // Alfvenic flux
]
```

#### 3. Induction Equation (SPH Form)

```
d(B/N)_i/dt = sum_j m_j * F_ij * B_parallel^* * [
    v_parallel^* * n_ij + v_perp^* - v_i^*
]
```

## Algorithm Overview

### Step 1: Pre-Interaction

1. Compute smoothing length `h_i` iteratively
2. Compute particle volume `V_p = 1 / sum_j W(r_ij, h)`
3. Compute lab-frame density `N = nu / V_p`
4. Recover primitive variables from conserved

### Step 2: Force Calculation

For each particle pair (i, j):

```python
# Compute separation
r_ij = x_i - x_j
r = |r_ij|
n_ij = r_ij / r  # Unit vector

# Decompose velocities and B-fields
v_parallel_i = dot(v_i, n_ij)
v_perp_i = v_i - v_parallel_i * n_ij
B_parallel_i = dot(B_i, n_ij)
B_perp_i = B_i - B_parallel_i * n_ij

# Total pressure
P_t_i = P_i + 0.5 * |B_perp_i|^2

# --- Riemann Problem (Compressive) ---
# Build left/right states
state_L = (rho_L, P_t_L, v_parallel_L, c_fast_L)
state_R = (rho_R, P_t_R, v_parallel_R, c_fast_R)

# Solve Riemann problem
P_t_star, v_parallel_star = riemann_solver(state_L, state_R)

# --- Method of Characteristics (Alfvenic) ---
B_perp_star, v_perp_star = moc_alfven(
    v_perp_L, B_perp_L, rho_L,
    v_perp_R, B_perp_R, rho_R,
    B_parallel_avg
)

# --- Compute forces ---
P_B_parallel = 0.25 * (B_parallel_i^2 + B_parallel_j^2)
B_parallel_star = 0.5 * (B_parallel_i + B_parallel_j)

# Kernel gradients
grad_W_i = dW(r_ij, sqrt(2)*h_i)
grad_W_j = dW(r_ij, sqrt(2)*h_j)
F_ij = V_i^2 * grad_W_i + V_j^2 * grad_W_j

# Accumulate momentum
force = -(P_t_star - P_B_parallel) * n_ij
        + B_parallel_star * B_perp_star
dS_i += F_ij * force

# Accumulate energy
v_star = v_parallel_star * n_ij + v_perp_star
de_i += dot(F_ij * force, v_star)

# Accumulate induction
dB_i += F_ij * B_parallel_star * (v_star - v_i_centered)
```

### Step 3: Time Integration

```python
# Update conserved variables
S_i = S_i + dS_i * dt
e_i = e_i + de_i * dt
x_i = x_i + v_i * dt
(B/N)_i = (B/N)_i + dB_i * dt

# CFL condition
dt = C_CFL * min(h_i / c_fast_i)
```

## Riemann Solver Strategy

### Decomposition

The MHD Riemann problem is split into:

1. **Compressive part** (fast magnetosonic): Solved with relativistic Riemann solver
2. **Alfvenic part** (Alfven waves): Solved with Method of Characteristics

### Relativistic Fast Riemann Solver

Input: Left/Right states with total pressure P_t = P + B_perp^2/2

For shock waves:
```
H_b = solve_quadratic(H_a, P_a, P_b, gamma)
j^2 = -(P_b - P_a) / (H_b/n_b - H_a/n_a)
V_s = shock_speed(j, N_a, v_a)
v_b = post_shock_velocity(H_a, gamma_a, P_a, P_b, j, V_s)
```

For rarefaction waves:
```
v_b = tanh(arctanh(v_a) + sign * integral(dP / (rho * h * c_s)))
```

### Method of Characteristics for Alfven Waves

```python
def moc_alfven(v_perp_L, B_perp_L, rho_L, v_perp_R, B_perp_R, rho_R, B_parallel):
    sqrt_rho_L = sqrt(rho_L * H_L * gamma_L^2)  # Relativistic impedance
    sqrt_rho_R = sqrt(rho_R * H_R * gamma_R^2)
    sign_B = sign(B_parallel)

    inv_L = 1 / sqrt_rho_L
    inv_R = 1 / sqrt_rho_R

    B_perp_star = (B_perp_L * inv_L + B_perp_R * inv_R
                   + sign_B * (v_perp_R - v_perp_L)) / (inv_L + inv_R)

    v_perp_star = (v_perp_L * sqrt_rho_L + v_perp_R * sqrt_rho_R
                   + sign_B * (B_perp_R - B_perp_L)) / (sqrt_rho_L + sqrt_rho_R)

    return B_perp_star, v_perp_star
```

## Primitive Variable Recovery

From conserved (S, e, N, B) to primitive (v, n, P, gamma):

```python
def recover_primitives(S, e, N, B, gamma_eos):
    X = gamma_eos / (gamma_eos - 1)
    S_mag = |S|

    # Solve quartic for Lorentz factor
    # (gamma^2 - 1)(X*e*gamma - 1)^2 - S^2*(X*gamma^2 - 1)^2 = 0
    gamma = newton_raphson_quartic(S_mag, e, X)

    # Recover velocity
    denom = gamma * (X * e * gamma - 1)
    v = S * (X * gamma^2 - 1) / denom

    # Recover other primitives
    n = N / gamma
    H = (X * e * gamma - 1) / (X * gamma^2 - 1)
    P = n * (H - 1) * (gamma_eos - 1) / gamma_eos
    c_s = sqrt((gamma_eos - 1) * (H - 1) / H)

    return v, n, P, gamma, H, c_s
```

## Divergence Cleaning (Powell Terms)

To maintain numerical stability and div(B) = 0:

```python
# Modified momentum equation
dS_i -= B_i * sum_j(m_j * B_parallel_star * F_ij)

# Modified energy equation
de_i -= dot(B_i, v_i_star) * sum_j(m_j * B_parallel_star * F_ij)
```

Note: Powell terms break exact conservation but prevent tensile instability.

## Test Problems

### 1D Tests
- Balsara shock tubes (1-5)
- Relativistic MHD blast waves
- Alfven wave propagation

### 2D Tests
- Relativistic Orszag-Tang vortex
- Relativistic rotor problem
- Kelvin-Helmholtz instability with magnetic fields

## References

1. Kitajima, K., Inutsuka, S., & Seno, I. (2025). Special Relativistic Smoothed Particle Hydrodynamics Based on Riemann Solver. arXiv:2510.18251.

2. Iwasaki, K., & Inutsuka, S. (2011). Smoothed particle magnetohydrodynamics with a Riemann solver and the method of characteristics. MNRAS, 418, 1668.

3. Inutsuka, S. (2002). Reformulation of Smoothed Particle Hydrodynamics with Riemann Solver. J. Comput. Phys., 179, 238.

4. Pons, J. A., Marti, J. M., & Muller, E. (2000). The exact solution of the Riemann problem with non-zero tangential velocities in relativistic hydrodynamics. J. Fluid Mech., 422, 125.

5. Mignone, A., & Bodo, G. (2005). An HLLC Riemann solver for relativistic flows. MNRAS, 364, 126.
