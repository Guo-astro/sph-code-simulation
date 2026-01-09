# SRMHD-GSPH Equation Summary

Quick reference for the final discretized equations.

## Notation

| Symbol | Description |
|--------|-------------|
| i, j | Particle indices |
| n_ij | Unit vector from j to i |
| r_ij | Position vector from j to i |
| V_i | Particle volume = nu_i / N_i |
| h_i | Smoothing length |
| W | Kernel function (Gaussian) |
| Omega_i | Grad-h correction factor |

## Kernel Function

Gaussian kernel:
```
W(r, h) = (1 / (pi^{d/2} h^d)) * exp(-r^2/h^2)
```

Gradient:
```
nabla W(r, h) = -2r/(h^2) * W(r, h)
```

## Density Estimation

Particle volume:
```
V_p(x_i) = [sum_j W(x_i - x_j, h_i)]^{-1}
```

Lab-frame density:
```
N_i = nu_i / V_p(x_i)
```

Grad-h correction:
```
Omega_i = 1 / [1 + (h_i / d) * (sum_j dW/dh) / (sum_j W)]
```

## Smoothing Length

Iterative solution:
```
h = eta * [V_p^*]^{1/d}

V_p^* = [sum_j W(r_ij, C_smooth * h)]^{-1}
```

Typical: eta = 1.0, C_smooth = 2.0

## Force Calculation

### Interface Decomposition

Normal/tangent decomposition:
```
v_parallel = dot(v, n_ij)
v_perp = v - v_parallel * n_ij
B_parallel = dot(B, n_ij)
B_perp = B - B_parallel * n_ij
```

Total pressure:
```
P_t = P + |B_perp|^2 / 2
```

### Riemann Solver (Compressive Part)

States:
```
State_L = (rho_L, P_t_L, v_parallel_L, c_fast_L)
State_R = (rho_R, P_t_R, v_parallel_R, c_fast_R)
```

Output:
```
P_t^*, v_parallel^* = Riemann(State_L, State_R)
```

### Method of Characteristics (Alfvenic Part)

```
Z_L = rho_L * H_L * gamma_L^2    (relativistic impedance)
Z_R = rho_R * H_R * gamma_R^2

sign_B = sign(B_parallel_avg)
B_parallel_avg = (B_parallel_i + B_parallel_j) / 2

B_perp^* = [B_perp_L/sqrt(Z_L) + B_perp_R/sqrt(Z_R)
            + sign_B * (v_perp_R - v_perp_L)]
           / [1/sqrt(Z_L) + 1/sqrt(Z_R)]

v_perp^* = [v_perp_L*sqrt(Z_L) + v_perp_R*sqrt(Z_R)
            + sign_B * (B_perp_R - B_perp_L)]
           / [sqrt(Z_L) + sqrt(Z_R)]
```

### Kernel Gradient Factor

```
F_ij = V_i^2 * Omega_i * nabla_i W(r_ij, sqrt{2}*h_i)
     + V_j^2 * Omega_j * nabla_j W(r_ij, sqrt{2}*h_j)
```

### Momentum Equation

```
dS_i/dt = -sum_j m_j * F_ij * [
    -(P_t^* - P_B_parallel) * n_ij
    + B_parallel^* * B_perp^*
]

P_B_parallel = (B_parallel_i^2 + B_parallel_j^2) / 4
B_parallel^* = (B_parallel_i + B_parallel_j) / 2
```

### Energy Equation

```
de_i/dt = -sum_j m_j * F_ij * [
    -(P_t^* - P_B_parallel) * v_parallel^*
    + B_parallel^* * dot(B_perp^*, v_perp^*)
]
```

### Induction Equation

```
d(B/N)_i/dt = sum_j m_j * F_ij * B_parallel^* * [
    v_parallel^* * n_ij + v_perp^* - v_i^*
]
```

where `v_i^* = v_i + a_i * dt/2` (time-centered).

## Powell Divergence Cleaning

Add to momentum:
```
dS_i/dt -= B_i * sum_j m_j * B_parallel^* * F_ij
```

Add to energy:
```
de_i/dt -= dot(B_i, v_i^*) * sum_j m_j * B_parallel^* * F_ij
```

## Relativistic Riemann Solver

### Shock Jump Conditions

Post-shock enthalpy (quadratic):
```
(1 + A) H_b^2 - A H_b + (B H_a - H_a^2) = 0

A = (Gamma - 1)(P_a - P_b) / (Gamma * P_b)
B = (P_a - P_b) / n_a
```

Mass flux:
```
j^2 = -(P_b - P_a) / (H_b/n_b - H_a/n_a)
```

Shock speed:
```
V_s = [N_a^2 v_a +/- j * sqrt(j^2 + N_a^2(1-v_a^2))] / [N_a^2 + j^2]
```

Post-shock velocity:
```
v_b = [H_a gamma_a v_a + gamma_s(P_b-P_a)/j] /
      [H_a gamma_a + (P_b-P_a)(gamma_s v_a/j + 1/N_a)]
```

### Rarefaction Wave

Isentropic relations:
```
n(P) = (P/K)^{1/Gamma}     where K = P_a / n_a^Gamma
H(P) = 1 + [Gamma/(Gamma-1)] * P/n(P)
c_s(P) = sqrt(Gamma * P / (n * H))
```

Velocity integration:
```
v_b = tanh(B)

B = arctanh(v_a) + sign * integral_{P_a}^{P_b} f(P) dP

f(P) = sqrt(H^2 + A^2(1-c_s^2)) / [(H^2 + A^2) * n * c_s]

A = H_a * gamma_a * v_t_a   (tangent momentum invariant)
```

### HLLC Approximation

Wave speeds:
```
S_L = min(lambda_L^-, lambda_HLL)
S_R = max(lambda_R^+, lambda_HLL)

lambda^+/- = (v +/- c_s) / (1 +/- v*c_s)  (relativistic addition)
```

Interface pressure (acoustic):
```
Z = rho * H * gamma^2 * c_s   (acoustic impedance)

P^* = (Z_L P_R + Z_R P_L + Z_L Z_R (v_L - v_R)) / (Z_L + Z_R)
v^* = (Z_L v_L + Z_R v_R + P_L - P_R) / (Z_L + Z_R)
```

## Primitive Recovery

### Quartic for Lorentz Factor

```
(gamma^2 - 1)(X e gamma - 1)^2 = S^2 (X gamma^2 - 1)^2

X = Gamma / (Gamma - 1)
```

### Newton-Raphson

```
f(gamma) = (gamma^2 - 1)(Xe*gamma - 1)^2 - S^2(X*gamma^2 - 1)^2

df/dgamma = 2*gamma*(Xe*gamma-1)^2 + (gamma^2-1)*2*(Xe*gamma-1)*Xe
          - S^2 * 2*(X*gamma^2-1) * 2*X*gamma

gamma_{n+1} = gamma_n - f(gamma_n) / df(gamma_n)
```

### Recover Primitives

```
H = (X e gamma - 1) / (X gamma^2 - 1)
v = S / (gamma * H)
n = N / gamma
P = n(H - 1)(Gamma - 1) / Gamma
c_s = sqrt[(Gamma - 1)(H - 1) / H]
```

## Time Integration

Euler:
```
x_i^{n+1} = x_i^n + v_i^n * dt
S_i^{n+1} = S_i^n + dS_i/dt * dt
e_i^{n+1} = e_i^n + de_i/dt * dt
(B/N)_i^{n+1} = (B/N)_i^n + d(B/N)_i/dt * dt
```

## CFL Condition

```
dt = C_CFL * min_i(h_i / c_fast_i)

c_fast = sqrt[(c_s^2 + c_A^2 + sqrt((c_s^2 + c_A^2)^2 - 4*c_s^2*c_A^2*cos^2(theta))) / 2]
```

Typical: C_CFL = 0.3

## Shock Detection

Disable 2nd order reconstruction when:
```
C_shock * dot(e_ij, v_i - v_j) > min(c_s_i, c_s_j)

OR

|log10(P_i / P_j)| > C_cd
```

Typical: C_shock = 3, C_cd = 1.0

## Physical Constants

| Constant | Description | Typical Value |
|----------|-------------|---------------|
| Gamma | Adiabatic index | 5/3 (non-rel) or 4/3 (rel) |
| c | Speed of light | 1 (natural units) |
| eta | Smoothing length ratio | 1.0 |
| C_smooth | Volume sampling factor | 2.0 |
| C_CFL | CFL number | 0.3 |
