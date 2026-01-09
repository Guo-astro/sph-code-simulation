# Step-by-Step Derivations for SRMHD-GSPH

This document provides detailed step-by-step derivations for the Special Relativistic Magnetohydrodynamics (SRMHD) formulation in Godunov SPH.

## Table of Contents

1. [Relativistic MHD Equations](#1-relativistic-mhd-equations)
2. [SPH Discretization via Convolution](#2-sph-discretization-via-convolution)
3. [Riemann Solver Derivation](#3-riemann-solver-derivation)
4. [Method of Characteristics for Alfven Waves](#4-method-of-characteristics-for-alfven-waves)
5. [Primitive Variable Recovery](#5-primitive-variable-recovery)
6. [Induction Equation Discretization](#6-induction-equation-discretization)

---

## 1. Relativistic MHD Equations

### 1.1 Four-Current Conservation

**Starting Point:** The baryon number 4-current conservation:

```
nabla_mu (n u^mu) = 0
```

where:
- `n` = rest-frame baryon number density
- `u^mu` = 4-velocity with `u^0 = gamma`, `u^i = gamma v^i`

**Step 1:** Expand in 3+1 form:

```
partial_t(n u^0) + partial_i(n u^i) = 0
partial_t(gamma n) + nabla . (gamma n v) = 0
partial_t(N) + nabla . (N v) = 0
```

where we define `N = gamma n` (lab-frame density).

**Step 2:** Convert to Lagrangian form using:

```
d/dt = partial_t + v . nabla
```

Therefore:

```
dN/dt + N nabla . v = 0
dN/dt = -N nabla . v
```

### 1.2 Energy-Momentum Tensor for SRMHD

**Starting Point:** The total energy-momentum tensor:

```
T^{mu nu} = T^{mu nu}_fluid + T^{mu nu}_EM
```

**Step 1:** Fluid contribution:

```
T^{mu nu}_fluid = rho h u^mu u^nu + P g^{mu nu}
```

where:
- `rho` = rest-mass density = `m_B n` (m_B = baryon mass)
- `h` = specific enthalpy per unit mass
- `H = h/c^2` = enthalpy per baryon (dimensionless with c=1)

**Step 2:** Electromagnetic contribution in ideal MHD:

```
T^{mu nu}_EM = (1/4pi) [F^{mu alpha} F_alpha^nu + (1/4) g^{mu nu} F^{alpha beta} F_{alpha beta}]
```

In ideal MHD with frozen-in condition `E + v x B = 0`:

```
T^{mu nu}_EM = b^2 u^mu u^nu + (b^2/2) g^{mu nu} - b^mu b^nu
```

where `b^mu` is the magnetic field 4-vector satisfying `b^mu u_mu = 0`:

```
b^0 = gamma (v . B) / c
b^i = B^i/gamma + gamma (v . B) v^i / c^2
b^2 = B^2/gamma^2 + (v . B)^2/c^2
```

**Step 3:** Combined tensor in lab frame (3+1 split, c=1):

For spatial components (i, j = 1, 2, 3):

```
T^{ij} = (rho H gamma^2 + B^2) v^i v^j + (P + B^2/2) delta^{ij} - B^i B^j
```

For momentum flux:

```
T^{0i} = (rho H gamma^2 + B^2) v^i
```

### 1.3 Momentum Equation Derivation

**Starting Point:** Energy-momentum conservation:

```
nabla_mu T^{mu nu} = 0
```

**Step 1:** Spatial components give momentum equation:

```
partial_t T^{0i} + nabla_j T^{ji} = 0
```

**Step 2:** Define canonical momentum per baryon:

```
S^i = (rho H gamma^2 + B^2/(4 pi)) v^i / N
    = gamma H v^i + (B^2/(4 pi N gamma)) v^i
```

For hydrodynamics only (B=0):

```
S = gamma H v
```

**Step 3:** Lagrangian form:

```
dS/dt = -(1/N) nabla_j T^{ji}
```

### 1.4 Energy Equation Derivation

**Starting Point:** Time component of energy-momentum conservation:

```
partial_t T^{00} + nabla_i T^{0i} = 0
```

**Step 1:** Define canonical energy per baryon:

For hydrodynamics:
```
e = gamma H - P/(N c^2)
```

With c=1:
```
e = gamma H - P/N
```

**Step 2:** Lagrangian form:

```
de/dt = -(1/N) nabla_i (T^{ij} v_j)
```

---

## 2. SPH Discretization via Convolution

### 2.1 Convolution Approach

**Key Insight (Kitajima et al.):** Instead of discretizing fields directly, we discretize the convolution of fields with the kernel.

**Definition:** For a field `f(x)`:

```
<f_i> = integral f(x) W(x - x_i, h(x)) d^3x
```

### 2.2 Volume-Based Density

**Step 1:** Define particle volume:

```
V_p(x) = [sum_j W(x - x_j, h(x))]^{-1}
```

**Step 2:** Number density from baryon number `nu`:

```
N(x) = nu(x) / V_p(x)
```

This is exact when each particle carries baryon number `nu_i`.

### 2.3 Variable Smoothing Length

**Step 1:** Define smoothing length relation:

```
h(x) = eta * [V_p^*(x)]^{1/d}
```

where:
```
V_p^*(x) = [sum_j W(x - x_j, C_smooth * h(x))]^{-1}
```

**Step 2:** Solve iteratively:

```python
for iteration in range(max_iter):
    sum_W_star = sum_j W(r_ij, C_smooth * h)
    V_p_star = 1 / sum_W_star
    h_new = eta * V_p_star^(1/d)
    if |h_new - h| / h < tolerance:
        break
    h = h_new
```

### 2.4 Momentum Equation Discretization

**Starting Point:** Lagrangian momentum equation:

```
dS/dt = -(1/N) nabla . T
```

**Step 1:** Take convolution:

```
<nu_i dS_i/dt> = integral nu(x) (dS/dt) W_i(x) d^3x
               = -integral (nu(x)/N(x)) nabla_j T^{muj} W_i(x) d^3x
```

**Step 2:** Use `nu/N = V_p` and integrate by parts:

```
<nu_i dS_i/dt> = integral V_p T^{muj} nabla_j W_i d^3x
```

**Step 3:** Insert resolution of identity:

```
1 = sum_j integral W(x - x_j, h) d^3x
```

Leads to:

```
<nu_i dS_i/dt> = sum_j integral T^{muj}(x) / N^2(x)
                 * [(nabla_i - nabla_j) W_i W_j] d^3x
```

**Step 4:** Saddle-point approximation at midpoint x_M:

```
<nu_i dS_i/dt> = -sum_j (T^{muj})_M (V_{ij})^2
                 * [nabla_i W(x_M - x_i, sqrt{2} h_i)
                    - nabla_j W(x_M - x_j, sqrt{2} h_j)]
```

The sqrt{2} factor comes from the convolution of two Gaussians.

### 2.5 Grad-h Correction

**Step 1:** Account for variable h:

```
Omega_i = 1 / (1 + (h_i / (d * N_i)) * dN/dh)
```

**Step 2:** Simplified form:

```
Omega_i = 1 / (1 + h_i * sum_j dW/dh / (d * sum_j W))
```

**Step 3:** Include in force calculation:

```
F_ij = V_i^2 * Omega_i * nabla_i W_i + V_j^2 * Omega_j * nabla_j W_j
```

---

## 3. Riemann Solver Derivation

### 3.1 Problem Setup

**Configuration:** Two uniform states separated at x=0:

```
        x < 0          |          x > 0
    --------------------|--------------------
    (rho_L, P_L, v_L)   |   (rho_R, P_R, v_R)
    (B_L)               |   (B_R)
```

### 3.2 Decomposition into Compressive and Alfvenic Parts

**Key Insight (Iwasaki & Inutsuka):** For B_parallel = 0, the Riemann problem simplifies.

**Step 1:** Project to interface normal `n_ij`:

```
v_parallel = v . n_ij
v_perp = v - v_parallel * n_ij
B_parallel = B . n_ij
B_perp = B - B_parallel * n_ij
```

**Step 2:** The Riemann problem depends only on:
- Compressive: `(rho, P_t, v_parallel)` where `P_t = P + B_perp^2/2`
- Alfvenic: `(v_perp, B_perp)` solved by MOC

### 3.3 Relativistic Shock Jump Conditions

**Starting Point:** Rankine-Hugoniot conditions in relativistic form.

**Step 1:** Mass flux conservation across shock:

```
[rho gamma (v - V_s)] = 0
```

where V_s = shock speed, brackets denote jump.

Define mass flux: `j = rho_a gamma_a (v_a - V_s) = rho_b gamma_b (v_b - V_s)`

**Step 2:** Momentum flux conservation:

```
[rho H gamma^2 v (v - V_s) + P] = 0
```

**Step 3:** Energy flux conservation:

```
[rho H gamma^2 (v - V_s) - P V_s] = 0
```

**Step 4:** Solve for post-shock enthalpy (Pons et al. 2000):

Quadratic for H_b:
```
(1 + A) H_b^2 - A H_b - H_a^2 + B H_a = 0
```

where:
```
A = (gamma_c - 1)(P_a - P_b) / (gamma_c * P_b)
B = (P_a - P_b) / n_a
```

Solution:
```
H_b = [-(-A) + sqrt(A^2 + 4(1+A)(H_a^2 - B*H_a))] / [2(1+A)]
```

**Step 5:** Mass flux from jump conditions:

```
j^2 = -(P_b - P_a) / (H_b/n_b - H_a/n_a)
```

**Step 6:** Shock speed:

```
V_s = [N_a^2 v_a +/- j * sqrt(j^2 + N_a^2(1-v_a^2))] / [N_a^2 + j^2]
```

Sign: `-` for left-going, `+` for right-going shock.

**Step 7:** Post-shock velocity:

```
v_b = [H_a gamma_a v_a + gamma_s (P_b - P_a) / j_signed] /
      [H_a gamma_a + (P_b - P_a)(gamma_s v_a / j_signed + 1/N_a)]
```

### 3.4 Relativistic Rarefaction Wave

**Step 1:** Riemann invariant:

```
J_+/- = (1/2) ln[(1+v)/(1-v)] +/- integral dP / (rho h c_s)
```

is constant along characteristics.

**Step 2:** Along isentrope (S = const):

```
P / n^gamma_c = const = K
```

Therefore:
```
n(P) = (P/K)^{1/gamma_c}
H(P) = 1 + (gamma_c/(gamma_c-1)) * P/n(P)
c_s(P) = sqrt(gamma_c * P / (n * H))
```

**Step 3:** Solve for v_b:

For zero tangent velocity (v_t = 0):
```
v_b = tanh[arctanh(v_a) + sign * integral_{P_a}^{P_b} dP / (rho h c_s)]
```

For nonzero tangent velocity:
```
B = arctanh(v_a) + sign * integral_{P_a}^{P_b} f(P) dP
v_b = tanh(B)
```

where:
```
f(P) = sqrt(H^2 + A^2(1 - c_s^2)) / [(H^2 + A^2) * rho * c_s]
A = H_a * gamma_a * v_t_a  (conserved tangent momentum invariant)
```

**Step 4:** Numerical integration using Gauss-Legendre quadrature:

```python
def rarefaction_integral(P_a, P_b, state_a, sign, A_squared):
    # Transform [P_a, P_b] to [-1, 1]
    P_mid = 0.5 * (P_a + P_b)
    half_dP = 0.5 * (P_b - P_a)

    integral = 0.0
    for k in range(N_gauss):
        t = gauss_nodes[k]
        P = P_mid + half_dP * t

        rho = isentrope_density(P)
        H = isentrope_enthalpy(P, rho)
        c_s = isentrope_sound_speed(P, rho, H)

        f = sqrt(H^2 + A_squared*(1 - c_s^2)) / ((H^2 + A_squared) * rho * c_s)
        integral += gauss_weights[k] * f

    return half_dP * integral
```

### 3.5 Complete Riemann Solver

**Step 1:** Define wave curves:

```
v_L(P*) = velocity behind left wave at pressure P*
v_R(P*) = velocity behind right wave at pressure P*
```

**Step 2:** Find P* such that:

```
f(P*) = v_L(P*) - v_R(P*) = 0
```

**Step 3:** Newton-Raphson with bisection fallback:

```python
def solve_riemann(state_L, state_R):
    # Initial bracket
    P_lo = min(P_L, P_R) * 1e-6
    P_hi = max(P_L, P_R) * 1e6

    # Geometric mean initial guess
    P_star = sqrt(P_lo * P_hi)

    for iter in range(max_iter):
        f = wave_curve_L(P_star) - wave_curve_R(P_star)

        if |f| < tolerance:
            break

        # Update bracket
        if f * f_lo < 0:
            P_hi = P_star
        else:
            P_lo = P_star

        # Try Newton step
        df = derivative(f, P_star)
        P_newton = P_star - f / df

        if P_lo < P_newton < P_hi:
            P_star = P_newton
        else:
            P_star = sqrt(P_lo * P_hi)  # Bisection fallback

    # Get interface velocity
    v_star = 0.5 * (wave_curve_L(P_star) + wave_curve_R(P_star))

    return P_star, v_star
```

---

## 4. Method of Characteristics for Alfven Waves

### 4.1 Alfven Wave Equations

**Starting Point:** Linearized MHD equations for transverse perturbations:

```
partial_t v_perp = (B_parallel / (4 pi rho)) nabla_parallel B_perp
partial_t B_perp = B_parallel nabla_parallel v_perp
```

### 4.2 Characteristic Form

**Step 1:** Combine into wave equation:

```
partial_t^2 v_perp = c_A^2 nabla_parallel^2 v_perp
```

where `c_A = |B_parallel| / sqrt(4 pi rho)` is Alfven speed.

**Step 2:** Characteristics:

```
dx/dt = +/- c_A
```

**Step 3:** Riemann invariants along characteristics:

```
R_+ = v_perp + B_perp / sqrt(4 pi rho)  (along dx/dt = +c_A)
R_- = v_perp - B_perp / sqrt(4 pi rho)  (along dx/dt = -c_A)
```

### 4.3 Solution at Interface

**Step 1:** Trace characteristics to interface:

- R_+ propagates from left state (if B_parallel > 0)
- R_- propagates from right state (if B_parallel > 0)

**Step 2:** At interface, both invariants must match:

```
R_+^* = v_perp^* + B_perp^* / sqrt(4 pi rho_L) = R_+^L
R_-^* = v_perp^* - B_perp^* / sqrt(4 pi rho_R) = R_-^R
```

**Step 3:** Solve for interface values:

Adding equations:
```
2 v_perp^* = R_+^L + R_-^R
v_perp^* = (v_perp_L sqrt(rho_L) + v_perp_R sqrt(rho_R)
            + sign(B_parallel)(B_perp_R - B_perp_L)) / (sqrt(rho_L) + sqrt(rho_R))
```

Subtracting equations:
```
B_perp^* (1/sqrt(rho_L) + 1/sqrt(rho_R)) = R_+^L - R_-^R
B_perp^* = (B_perp_L/sqrt(rho_L) + B_perp_R/sqrt(rho_R)
            + sign(B_parallel)(v_perp_R - v_perp_L)) / (1/sqrt(rho_L) + 1/sqrt(rho_R))
```

### 4.4 Relativistic Extension

**Step 1:** Replace `rho` with relativistic wave impedance:

```
Z = rho H gamma^2
```

**Step 2:** Modified formulas:

```
v_perp^* = (v_perp_L sqrt(Z_L) + v_perp_R sqrt(Z_R)
            + sign(B_parallel)(B_perp_R - B_perp_L)) / (sqrt(Z_L) + sqrt(Z_R))

B_perp^* = (B_perp_L/sqrt(Z_L) + B_perp_R/sqrt(Z_R)
            + sign(B_parallel)(v_perp_R - v_perp_L)) / (1/sqrt(Z_L) + 1/sqrt(Z_R))
```

---

## 5. Primitive Variable Recovery

### 5.1 Problem Statement

**Given:** Conserved variables (S, e, N, B)
**Find:** Primitive variables (v, n, P, gamma)

### 5.2 Derivation of Quartic for gamma

**Step 1:** Start with definitions:

```
S = gamma H v
e = gamma H - P/(N c^2)
```

**Step 2:** From EOS:

```
H = 1 + gamma_c/(gamma_c - 1) * P/(n c^2)
  = 1 + X * P/(n c^2)
```

where `X = gamma_c/(gamma_c - 1)`.

**Step 3:** Express P in terms of H and n:

```
P = n c^2 (H - 1) / X
```

**Step 4:** Use `n = N/gamma`:

```
P = N c^2 (H - 1) / (X gamma)
```

**Step 5:** Substitute into energy equation:

```
e = gamma H - (N c^2 (H-1)) / (X gamma N c^2)
  = gamma H - (H-1)/(X gamma)
  = [X gamma^2 H - (H-1)] / (X gamma)
```

Solve for H:
```
X e gamma = X gamma^2 H - (H-1)
X e gamma = H (X gamma^2 - 1) + 1
H = (X e gamma - 1) / (X gamma^2 - 1)
```

**Step 6:** From momentum equation:

```
|v| = |S| / (gamma H)
```

**Step 7:** Use Lorentz factor constraint:

```
gamma^2 = 1 / (1 - v^2)
v^2 = 1 - 1/gamma^2 = (gamma^2 - 1)/gamma^2
```

Therefore:
```
|S|^2 = gamma^2 H^2 v^2 = H^2 (gamma^2 - 1)
```

**Step 8:** Substitute H:

```
|S|^2 = [(X e gamma - 1)^2 / (X gamma^2 - 1)^2] * (gamma^2 - 1)
```

Rearrange:
```
|S|^2 (X gamma^2 - 1)^2 = (X e gamma - 1)^2 (gamma^2 - 1)
```

This is the quartic equation in gamma.

### 5.3 Newton-Raphson Solution

```python
def solve_gamma(S_mag, e, N, gamma_c):
    X = gamma_c / (gamma_c - 1)

    def f(gamma):
        term1 = (gamma^2 - 1) * (X * e * gamma - 1)^2
        term2 = S_mag^2 * (X * gamma^2 - 1)^2
        return term1 - term2

    def df(gamma):
        A = gamma^2 - 1
        B = X * e * gamma - 1
        dA = 2 * gamma
        dB = X * e

        C = X * gamma^2 - 1
        dC = 2 * X * gamma

        d_term1 = dA * B^2 + A * 2 * B * dB
        d_term2 = S_mag^2 * 2 * C * dC
        return d_term1 - d_term2

    # Initial guess
    v_guess = S_mag / e if e > 0 else 0.5
    v_guess = min(v_guess, 0.99)
    gamma = 1 / sqrt(1 - v_guess^2)

    for iter in range(max_iter):
        delta = f(gamma) / df(gamma)
        gamma -= delta
        gamma = max(gamma, 1.0)  # Physical constraint

        if |delta| < tolerance:
            break

    return gamma
```

### 5.4 With Tangent Velocity

**Modified constraint:**

```
gamma = 1 / sqrt(1 - v_x^2 - v_t^2)
```

**Coupled system:**

```
v_x = S_x / (gamma H)
v_t = S_t / (gamma H)
gamma = 1 / sqrt(1 - (S_x^2 + S_t^2) / (gamma^2 H^2))
```

**Iterative solution:**

```python
def recover_with_tangent(S_x, S_t, e, N, gamma_c):
    S_total = sqrt(S_x^2 + S_t^2)

    # Initial estimate
    gamma = solve_gamma(S_total, e, N, gamma_c)

    for iter in range(max_iter):
        H = (X * e * gamma - 1) / (X * gamma^2 - 1)
        H = max(H, 1.0 + epsilon)

        gamma_H = gamma * H
        v_x = S_x / gamma_H
        v_t = S_t / gamma_H

        v2 = v_x^2 + v_t^2
        v2 = min(v2, 0.9999)  # Subluminal

        gamma_new = 1 / sqrt(1 - v2)

        if |gamma_new - gamma| / gamma < tolerance:
            break

        gamma = gamma_new

    return v_x, v_t, gamma, H
```

---

## 6. Induction Equation Discretization

### 6.1 Continuum Form

**Ideal MHD:** Magnetic field frozen into fluid.

```
partial_t B + nabla x E = 0
E = -v x B  (ideal MHD)
```

Therefore:
```
partial_t B = nabla x (v x B)
            = (B . nabla) v - (v . nabla) B - B (nabla . v) + v (nabla . B)
```

For `nabla . B = 0`:
```
partial_t B = (B . nabla) v - (v . nabla) B - B (nabla . v)
```

### 6.2 Lagrangian Form

**Step 1:** Using `d/dt = partial_t + v . nabla`:

```
dB/dt = (B . nabla) v - B (nabla . v)
```

**Step 2:** For `B/rho` (or `B/N`):

```
d(B/N)/dt = (1/N) dB/dt - B/N^2 dN/dt
          = (1/N)[(B . nabla) v - B (nabla . v)] + (B/N)(nabla . v)
          = (1/N)(B . nabla) v
```

### 6.3 SPH Discretization

**Step 1:** Convolution form:

```
<nu_i d(B/N)_i/dt> = integral nu(x)/N(x) (B . nabla) v W_i(x) d^3x
```

**Step 2:** Discretize:

```
d(B/N)_i/dt = sum_j m_j (B . n_ij)^* (v^* - v_i^*) F_ij
```

where:
- `F_ij` = kernel gradient factor
- `(B . n_ij)^* = (B_parallel,i + B_parallel,j)/2`
- `v^*` = interface velocity from Riemann/MOC
- `v_i^*` = time-centered velocity for particle i

### 6.4 Powell Correction

**Motivation:** Numerical errors can create `nabla . B != 0`, leading to:
1. Tensile instability in low-beta plasma
2. Incorrect wave propagation
3. Artificial forces along field lines

**Powell source terms:** Add to evolution equations:

```
dS_i/dt -= B_i * (nabla . B)_i
de_i/dt -= (B_i . v_i) * (nabla . B)_i
```

**SPH form of div(B):**

```
(nabla . B)_i = sum_j m_j B_parallel^* F_ij
```

This breaks exact conservation but maintains stability.

---

## Summary

This document derived:

1. **Conservation laws** from relativistic MHD in Lagrangian form
2. **SPH discretization** using convolution integrals and saddle-point approximation
3. **Riemann solver** for relativistic shocks and rarefactions with tangent velocity
4. **Method of Characteristics** for Alfven waves with relativistic extension
5. **Primitive recovery** algorithm including quartic solver
6. **Induction equation** discretization with Powell corrections

The complete algorithm enables simulation of relativistic MHD flows with:
- Shock capturing without artificial viscosity
- Proper Alfven wave dynamics
- Conservation of magnetic flux
- Numerical stability in low-beta plasmas
