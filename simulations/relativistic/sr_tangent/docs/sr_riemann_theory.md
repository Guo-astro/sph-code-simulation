# Special Relativistic Riemann Solvers: Theory and Implementation

## Overview

This document describes the mathematical theory behind the Exact Riemann solver and HLLC approximate solver for Special Relativistic Hydrodynamics (SRHD), as implemented in the SR-GSPH code. The derivations follow Pons et al. (2000), Rezzolla & Zanotti (2003), Mignone & Bodo (2005), and Kitajima et al. (2025).

---

## 1. Relativistic Hydrodynamics Fundamentals

### 1.1 Basic Variables and Notation

| Symbol | Description |
|--------|-------------|
| $n$ | Rest-frame baryon number density |
| $N = \gamma n$ | Lab-frame baryon number density |
| $P$ | Pressure |
| $\rho$ | Rest-mass density (often used interchangeably with $n$) |
| $h$ | Specific enthalpy per baryon |
| $u$ | Specific internal energy per baryon |
| $\gamma$ (or $\Gamma$) | Adiabatic index (ratio of specific heats) |
| $v^x$ | Normal velocity component |
| $v^t$ | Tangential velocity component |
| $W = \gamma_L$ | Lorentz factor |
| $c_s$ | Sound speed |
| $c$ | Speed of light (set to 1 in code) |

### 1.2 Equation of State (Ideal Gas)

For an ideal relativistic gas with adiabatic index $\Gamma$:

**Step 1: Internal energy**
$$u = \frac{P}{(\Gamma - 1) n}$$

**Step 2: Specific enthalpy**
$$h = 1 + u + \frac{P}{n} = 1 + \frac{\Gamma}{\Gamma - 1} \frac{P}{n}$$

**Step 3: Sound speed**
$$c_s = \sqrt{\frac{\Gamma P}{n h}}$$

### 1.3 Lorentz Factor

The Lorentz factor includes ALL velocity components:

$$W = \frac{1}{\sqrt{1 - v^2}} = \frac{1}{\sqrt{1 - (v^x)^2 - (v^t)^2}}$$

This is critical for tangent velocity problems.

---

## 2. The Relativistic Riemann Problem

### 2.1 Problem Setup

Given two initial states separated at $x = 0$:
- **Left state (L)**: $(P_L, n_L, v^x_L, v^t_L)$
- **Right state (R)**: $(P_R, n_R, v^x_R, v^t_R)$

Find the solution at time $t > 0$, which consists of:
1. Left wave (shock or rarefaction)
2. Contact discontinuity
3. Right wave (shock or rarefaction)

### 2.2 Wave Structure

The solution is self-similar: all quantities depend only on $\xi = x/t$.

```
        t
        ^
        |     /   Left       | Contact |   Right     \
        |    /    Wave       |   (CD)  |    Wave      \
        |   /                |         |               \
        |  /    State L*     |         |    State R*    \
        | /                  |         |                 \
        |/___________________v_________v__________________\___> x
       Left State (L)                               Right State (R)
```

**Star region properties:**
- $P^* = P_L^* = P_R^*$ (pressure continuous across CD)
- $v^{x*} = v^{x*}_L = v^{x*}_R$ (normal velocity continuous across CD)
- $v^{t*}_L \neq v^{t*}_R$ (tangent velocity discontinuous across CD)

---

## 3. Exact Riemann Solver

### 3.1 Root-Finding Formulation

The exact solution requires finding $P^*$ such that:

$$f(P^*) = v^x_L(P^*) - v^x_R(P^*) = 0$$

where $v^x_L(P^*)$ and $v^x_R(P^*)$ are the normal velocities obtained by tracing the left and right wave curves from the initial states to pressure $P^*$.

### 3.2 Shock Wave Relations (Taub Adiabat)

When $P^* > P_a$ (ahead state), a shock forms. The relations follow from the relativistic Rankine-Hugoniot conditions.

**Step 1: Enthalpy behind shock (Taub Adiabat)**

From Pons et al. (2000) Eq. 27, solve the quadratic:

$$a H_b^2 + b H_b + c = 0$$

where:
- $a = 1 + \frac{(\Gamma - 1)(P_a - P_b)}{\Gamma P_b}$
- $b = -\frac{(\Gamma - 1)(P_a - P_b)}{\Gamma P_b}$
- $c = \frac{H_a (P_a - P_b)}{n_a} - H_a^2$

Solution:
$$H_b = \frac{-b + \sqrt{b^2 - 4ac}}{2a}$$

**Step 2: Density behind shock**

From the equation of state:
$$n_b = \frac{\Gamma}{\Gamma - 1} \frac{P_b}{H_b - 1}$$

**Step 3: Mass flux**

$$j^2 = -\frac{P_b - P_a}{\frac{H_b}{n_b} - \frac{H_a}{n_a}}$$

**Step 4: Shock speed**

From Pons et al. (2000) Eq. 26:
$$V_s = \frac{N_a^2 v^x_a \pm j\sqrt{j^2 + N_a^2(1 - (v^x_a)^2)}}{N_a^2 + j^2}$$

where $N_a = n_a W_a$ is the lab-frame density, and $\pm$ depends on wave direction.

**Step 5: Normal velocity behind shock**

$$v^x_b = \frac{H_a W_a v^x_a + W_s \frac{P_b - P_a}{\pm j}}{H_a W_a + (P_b - P_a)\left(\frac{W_s v^x_a}{\pm j} + \frac{1}{N_a}\right)}$$

**Step 6: Tangent velocity behind shock (Pons et al. Eq. 25)**

The tangent velocity changes across a shock:
$$v^t_b = v^t_a \cdot H_a W_a \sqrt{\frac{1 - (v^x_b)^2}{H_b^2 + (H_a W_a v^t_a)^2}}$$

### 3.3 Rarefaction Wave Relations

When $P^* < P_a$, a rarefaction wave forms. The solution uses Riemann invariants.

#### 3.3.1 The K-Invariant (Pons et al. 2000)

**Critical concept:** The quantity
$$A = h W v^t$$
is conserved along particle paths through a rarefaction wave.

This is the "K-invariant" or conserved tangential momentum.

#### 3.3.2 Isentropic Relations

Entropy is constant across rarefaction:
$$K_s = \frac{P}{n^\Gamma} = \text{const}$$

Therefore:
$$n_b = \left(\frac{P_b}{K_s}\right)^{1/\Gamma}$$

#### 3.3.3 Normal Velocity via Arctanh Transformation

**This is where Gauss-Legendre quadrature is used.**

From Pons et al. (2000) and the SRRP library:

**Step 1: Define the integrand**

$$I(P) = \frac{\sqrt{h^2 + A^2(1 - c_s^2)}}{(h^2 + A^2) \cdot n \cdot c_s}$$

where all quantities ($h$, $n$, $c_s$) are functions of $P$ via the isentropic relations.

**Step 2: Compute the integral**

$$\int_{P_a}^{P_b} I(P) \, dP$$

**Step 3: Apply arctanh transformation**

$$B = \frac{1}{2} \ln\left(\frac{1 + v^x_a}{1 - v^x_a}\right) + \text{sign} \cdot \int_{P_a}^{P_b} I(P) \, dP$$

where $\text{sign} = -1$ for left wave, $+1$ for right wave.

**Step 4: Recover normal velocity**

$$v^x_b = \tanh(B)$$

#### 3.3.4 Gauss-Legendre Quadrature

The integral is computed using 16-point Gauss-Legendre quadrature:

$$\int_a^b f(x) dx \approx \frac{b-a}{2} \sum_{i=1}^{16} w_i \cdot f\left(\frac{b-a}{2} t_i + \frac{a+b}{2}\right)$$

where $t_i$ are the Gauss nodes on $[-1, 1]$ and $w_i$ are the weights.

**Gauss-Legendre nodes (16-point):**
```
t = [-0.9894, -0.9446, -0.8656, -0.7554, -0.6179, -0.4580, -0.2816, -0.0950,
      0.0950,  0.2816,  0.4580,  0.6179,  0.7554,  0.8656,  0.9446,  0.9894]
```

**Gauss-Legendre weights (16-point):**
```
w = [0.0272, 0.0623, 0.0952, 0.1246, 0.1496, 0.1692, 0.1826, 0.1895,
     0.1895, 0.1826, 0.1692, 0.1496, 0.1246, 0.0952, 0.0623, 0.0272]
```

#### 3.3.5 Tangent Velocity from K-Invariant

After finding $v^x_b$, recover $v^t_b$:

$$v^t_b = A \sqrt{\frac{1 - (v^x_b)^2}{h_b^2 + A^2}}$$

### 3.4 Root-Finding Algorithm

**Step 1: Initial bracket**

Set $P_{lo} = P_{min} \times 10^{-6}$ and $P_{hi} = P_{max} \times 10^6$.

**Step 2: Evaluate residual**

$$f(P) = v^x_L(P) - v^x_R(P)$$

**Step 3: Hybrid Newton-Bisection**

```python
P_star = sqrt(P_lo * P_hi)  # Geometric mean (log-scale)
for iter in range(max_iter):
    f_mid = f(P_star)
    if |f_mid| < tol:
        break

    # Update bracket
    if f_mid * f_lo < 0:
        P_hi = P_star
    else:
        P_lo = P_star

    # Newton step
    df = (f(P_star + dP) - f(P_star - dP)) / (2*dP)
    P_newton = P_star - f_mid / df

    # Accept Newton if in bracket, else bisect
    if P_lo < P_newton < P_hi:
        P_star = P_newton
    else:
        P_star = sqrt(P_lo * P_hi)
```

### 3.5 Contact Discontinuity

Once $P^*$ and $v^{x*}$ are found:

**Tangent velocity selection (upwinding):**
- If $v^{x*} > 0$: use $v^{t*} = v^{t*}_L$
- If $v^{x*} < 0$: use $v^{t*} = v^{t*}_R$

---

## 4. HLLC Approximate Riemann Solver

### 4.1 Overview

HLLC (Harten-Lax-van Leer-Contact) is a 3-wave approximate solver based on Mignone & Bodo (2005).

**Key advantage:** Explicit wave speed estimates, which may provide better shock capturing in SPH.

### 4.2 Wave Speed Estimates

**Step 1: Relativistic characteristic speeds**

$$\lambda_\pm = \frac{v \pm c_s}{1 \pm v \cdot c_s}$$

**Step 2: Davis-type estimates**

$$S_L = \min\left(\lambda^-_L, \frac{v_L + v_R}{2} - \max(c_{s,L}, c_{s,R})\right)$$

$$S_R = \max\left(\lambda^+_R, \frac{v_L + v_R}{2} + \max(c_{s,L}, c_{s,R})\right)$$

### 4.3 Acoustic Approximation for Star State

**Step 1: Acoustic impedances**

$$Z_L = n_L h_L W_L^2 c_{s,L}$$
$$Z_R = n_R h_R W_R^2 c_{s,R}$$

**Step 2: Contact speed (normal velocity)**

$$v^{x*} = \frac{Z_L v^x_L + Z_R v^x_R + P_L - P_R}{Z_L + Z_R}$$

**Step 3: Star pressure**

$$P^* = \frac{Z_L P_R + Z_R P_L + Z_L Z_R (v^x_L - v^x_R)}{Z_L + Z_R}$$

### 4.4 Tangent Velocity in HLLC

**Step 1: Compute K-invariants**

$$K_L = h_L W_L v^t_L$$
$$K_R = h_R W_R v^t_R$$

**Step 2: Upwind selection**

$$K^* = \begin{cases} K_L & \text{if } v^{x*} > 0 \\ K_R & \text{if } v^{x*} \leq 0 \end{cases}$$

**Step 3: Estimate star enthalpy**

$$n^* \approx n_{upwind} \left(\frac{P^*}{P_{upwind}}\right)^{1/\Gamma}$$

$$h^* = 1 + \frac{\Gamma}{\Gamma - 1} \frac{P^*}{n^*}$$

**Step 4: Solve for tangent velocity**

From $K^* = h^* W^* v^{t*}$ and $W^* = 1/\sqrt{1 - (v^{x*})^2 - (v^{t*})^2}$:

$$(v^{t*})^2 = \frac{(K^*)^2 (1 - (v^{x*})^2)}{(K^*)^2 + (h^*)^2}$$

$$v^{t*} = \text{sign}(K^*) \sqrt{(v^{t*})^2}$$

---

## 5. Comparison: Exact vs HLLC

| Aspect | Exact Solver | HLLC Solver |
|--------|--------------|-------------|
| Accuracy | Exact (to numerical precision) | Approximate |
| Cost | Iterative (many function evaluations) | Direct (few evaluations) |
| Robustness | May fail in extreme cases | Generally robust |
| Shock capturing | Depends on SPH discretization | Explicit wave speeds may help |
| Rarefaction | Uses Gauss quadrature | Acoustic approximation |

### 5.1 When HLLC May Outperform Exact

At **low resolution**, HLLC's explicit wave speed estimates can provide better shock capturing because:
1. The acoustic impedance-based estimate is tuned for numerical stability
2. SPH discretization errors may dominate Riemann solver accuracy anyway

This explains why HLLC showed +15.2% error vs Exact's +22.9% at 800+800 resolution in the Kitajima tangent velocity test.

---

## 6. Implementation Notes

### 6.1 C++ Implementation (`sr_exact_riemann.cpp`)

Key functions:
- `solve_shock()`: Shock wave relations
- `solve_rarefaction()`: Rarefaction with Gauss quadrature
- `wave_curve()`: Unified interface for either wave type
- `exact_riemann_solver()`: Main root-finding loop

### 6.2 Python Implementation (SRRP library)

Key classes:
- `Shock`: Shock wave computations
- `Rarefaction`: Rarefaction with scipy.integrate.quad
- `Solver`: Wave pattern determination and root finding

### 6.3 Critical Implementation Details

1. **Always include tangent velocity in Lorentz factor**
2. **Use arctanh transformation for numerical stability**
3. **16-point Gauss-Legendre for accurate rarefaction integration**
4. **Geometric mean for pressure bracket (log-scale)**
5. **Upwind tangent velocity based on contact wave direction**

---

## References

1. Pons, J. A., Martí, J. M., & Müller, E. (2000). "The exact solution of the Riemann problem with non-zero tangential velocities in relativistic hydrodynamics." *J. Fluid Mech.*, 422, 125-139.

2. Rezzolla, L., & Zanotti, O. (2003). "An improved exact Riemann solver for relativistic hydrodynamics." *J. Fluid Mech.*, 479, 199-219.

3. Mignone, A., & Bodo, G. (2005). "An HLLC Riemann solver for relativistic flows - I. Hydrodynamics." *MNRAS*, 364, 126-136.

4. Kitajima, Y., et al. (2025). "Special Relativistic Godunov SPH." *arXiv:2510.18251v1*.
