# SR-GSPH Fluid Force Formulas - Complete Reference

This document provides explicit formula references for every equation in `sr_fluid_force.cpp`.

---

## 1. SPH EVOLUTION EQUATIONS

### Momentum Evolution
**Kitajima Eq. 209**:
```
⟨dS_i/dt⟩ = -Σ_j ν P*_ij V²_ij,interp [∇_i W - ∇_j W]
```

**Plain text**:
```
dS_i/dt = -sum_j { nu_j * P_star * V_sq_interp * (grad_i_W - grad_j_W) }
```

**Code location**: sr_fluid_force.cpp:461
```cpp
dS_dt -= f;  // where f = dw_ij * (nu_j * pstar * N_prod_inv)
```

### Energy Evolution
**Kitajima Eq. 212**:
```
⟨de_i/dt⟩ = -Σ_j ν P*_ij v*_ij · V²_ij,interp [∇_i W - ∇_j W]
```

**Plain text**:
```
de_i/dt = -sum_j { nu_j * P_star * v_star · V_sq_interp * (grad_i_W - grad_j_W) }
```

**Code location**: sr_fluid_force.cpp:465
```cpp
de_dt -= inner_product(f, v_star);  // where v_star = e_ij * vstar
```

### Volume Interpolation (Numerical Stability Approximation)
**Kitajima Eq. 209 (strict GSPH)**:
```
V²_ij,interp = 0.5 * (1/N_i² + 1/N_j²)
```

**Code uses intermediate scaling** (for numerical stability):
```
V²_ij,interp ≈ 1/(N_i * N_j)
```

**Code location**: sr_fluid_force.cpp:426-428
```cpp
const real N_prod_inv = 1.0 / (N_i * N_j);  // Intermediate scaling
```

**Rationale**: The strict formula causes numerical instabilities. The intermediate scaling provides better stability while still maintaining correct physics scaling.

---

## 2. RIEMANN PROBLEM STRUCTURE

### Wave Speeds
**Pons et al. (1999) Section 3**:
- λ₁ = v - c_s/(1 + vc_s/c²)     (left-going rarefaction head)
- λ₂ = v + c_s/(1 - vc_s/c²)     (right-going rarefaction head)
- λ_shock                         (shock velocity)

### Star Region Variables
**Pons Eq. 3.1-3.2**:
- P* : pressure in star region
- v* : normal velocity in star region
- v^t : tangential velocity (changes across waves due to Lorentz coupling)

### Four-Wave Structure
```
Left State → Rarefaction/Shock → Left Star → Contact → Right Star → Rarefaction/Shock → Right State
    L               wave 1           L*        wave 2      R*           wave 3              R
```

---

## 3. RAREFACTION WAVE FORMULAS

### Isentropic Relations
**Pons Eq. 3.12-3.14**:

Pressure-density relation:
```
P/P_a = (n/n_a)^γ_c
```

Sound speed (relativistic):
```
c_s² = (γ_c - 1)(H - 1)/H    [Kitajima Eq. 121]
```

### Rarefaction Velocity Formula
**Pons Eq. 3.20** (adapted to SRGSPH notation):
```
v* = c · (A - 1)/(A + 1)
```

where:
```
A = [(1 + v_a/c)/(1 - v_a/c)] · [(√(γ-1) + c_s,a/c)/(√(γ-1) - c_s,a/c)]^(-sign·2/√(γ-1))
    · [(√(γ-1) - c_s/c)/(√(γ-1) + c_s/c)]^(sign·2/√(γ-1))
```

**Code location**: sr_fluid_force.cpp:981-999
```cpp
const real sqgl1 = std::sqrt(gamma_c - 1.0);
const real term1 = (1.0 + va / c) / (1.0 - va / c);
const real term2 = (sqgl1 + csa / c) / (sqgl1 - csa / c);
const real term3 = (sqgl1 - cs / c) / (sqgl1 + cs / c);
const real exponent = -sign * 2.0 / sqgl1;
const real A = term1 * std::pow(term2 * term3, exponent);
vel = c * (A - 1.0) / (A + 1.0);
```

---

## 4. SHOCK WAVE FORMULAS

### Taub Adiabat (Ideal Gas)
**Pons et al. (1999) riem.tex Eq. 4.16** (LaTeX source):
```
h_b² [1 + (γ-1)(p_a-p_b)/(γp_b)] - [(γ-1)(p_a-p_b)/(γp_b)]h_b + [h_a(p_a-p_b)/ρ_a - h_a²] = 0
```

**Converted to Kitajima notation** (H = enthalpy per baryon, n = rest density):
```
H² [1 + (γ_c-1)(P_a-P)/(γ_c·P)] - [(γ_c-1)(P_a-P)/(γ_c·P)]H + [H_a(P_a-P)/n_a - H_a²] = 0
```

**Quadratic form**: a·H² + b·H + c = 0

where:
```
a = 1 + (γ_c-1)(P_a - P)/(γ_c·P)
b = -(γ_c-1)(P_a - P)/(γ_c·P) = -(a-1)
c = H_a(P_a - P)/n_a - H_a²
```

**Solution**:
```
H = [-b + √(b² - 4ac)] / (2a)
```

**Code location**: sr_fluid_force.cpp:870-887
```cpp
const real a = 1.0 + (gamma_c - 1.0) * (Pa - P) / (gamma_c * P);
const real b_term = -(gamma_c - 1.0) * (Pa - P) / (gamma_c * P);  // = -(a-1)
const real c_term = Ha * (Pa - P) / na - Ha * Ha;
const real discriminant = b_term * b_term - 4.0 * a * c_term;
H = (-b_term + std::sqrt(discriminant)) / (2.0 * a);
```

### Mass Flux
**Pons Eq. 4.7**:
```
j² = -[P]/[H/n] = -(P - P_a)/(H/n - H_a/n_a)
```

**Code location**: sr_fluid_force.cpp:904-913
```cpp
const real j2 = -(P - Pa) / (H / n - Ha / na);
```

### Shock Velocity
**Pons Eq. 5.22**:
```
V_s = [N_a² v_a + sign·|j|·√(j² + N_a²(1 - v_a²/c²))] / (j² + N_a²)
```

where N_a = W_a·n_a (lab-frame density)

**Code location**: sr_fluid_force.cpp:924-930
```cpp
const real Na = wa * na;
const real Na2 = Na * Na;
const real sqrt_arg = j2 + Na2 * (1.0 - va * va / (c * c));
const real vshock = (Na2 * va + sign * std::abs(j) * std::sqrt(sqrt_arg)) / (j2 + Na2);
```

### Post-Shock Normal Velocity
**Pons Eq. 6.1** (adapted):
```
v* = [γ_s(P-P_a)/j + H_a·W_a·v_a] / [H_a·W_a + (P-P_a)(γ_s·v_a/j + 1/(n_a·W_a))]
```

**Code location**: sr_fluid_force.cpp:937-940
```cpp
const real a_vel = gamma_shock * (P - Pa) / j + Ha * wa * va;
const real b_vel = Ha * wa + (P - Pa) * (gamma_shock * va / j + 1.0 / (na * wa));
vel = a_vel / b_vel;
```

---

## 5. ITERATIVE RIEMANN SOLVER

### Secant Method Iteration
**Pons Section 4.2**:

Initial guess:
```
P_guess = max(P_L, P_R)
```

Iteration:
```
f(P) = v*_L(P) - v*_R(P)
P_new = P - f(P) · (P - P_old) / [f(P) - f(P_old)]
```

Convergence:
```
|v*_L - v*_R| < tol
```

**Code location**: sr_fluid_force.cpp:661-766

---

## 6. TANGENTIAL VELOCITY COUPLING

### Lorentz Factor with Tangential Velocity
**Pons Eq. 2.8**:
```
γ = 1/√(1 - v²/c²) = 1/√(1 - (v^n)²/c² - (v^t)²/c²)
```

where:
- v^n = normal velocity component
- v^t = tangential velocity magnitude

### Post-Wave Tangential Velocity
**Pons Eq. 2.15**:
```
v^t_* = v^t_a · W_a/W_*
```

where W = γ (Lorentz factor)

**Code location**: sr_fluid_force.cpp:959-960
```cpp
const real vt_star = vt_mag_a * wa / w;  // v^t_* = v^t_a · γ_a/γ_*
```

---

## 7. PRIMITIVE VARIABLES

### Conserved Variables (per baryon)
**Kitajima**:
```
S = γ H v          [Eq. 104]
e = γ H - P/(Nc²)  [Eq. 107]
N = γ n            [Eq. 91]
```

### Ideal Gas EOS
**Kitajima**:
```
P = (γ_c - 1) n u              [Eq. 119]
H = 1 + u/c² + P/(nc²)         [Eq. 112]
c_s² = (γ_c - 1)(H - 1)/H      [Eq. 121]
```

---

## 8. FORMULA VERIFICATION CHECKLIST

For each equation in sr_fluid_force.cpp, verify:

1. **Sign consistency**:
   - Force evolution has correct ± sign
   - Shock jump conditions use correct [P] = P - P_a direction

2. **Units**:
   - All velocities in units of c (set c=1 in code)
   - Pressure, density, enthalpy dimensionally consistent
   - Factor of c² appears in energy/enthalpy where needed

3. **Lorentz factor**:
   - γ = 1/√(1 - v²/c²) includes ALL velocity components
   - Tangential velocities correctly coupled through γ

4. **Riemann solver**:
   - Left/right state assignment correct
   - Wave speeds have correct signs
   - Secant method converges to v*_L = v*_R

5. **Taub adiabat**:
   - Quadratic coefficients a, b, c match Pons
   - Positive root selected (H > H_a for shock)
   - Discriminant check prevents unphysical states
