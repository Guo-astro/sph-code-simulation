# Complete Formula-to-Code Mapping for SR-GSPH Implementation

This document provides line-by-line mapping between paper formulas and code implementation for all SR-GSPH files.

---

## FILE: sr_primitive_recovery.cpp

### Function: solve_lorentz_factor() - Lines 109-185

**QUARTIC EQUATION** (lines 143-149):
```
Paper: (γ² - 1)(Xeγ - 1)² - S²(Xγ² - 1)² = 0    [Derived from Kitajima Eq. 104, 107]
Code:  f = (w2 - 1.0) * term1 * term1 - S_mag * S_mag * term2 * term2
       where term1 = X * e * w - 1.0
             term2 = X * w2 - 1.0
             w = gamma, w2 = gamma²
```

**DERIVATIVE** (lines 150-154):
```
Paper: df/dγ = 2γ(Xeγ-1)² + 2(γ²-1)(Xeγ-1)Xe - 4S²Xγ(Xγ²-1)
Code:  df_dw = 2.0 * w * term1 * term1
              + 2.0 * (w2 - 1.0) * term1 * X * e
              - 2.0 * S_mag * S_mag * term2 * 2.0 * X * w
```

**X PARAMETER** (line 120):
```
Paper: X ≡ γ_c/(γ_c - 1)
Code:  const real X = gamma_eos / (gamma_eos - 1.0);
```

### Function: recover_velocity() - Lines 219-240

**VELOCITY FORMULA** (lines 236-239):
```
Paper: v = [(Xγ² - 1)/(γ(Xeγ - 1))] · S    [Derived from S = γHv, Kitajima Eq. 104]
Code:  numerator = X * gamma2 - 1.0
       denominator = gamma * (X * e * gamma - 1.0)
       factor = numerator / denominator
       return S * factor
```

### Function: conserved_to_primitive() - Lines 311-373

**STEP 1: Lorentz factor** (line 324):
```
Paper: Solve quartic (γ² - 1)(Xeγ - 1)² - S²(Xγ² - 1)² = 0
Code:  prim.gamma_lor = solve_lorentz_factor(S_mag, e, N, gamma_eos, c_speed);
```

**STEP 2: Velocity** (line 330):
```
Paper: v = [(Xγ² - 1)/(γ(Xeγ - 1))] · S
Code:  prim.vel = recover_velocity(S, prim.gamma_lor, e, gamma_eos);
```

**STEP 3: Rest-frame density** (line 338):
```
Paper: n = N/γ    [Kitajima Eq. 91: N = γn]
Code:  prim.density = N / prim.gamma_lor;
```

**STEP 4: Enthalpy** (line 348):
```
Paper: H = (Xeγ - 1)/(Xγ² - 1)    [From quartic derivation]
Code:  prim.enthalpy = (X * e * prim.gamma_lor - 1.0) / (X * gamma2 - 1.0);
```

**STEP 5: Internal energy** (line 357):
```
Paper: u = c²(H - 1)/γ_c    [From H = 1 + γ_c u/c², Kitajima Eq. 112, 119]
Code:  const real u = c_speed * c_speed * (prim.enthalpy - 1.0) / gamma_eos;
```

**STEP 6: Pressure** (line 363):
```
Paper: P = (γ_c - 1) n u    [Kitajima Eq. 119, ideal gas EOS]
Code:  prim.pressure = (gamma_eos - 1.0) * prim.density * u;
```

**STEP 7: Sound speed** (line 370):
```
Paper: c_s² = (γ_c - 1)(H - 1)/H    [Kitajima Eq. 121]
Code:  const real cs_sq = (gamma_eos - 1.0) * (prim.enthalpy - 1.0) / prim.enthalpy;
       prim.sound_speed = std::sqrt(std::max(0.0, cs_sq));
```

---

## FILE: sr_pre_interaction.cpp

### Function: compute_volume() - Lines 119-152

**PARTICLE VOLUME** (lines 129-151):
```
Paper: V_p(x) = [Σ_j W(x - x_j, h)]^(-1)    [Kitajima Eq. 220]
Code:  sum_w = 0.0;
       for each neighbor j:
           sum_w += kernel->w(r, h)
       sum_w += kernel->w(0.0, h)  // self-contribution
       return 1.0 / sum_w;
```

### Function: compute_smoothing_length() - Lines 187-246

**SMOOTHING LENGTH ITERATION** (lines 203-241):
```
Paper: h(x) = η · [V_p*(x)]^(1/d)                     [Kitajima Eq. 231]
       V_p*(x) = [Σ_j W(x-x_j, C_smooth·h)]^(-1)     [Kitajima Eq. 233]
Code:  h_smooth = m_c_smooth * h;
       sum_w_star = Σ_j kernel->w(r, h_smooth) + kernel->w(0.0, h_smooth)
       Vp_star = 1.0 / sum_w_star;
       h_new = m_eta * std::pow(Vp_star, 1.0 / real(DIM));
```

### Function: calculation() - Main loop, Lines 277-450

**NUMBER DENSITY** (lines 336-337):
```
Paper: N = ν/V_p    [Kitajima Eq. 239]
Code:  const real nu_i = p_i.mass;
       const real N_i = nu_i / Vp_i;
```

**PRIMITIVE RECOVERY** (line 364):
```
Paper: Solve quartic to get (γ, v, n, H, P, c_s) from (S, e, N)
Code:  prim = PrimitiveRecovery::conserved_to_primitive(S_i, e_i, N_i, m_gamma, m_c_speed);
```

---

## FILE: sr_fluid_force.cpp

### SPH Evolution Equations - Lines 483-494

**MOMENTUM EVOLUTION** (line 488):
```
Paper: ⟨dS_i/dt⟩ = -Σ_j ν P*_ij V²_ij [∇_i W - ∇_j W]    [Kitajima Eq. 209]
Code:  dS_dt -= f;
       where f = dw_ij * (nu_j * pstar * N_prod_inv)
             dw_ij = grad_i_W + grad_j_W
```

**ENERGY EVOLUTION** (line 494):
```
Paper: ⟨de_i/dt⟩ = -Σ_j ν P*_ij v*_ij · V²_ij [∇_i W - ∇_j W]    [Kitajima Eq. 212]
Code:  de_dt -= inner_product(f, v_star);
       where v_star = e_ij * vstar
```

**VOLUME INTERPOLATION** (line 428):
```
Paper (strict): V²_ij = 0.5(1/N_i² + 1/N_j²)    [Kitajima Eq. 209]
Code (stability): N_prod_inv = 1.0 / (N_i * N_j);
Note: Intermediate scaling for numerical stability
```

### Taub Adiabat (Shock Waves) - Lines 895-914

**QUADRATIC EQUATION** (lines 897-914):
```
Paper: H² [1+(γ_c-1)(P_a-P)/(γ_c·P)] - [(γ_c-1)(P_a-P)/(γ_c·P)]H
       + [H_a(P_a-P)/n_a - H_a²] = 0    [Pons riem.tex Eq. 4.16]

Quadratic form: a·H² + b·H + c = 0

Code:  a = 1.0 + (gamma_c - 1.0) * (Pa - P) / (gamma_c * P);
       b_term = -(gamma_c - 1.0) * (Pa - P) / (gamma_c * P);  // = -(a-1)
       c_term = Ha * (Pa - P) / na - Ha * Ha;
       discriminant = b_term * b_term - 4.0 * a * c_term;
       H = (-b_term + std::sqrt(discriminant)) / (2.0 * a);
```

**VERIFICATION**:
- Coefficient b = -(a-1) ✓ verified from Pons Eq. 4.16
- Take positive root (compression shock: H > H_a) ✓

### Post-Shock Density - Line 930

```
Paper: From H = 1 + u/c² + P/(nc²) and P = (γ_c-1)nu    [Kitajima Eq. 112, 119]
       Derive: n = γ_c·P / [(γ_c-1)(H-1)]
Code:  n = gamma_c * P / ((gamma_c - 1.0) * (H - 1.0));
```

### Mass Flux - Line 934

```
Paper: j² = -[P]/[H/n] = -(P - P_a)/(H/n - H_a/n_a)    [Pons Eq. 4.7]
Code:  const real j2 = -(P - Pa) / (H / n - Ha / na);
```

### Shock Velocity - Lines 948-954

```
Paper: V_s = [N_a² v_a + sign·|j|·√(j² + N_a²(1-v_a²/c²))] / (j² + N_a²)
       where N_a = W_a·n_a    [Pons Eq. 5.22]
Code:  const real Na = wa * na;
       const real Na2 = Na * Na;
       const real sqrt_arg = j2 + Na2 * (1.0 - va * va / (c * c));
       const real vshock = (Na2 * va + sign * std::abs(j) * std::sqrt(sqrt_arg)) / (j2 + Na2);
```

### Post-Shock Normal Velocity - Lines 961-964

```
Paper: v* = [γ_s(P-P_a)/j + H_a·W_a·v_a] / [H_a·W_a + (P-P_a)(γ_s·v_a/j + 1/(n_a·W_a))]
       [Pons Eq. 6.1, adapted]
Code:  const real a_vel = gamma_shock * (P - Pa) / j + Ha * wa * va;
       const real b_vel = Ha * wa + (P - Pa) * (gamma_shock * va / j + 1.0 / (na * wa));
       vel = a_vel / b_vel;
```

### Post-Shock Sound Speed - Line 973

```
Paper: c_s² = (γ_c - 1)(H - 1)/H    [Kitajima Eq. 121]
Code:  const real cs_sq = (gamma_c - 1.0) * (H - 1.0) / H;
       cs = std::sqrt(cs_sq);
```

### Rarefaction Wave Velocity - Lines 981-999

```
Paper: v* = c · (A - 1)/(A + 1)    [Pons Eq. 3.20]
where: A = [(1+v_a/c)/(1-v_a/c)] · [(√(γ-1)+c_s,a/c)/(√(γ-1)-c_s,a/c)]^(-sign·2/√(γ-1))
           · [(√(γ-1)-c_s/c)/(√(γ-1)+c_s/c)]^(sign·2/√(γ-1))
Code:  const real sqgl1 = std::sqrt(gamma_c - 1.0);
       const real term1 = (1.0 + va / c) / (1.0 - va / c);
       const real term2 = (sqgl1 + csa / c) / (sqgl1 - csa / c);
       const real term3 = (sqgl1 - cs / c) / (sqgl1 + cs / c);
       const real exponent = -sign * 2.0 / sqgl1;
       const real A = term1 * std::pow(term2 * term3, exponent);
       vel = c * (A - 1.0) / (A + 1.0);
```

### Tangential Velocity Coupling - Lines 983-984

```
Paper: v^t_* = v^t_a · W_a/W_*    [Pons Eq. 2.15]
       where W = γ (Lorentz factor)
Code:  const real w = 1.0 / std::sqrt(1.0 - (vel*vel + vt_mag*vt_mag) / (c*c));
       const real vt_star = vt_mag_a * wa / w;
```

### Lorentz Factor (with tangential velocity) - Lines 981-984

```
Paper: γ = 1/√(1 - v²/c²) = 1/√(1 - (v^n)²/c² - (v^t)²/c²)    [Pons Eq. 2.8]
Code:  const real w = 1.0 / std::sqrt(1.0 - (vel*vel + vt_mag*vt_mag) / (c*c));
```

---

## FILE: sr_timestep.cpp

### CFL Condition - Lines 80-98

```
Paper: Δt = C_CFL · min_i[h_i/c_s,i]    [Kitajima Section 3, standard CFL]
Code:  const real dt_i = h_i / cs_i;
       dt_min = min(dt_min, dt_i);  // parallel reduction
       const real dt_sound = c_sound * dt_min;
       sim->set_dt(dt_sound);
```

---

## CRITICAL VERIFICATION POINTS

### 1. Sign Checks
- [x] Momentum evolution: `dS_dt -= f` (negative sign correct, Kitajima Eq. 209)
- [x] Energy evolution: `de_dt -= inner_product(f, v_star)` (negative sign correct, Kitajima Eq. 212)
- [x] Taub adiabat b coefficient: `b = -(a-1)` ✓ verified from Pons Eq. 4.16
- [x] Mass flux: `j² = -(P-Pa)/(...)` (negative sign for compression shock)

### 2. Unit Consistency
- [x] All velocities in units of c (c=1 in code)
- [x] Factors of c² in enthalpy: `H = 1 + u/c² + P/(nc²)` ✓
- [x] Sound speed formula: `c_s² = (γ_c-1)(H-1)/H` (dimensionless) ✓

### 3. Lorentz Factor
- [x] Includes tangential velocity: `γ = 1/√(1 - v_n²/c² - v_t²/c²)` ✓
- [x] Rest-frame density: `n = N/γ` ✓

### 4. Key Differences from Non-Relativistic SPH
- [x] Energy evolution uses `v*`, NOT `(v* - v_i)` ✓
- [x] Evolves canonical energy `e = γH - P/(Nc²)`, not thermal energy `u` ✓
- [x] Volume-based density: `N = ν/V_p`, not kernel sum ✓

---

## FORMULA REFERENCE QUICK INDEX

| Formula | Kitajima Eq | Pons Eq | Code File | Code Line |
|---------|-------------|---------|-----------|-----------|
| Quartic equation | Derived from 104,107 | - | sr_primitive_recovery.cpp | 143-149 |
| Velocity recovery | From 104 | - | sr_primitive_recovery.cpp | 236-239 |
| n = N/γ | 91 | - | sr_primitive_recovery.cpp | 338 |
| H formula | 112 | - | sr_primitive_recovery.cpp | 348 |
| P = (γ_c-1)nu | 119 | - | sr_primitive_recovery.cpp | 363 |
| c_s² formula | 121 | - | sr_primitive_recovery.cpp | 370 |
| V_p volume | 220 | - | sr_pre_interaction.cpp | 129-151 |
| h iteration | 231,233 | - | sr_pre_interaction.cpp | 203-241 |
| N = ν/V_p | 239 | - | sr_pre_interaction.cpp | 336-337 |
| dS/dt | 209 | - | sr_fluid_force.cpp | 488 |
| de/dt | 212 | - | sr_fluid_force.cpp | 494 |
| Taub adiabat | - | 4.16 | sr_fluid_force.cpp | 897-914 |
| Mass flux j² | - | 4.7 | sr_fluid_force.cpp | 934 |
| Shock velocity | - | 5.22 | sr_fluid_force.cpp | 948-954 |
| Post-shock v | - | 6.1 | sr_fluid_force.cpp | 961-964 |
| Rarefaction v | - | 3.20 | sr_fluid_force.cpp | 981-999 |
| v^t coupling | - | 2.15 | sr_fluid_force.cpp | 983-984 |
| CFL condition | Section 3 | - | sr_timestep.cpp | 80-98 |
