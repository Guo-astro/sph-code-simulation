# Theoretical Derivation: IMBH-Molecular Cloud Scattering Simulation

## Executive Summary

This document provides first-principles derivations for simulating tidal deformation and shock formation during intermediate-mass black hole (IMBH) molecular cloud scattering, as observed in CO-0.40-0.22 (Oka et al. 2017).

**Key Conclusions:**
1. **MHD effects: NOT required** — super-Alfvénic flow ($M_A \sim 250$)
2. **Two-fluid approximation: NOT required** — $\tau_\mathrm{AD} \gg \tau_\mathrm{dyn}$
3. **Shock formation: YES** — Mach number $\mathcal{M} \sim 10$–15
4. **Scientific goal:** Test whether IMBH hypothesis survives rigorous hydrodynamical treatment

**Minimum Requirements:**
- Single-fluid SPH with self-gravity and IMBH point mass
- Shock-capturing artificial viscosity
- $N \geq 10^6$ particles for $\beta \approx 3$

---

## 1. Observational Constraints

### 1.1 System Parameters

| Parameter | IMBH | Cloud (CO-0.40-0.22) |
|-----------|------|----------------------|
| Mass | $10^5\,M_\odot$ | 40 $M_\odot$ (HCN), 1000 $M_\odot$ (model) |
| Size | $< 0.022$ pc | 0.3 pc (dense), 5 pc (total) |
| Velocity | — | 100 km/s width, $\sigma_v = 22$ km/s |
| Density | — | $10^{6.5}$ cm$^{-3}$ (dense clump) |
| Temperature | — | 60 K |
| Distance from GC | ~60 pc | — |

### 1.2 Mass Derivation from Radio Observations

**Dense clump mass from HCN J=3-2 (Non-LTE):**

Critical density: $n_\mathrm{crit} = A_{ul}/\gamma_{ul} \approx 3 \times 10^6$ cm$^{-3}$

Column density from integrated intensity:
$$N(\mathrm{HCN}) = \frac{8\pi k_B \nu^2}{h c^3 A_{ul}} \frac{Q(T_\mathrm{ex})}{g_u} \frac{e^{E_u/k_B T_\mathrm{ex}}}{1 - e^{-h\nu/k_B T_\mathrm{ex}}} \int T_b\,dv$$

Mass: $M = \mu m_H \cdot N(\mathrm{HCN})/X(\mathrm{HCN}) \cdot \Omega D^2$

With $T_k = 60$ K, $n = 10^{6.5}$ cm$^{-3}$, $X(\mathrm{HCN}) \approx 10^{-8}$: **$M_\mathrm{clump} \approx 40\,M_\odot$**

**Virial mass:**
$$M_\mathrm{vir} = \frac{5\sigma_V^2 R}{G} = \frac{5 \times (22)^2 \times 0.15}{4.30 \times 10^{-3}} \approx 4100\,M_\odot$$

### 1.3 Critical Insight: Cloud is NOT Gravitationally Bound

$$\frac{M_\mathrm{vir}}{M_\mathrm{gas}} = \frac{4100}{40} \approx 100$$

This proves:
1. Velocity dispersion is externally imposed (tidal acceleration)
2. Cloud is being tidally disrupted
3. Observed structure is POST-encounter, not equilibrium

---

## 2. Why the Cloud is NOT in Thermal Equilibrium

### 2.1 Jeans Length Analysis

$$\lambda_J = c_s \sqrt{\frac{\pi}{G\rho}} \approx 0.4 \times \sqrt{\frac{T}{n}}\,\mathrm{pc}$$

| Phase | $n$ (cm$^{-3}$) | $T$ (K) | $\lambda_J$ (pc) |
|-------|-----------------|---------|------------------|
| WNM | 0.5 | 6000 | 44 |
| CNM | 50 | 80 | 0.5 |
| Molecular | 1000 | 10 | 0.04 |
| Dense | $10^{6.5}$ | 8 | 0.003 |

**Problem:** Observed clump (0.3 pc) is 100× larger than equilibrium Jeans length at observed density.

### 2.2 Three Independent Proofs

1. **Virial ratio:** $M_\mathrm{vir}/M_\mathrm{gas} \approx 100$ (should be ~1 for equilibrium)
2. **Size:** Observed 0.3 pc >> equilibrium $\lambda_J \sim 0.003$ pc
3. **Velocity:** 100 km/s >> sound speed 0.3 km/s (supersonic, non-thermal)

### 2.3 Why K&I 2000 Equilibrium Produces Small Clouds

Bonnor-Ebert mass scaling: $M_\mathrm{BE} \propto T^2/\sqrt{P}$

For cold molecular gas ($T \sim 10$ K) vs warm gas ($T \sim 6000$ K):
$$\frac{M_\mathrm{BE,cold}}{M_\mathrm{BE,warm}} \sim \left(\frac{10}{6000}\right)^2 \sim 3 \times 10^{-6}$$

Result: K&I equilibrium clouds are limited to $R \sim 0.05-0.25$ pc for molecular densities.

---

## 3. SPH vs N-body: Scientific Purpose

### 3.1 What Oka's N-body Cannot Capture

| Physics | N-body | SPH |
|---------|--------|-----|
| Pressure forces | ❌ | ✅ |
| Shock formation | ❌ | ✅ |
| Thermal physics | ❌ | ✅ |
| Compression heating | ❌ | ✅ |
| Energy dissipation | ❌ | ✅ |

### 3.2 Possible Outcomes

| Outcome | Implication |
|---------|-------------|
| SPH matches observations | IMBH hypothesis validated with rigorous physics |
| SPH differs from observations | IMBH parameters need adjustment OR alternative needed |
| SPH reveals new physics | New predictions for follow-up observations |

### 3.3 Simulation Strategy

1. Use physically motivated IC (Lane-Emden + K&I temperatures)
2. Vary IMBH parameters: $M = 10^4 - 10^6\,M_\odot$, $Y_0 = 0.5 - 3$ pc
3. Compare to **observations**, not Oka's N-body
4. Report: Confirmed / Refined / Challenged

---

## 4. Physics Requirements

### 4.1 MHD: NOT Required

- Alfvén Mach number: $M_A = v_\mathrm{bulk}/v_A = 100/0.39 \approx 250$
- Force ratio: $F_\mathrm{tidal}/F_\mathrm{magnetic} \approx 300$
- Kinetic energy dominates by $M_A^2 \sim 6 \times 10^4$

### 4.2 Two-Fluid: NOT Required

| Timescale | Value | Ratio to $\tau_\mathrm{dyn}$ |
|-----------|-------|------------------------------|
| Neutral-ion coupling | 5000 yr | 0.5 |
| Ambipolar diffusion | $10^7$ yr | 1000 |
| Dynamical | $10^4$ yr | 1 |

Single-fluid hydrodynamics is sufficient.

### 4.3 Shock Formation: YES

- Compression Mach number: $\mathcal{M} \approx v_\mathrm{compress}/c_s = 7/0.47 \approx 15$
- Post-shock temperature: $T_\mathrm{post} \approx 2600$ K (before cooling)

---

## 5. Tidal Compression Theory

### 5.1 Penetration Parameter

$$\beta \equiv \frac{r_t}{r_p}, \quad r_t = R_\mathrm{cloud}\left(\frac{M_\mathrm{BH}}{M_\mathrm{cloud}}\right)^{1/3}$$

For CO-0.40-0.22: $\beta \approx 3.1$ (shock formation threshold)

### 5.2 Compression Regimes (Coughlin & Nixon 2021)

| $\beta$ | Physics | Scaling |
|---------|---------|---------|
| $< 1$ | Partial disruption | Survives |
| $1-3$ | Adiabatic compression | $\rho_\mathrm{max} \propto \beta^3$ |
| $3-10$ | Shock forms, stalls | $\rho_\mathrm{max} \propto \beta^3$ (reduced) |
| $> 10$ | Shock-dominated | $\rho_\mathrm{max} \propto \beta^{1.62}$ |

### 5.3 Homologous Compression

With pressure and self-gravity:
$$\mathcal{L}[H] - \frac{2}{\beta^3}\frac{\rho_c}{\rho_\star}\left(H^{-\gamma} - 1\right)\cosh^6(\tau) = 0$$

Key result: Physics determined by $\beta$ and $\rho_c/\rho_\star$ alone.

---

## 6. Simulation Configuration

### 6.1 Required Physics

| Physics | Include | Justification |
|---------|---------|---------------|
| Hydrodynamics | YES | Primary dynamics |
| Self-gravity | YES | Cloud structure |
| Point-mass IMBH | YES | Tidal field |
| Shock capturing | YES | $\mathcal{M} \sim 15$ |
| Sink accretion | YES | Remove accreted particles |
| Radiative cooling | OPTIONAL | Post-shock structure |
| MHD | NO | $M_A \sim 250$ |
| Two-fluid | NO | $\tau_\mathrm{AD} \gg \tau_\mathrm{dyn}$ |

### 6.2 IMBH Implementation

Plummer-softened gravity:
$$\vec{a}_\mathrm{BH} = -\frac{G M_\mathrm{BH} \vec{r}}{(r^2 + \varepsilon^2)^{3/2}}$$

| Parameter | Value |
|-----------|-------|
| Softening $\varepsilon$ | 0.001 pc |
| Sink radius | 0.005 pc |
| Sink criteria | $v_r < 0$ AND $E < 0$ |

### 6.3 Oka et al. Orbital Parameters

| Parameter | Value |
|-----------|-------|
| Initial position | $(X_0, Y_0) = (9.8, -0.65)$ pc |
| Initial velocity | $(v_X, v_Y) = (-8.19, 0.4)$ km/s |
| Cloud mass | 1000 $M_\odot$ |
| Cloud $\sigma_r$ | 0.2 pc |
| Velocity dispersion | 1.43 km/s |
| Inclination | 70° |
| Position angle | 41.6° |
| $V_\mathrm{LSR}$ | −120 km/s |
| Best-fit time | $7.2 \times 10^5$ yr |

### 6.4 Resolution Requirements

| $\beta$ | Min Particles | Recommended |
|---------|---------------|-------------|
| 1 | $10^4$ | $10^5$ |
| 3 | $5 \times 10^5$ | $2 \times 10^6$ |
| 8 | $2 \times 10^7$ | $10^8$ |

For $\beta \approx 3$: $N_\mathrm{min} = 10^6$, $N_\mathrm{rec} = 10^7$

Critical smoothing length: $h_\mathrm{max} = 0.004$ pc

---

## 7. Initial Condition Design

### 7.1 IC Options Summary

| IC Type | Use Case | Pros | Cons |
|---------|----------|------|------|
| **Lane-Emden + K&I** | **RECOMMENDED** | Smooth profile, thermal physics | Initial transient |
| Gaussian Virialized | Oka comparison | Matches N-body | Not hydrostatic |
| K&I Bonnor-Ebert | Phase structure | True equilibrium | Too small (~0.2 pc) |

### 7.2 Recommended: Lane-Emden + K&I Temperatures

- Density: Lane-Emden n=3/2 polytrope
- Temperature: $T(r) = T_\mathrm{eq}(n(r))$ from K&I 2000 tables
- Represents realistic pre-encounter molecular cloud
- NOT in exact equilibrium — will adjust during initial phase

### 7.3 Gaussian IC (For Oka Comparison Only)

Density: $\rho(r) = \rho_c \exp(-r^2/2\sigma_r^2)$

Mass-radius: $M = (2\pi)^{3/2} \rho_c \sigma_r^3$

Virial velocity: $\sigma_v = \sqrt{GM/\alpha\sigma_r}$ where $\alpha \approx 2.4$

With $M = 1000\,M_\odot$, $\sigma_r = 0.2$ pc: $n_c \approx 2000$ cm$^{-3}$, $\sigma_v \approx 1.4$ km/s

### 7.4 Bonnor-Ebert with K&I 2000 EOS

Barotropic EOS: $P(\rho) = n k_B T_\mathrm{eq}(n)$

Effective sound speed: $c_\mathrm{eff}^2 = (k_B T/\mu m_H)(1 + d\ln T/d\ln n)$

Hydrostatic ODE:
$$\frac{d\rho}{dr} = -\rho\frac{GM(r)}{r^2 c_\mathrm{eff}^2}, \quad \frac{dM}{dr} = 4\pi r^2 \rho$$

Integrate until $P = P_\mathrm{ext}$. Tune $\rho_c$ and $P_\mathrm{ext}$ to hit target $R$ and $M$.

**Limitation:** Cold phase ($T \sim 10$ K) produces small clouds ($R \sim 0.05-0.25$ pc).

---

## 8. Post-Processing

### 8.1 P-V Diagram Generation

1. **Transform to observer frame:** $R = R_x(i=70°) \cdot R_z(PA=41.6°)$
2. **Line-of-sight velocity:** $v_\mathrm{LOS} = \vec{v} \cdot \hat{z}_\mathrm{obs} + V_\mathrm{LSR}$
3. **Grid and weight** by mass or emissivity
4. **Convolve** with ALMA beam

### 8.2 Key Features to Match

| Feature | Observation | Target |
|---------|-------------|--------|
| Velocity range | −105 to −5 km/s | ~100 km/s width |
| P-V shape | Parallelogram | Tidal stretching |
| Dense clump offset | 0.2 pc | Density maximum |
| Clump velocity | ~−60 km/s | Main body |

### 8.3 Chemistry

- CO: $X(\mathrm{CO}) = 10^{-4}$ (equilibrium in shielded gas)
- HCN: $X(\mathrm{HCN}) = 10^{-8} \times (n/10^4)^{0.5}$
- Use RADEX for non-LTE emission modeling

---

## 9. Verification

### 9.1 Convergence Test

Run at $N$, $2N$, $4N$ particles. Criterion:
$$\frac{|\rho_\mathrm{max}^{(2N)} - \rho_\mathrm{max}^{(N)}|}{\rho_\mathrm{max}^{(2N)}} < 0.1$$

### 9.2 Shock Detection

1. Velocity divergence: $\nabla \cdot \vec{v} < 0$
2. Artificial viscosity: High $\Pi_{ij}$
3. Entropy: Non-adiabatic increase

### 9.3 Theory Comparison

| Plot | Theory |
|------|--------|
| $z/z_0$ vs $t$ | $H(\tau)$ solution |
| $\rho_c/\rho_{c,0}$ vs $t$ | $1/H(\tau)$ |
| $v_z$ vs $z$ | Linear (homologous) |

---

## 10. Key Formulas

| Formula | Expression |
|---------|------------|
| Tidal radius | $r_t = R_\mathrm{cloud}(M_\mathrm{BH}/M_\mathrm{cloud})^{1/3}$ |
| Virial mass | $M_\mathrm{vir} = 5\sigma_v^2 R/G$ |
| Jeans length | $\lambda_J = c_s\sqrt{\pi/G\rho}$ |
| Alfvén speed | $v_A = B/\sqrt{4\pi\rho}$ |
| Plummer potential | $\Phi = -GM/\sqrt{r^2 + \varepsilon^2}$ |

---

## 11. Physical Constants

| Constant | Value |
|----------|-------|
| $G$ | $6.67 \times 10^{-8}$ cm³ g⁻¹ s⁻² |
| $k_B$ | $1.38 \times 10^{-16}$ erg K⁻¹ |
| $m_H$ | $1.67 \times 10^{-24}$ g |
| $M_\odot$ | $2 \times 10^{33}$ g |
| 1 pc | $3.09 \times 10^{18}$ cm |
| 1 km/s | $10^5$ cm/s |

---

## References

1. Oka, T., et al. (2017). Nature Astronomy. (IMBH candidate in CO-0.40-0.22)
2. Coughlin, E. R., & Nixon, C. J. (2021). ApJ. (Tidal compression theory)
3. Koyama, H., & Inutsuka, S. (2000). ApJ 532, 980. (ISM cooling)
4. Inoue, T., & Inutsuka, S. (2008). ApJ. (Two-fluid MHD)
