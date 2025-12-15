# Theoretical Derivation: IMBH-Molecular Cloud Scattering Simulation

## Executive Summary

This document provides first-principles derivations for simulating the tidal deformation and shock formation during intermediate-mass black hole (IMBH) molecular cloud scattering, as observed in CO-0.40-0.22 (Oka et al. 2017).

**Key Conclusions:**
1. **MHD effects: NOT required** — super-Alfvénic flow ($M_A \sim 250$)
2. **Two-fluid approximation: NOT required** — ambipolar diffusion timescale $\gg$ dynamical timescale
3. **Shock formation: YES** — Mach number $\mathcal{M} \sim 10$–15 during tidal compression
4. **Chemistry: Use equilibrium** — chemical timescales $\lesssim$ dynamical timescale

**Minimum Simulation Requirements:**
- Single-fluid hydrodynamics with self-gravity and IMBH point mass
- Shock-capturing artificial viscosity
- $N \geq 10^6$ particles for $\beta \approx 3$

---

## 1. Physical Parameters from Observations

### 1.1 IMBH Properties (Oka et al. 2017)

| Parameter | Value | Symbol |
|-----------|-------|--------|
| Black hole mass | $10^5\,M_\odot$ | $M_\mathrm{BH}$ |
| Distance from Galactic center | $\sim 60$ pc | — |
| Point source size upper limit | $< 0.022$ pc | — |

### 1.2 Molecular Cloud Properties (CO-0.40-0.22)

| Parameter | Value | Symbol |
|-----------|-------|--------|
| Cloud size | $\sim 5$ pc | $R_\mathrm{cloud}$ |
| Dense clump size | $\sim 0.3$ pc | $r_\mathrm{clump}$ |
| Velocity width | $\sim 100$ km/s | $\Delta v$ |
| Velocity dispersion | 22 km/s | $\sigma_v$ |
| Dense clump mass | $\sim 40\,M_\odot$ | $M_\mathrm{clump}$ |
| Virial mass | $\sim 4100\,M_\odot$ | $M_\mathrm{vir}$ |
| Temperature | 60 K | $T$ |
| Number density (dense) | $10^{6.5}$ cm$^{-3}$ | $n_\mathrm{dense}$ |

### 1.3 Cloud Mass Estimates

| Estimate | Value | Method | Use |
|----------|-------|--------|-----|
| Dense clump (HCN) | 40 $M_\odot$ | Non-LTE excitation | Dense core only |
| Virial mass | 4100 $M_\odot$ | $5\sigma_v^2 R/G$ | Upper limit |
| Tidal mass | 800 $M_\odot$ | $M_\mathrm{BH}(r/r_\mathrm{peri})^3$ | Consistency check |
| **N-body model** | **1000 $M_\odot$** | Best fit to P-V | **Recommended** |

**Physical interpretation:** $M_\mathrm{vir} \gg M_\mathrm{clump}$ confirms the clump is **not gravitationally bound** — the large velocity dispersion is due to tidal acceleration by the IMBH.

### 1.4 Key Observational Diagnostics

The characteristic "parallelogram" shape in the P-V diagram indicates tidal stretching:
- Velocity range: $V_\mathrm{LSR} = -105$ to $-5$ km/s ($\Delta v \sim 100$ km/s)
- Dense clump offset: 0.2 pc from continuum source
- Spectral index $\alpha = 1.18$ (inverted spectrum → compact accretion source)

---

## 2. Physics Requirements Analysis

### 2.1 MHD Effects: NOT Required

**Plasma beta** (thermal/magnetic pressure ratio):
$$\beta_\mathrm{bulk} \approx 1.2, \quad \beta_\mathrm{dense} \approx 1.5$$

**Alfvén Mach number:**
$$M_A = \frac{v_\mathrm{bulk}}{v_A} = \frac{100\,\mathrm{km/s}}{0.39\,\mathrm{km/s}} \approx 250$$

**Force ratio:**
$$\frac{F_\mathrm{tidal}}{F_\mathrm{magnetic}} \approx 300$$

**Conclusion:** Flow is highly super-Alfvénic. Kinetic energy dominates magnetic energy by $M_A^2 \sim 6 \times 10^4$. MHD is unnecessary for primary dynamics.

### 2.2 Two-Fluid Effects: NOT Required

**Ionization fraction:** $\chi_i \approx 4 \times 10^{-8}$

**Timescales:**
| Timescale | Value | Ratio to $\tau_\mathrm{dyn}$ |
|-----------|-------|------------------------------|
| Neutral-ion coupling $\tau_{ni}$ | 5000 yr | 0.5 |
| Ambipolar diffusion $\tau_\mathrm{AD}$ | $10^7$ yr | 1000 |
| Dynamical $\tau_\mathrm{dyn}$ | $10^4$ yr | 1 |

**Conclusion:** Neutrals and ions are well-coupled ($\tau_{ni} < \tau_\mathrm{dyn}$). Ambipolar diffusion is negligible. Single-fluid hydrodynamics is sufficient.

### 2.3 Shock Formation: YES

**Tidal compression Mach number:**
$$\mathcal{M}_\mathrm{compress} = \frac{v_\mathrm{compress}}{c_s} \approx \frac{7\,\mathrm{km/s}}{0.47\,\mathrm{km/s}} \approx 15$$

**Post-shock temperature** (before cooling):
$$T_\mathrm{post} \approx 2600\,\mathrm{K}$$

**Conclusion:** Strong shocks ($\mathcal{M} \sim 10$–15) will form. Simulation must use shock-capturing scheme.

---

## 3. Tidal Compression Theory (Coughlin & Nixon 2021)

### 3.1 The Penetration Parameter β

The key dimensionless parameter controlling tidal compression strength:
$$\beta \equiv \frac{r_t}{r_p}$$

where $r_t = R_\mathrm{cloud}(M_\mathrm{BH}/M_\mathrm{cloud})^{1/3}$ is the tidal radius and $r_p$ is the pericentre distance.

**For CAT_OKA:** $\beta = 5.25/1.69 \approx 3.1$

### 3.2 Homologous Compression

In the pressureless limit, all fluid elements compress homologously:
$$\frac{z}{z_0} = H(\tau)$$

where $\tau$ is dimensionless time related to orbital phase.

**With pressure and self-gravity**, $H(\tau)$ satisfies:
$$\boxed{\mathcal{L}[H] - \frac{2}{\beta^3}\frac{\rho_c}{\rho_\star}\left(H^{-\gamma} - 1\right)\cosh^6(\tau) = 0}$$

**Key result:** The black hole mass does not explicitly appear — the physics is determined by $\beta$ and $\rho_c/\rho_\star$ alone.

### 3.3 Shock Formation Regimes

| $\beta$ Range | Physics | $\rho_\mathrm{max}$ Scaling |
|---------------|---------|----------------------------|
| $\beta < 1$ | Partial disruption | Cloud survives |
| $1 < \beta < 3$ | Adiabatic compression | $\propto \beta^3$ (reduced coeff.) |
| $3 < \beta < 10$ | Shock forms, stalls | $\propto \beta^3$ (reduced coeff.) |
| $\beta > 10$ | Shock-dominated | $\propto \beta^{1.62}$ |

**For CAT_OKA ($\beta \approx 3.1$):** At the threshold of shock formation. Compression primarily adiabatic but shocks may form in outer layers.

### 3.4 Key Insight: Weak Shocks

Coughlin & Nixon found that **internal shocks are weak** ($\mathcal{M} \sim 1.2$), even for large $\beta$. This is because pre-shock gas is simultaneously compressing, raising its sound speed.

---

## 4. Simulation Configuration

### 4.1 Required Physics

| Physics | Include? | Justification |
|---------|----------|---------------|
| **Hydrodynamics** | YES | Primary dynamics |
| **Self-gravity** | YES | Cloud structure |
| **Point-mass gravity** | YES | IMBH tidal field |
| **Shock capturing** | YES | $\mathcal{M} \sim 10$–15 |
| **Sink accretion** | YES | Remove particles that fall into BH |
| MHD | OPTIONAL | Not dynamically important |
| Two-fluid | NO | $\tau_\mathrm{AD} \gg \tau_\mathrm{dyn}$ |
| Radiative cooling | OPTIONAL | For post-shock structure |

### 4.1.1 Fixed Black Hole Implementation

The IMBH is implemented as a **fixed point mass** at the origin with Plummer softening:

**BH Parameters:**
| Parameter | Symbol | Typical Value | Description |
|-----------|--------|---------------|-------------|
| BH mass | $M_\mathrm{BH}$ | $10^5\,M_\odot$ | Fixed, does not change on accretion |
| BH position | $\vec{r}_\mathrm{BH}$ | $(0,0,0)$ | Fixed at origin |
| Softening | $\varepsilon$ | 0.001 pc (~200 AU) | Prevents singularity |
| Sink radius | $r_\mathrm{sink}$ | 0.005 pc (~1000 AU) | Accretion boundary |

**Plummer-Softened Gravity:**
$$\vec{a}_\mathrm{BH}(\vec{r}) = -\frac{G M_\mathrm{BH} \vec{r}}{(r^2 + \varepsilon^2)^{3/2}}$$

$$\Phi_\mathrm{BH}(r) = -\frac{G M_\mathrm{BH}}{\sqrt{r^2 + \varepsilon^2}}$$

**Sink Accretion Algorithm:**
For each particle inside $r < r_\mathrm{sink}$:
1. **Inflow check:** $v_\mathrm{rad} = \vec{v} \cdot \hat{r} < 0$ (moving toward BH)
2. **Bound check:** $E = \frac{1}{2}|\vec{v}|^2 + \Phi_\mathrm{BH} < 0$ (gravitationally bound)

If both conditions satisfied, particle is removed. BH mass remains fixed (fixed-potential approximation).

**Time Stepping:**
Two constraints must be satisfied:
1. **Hydro CFL:** $\Delta t_\mathrm{hydro} \leq C_\mathrm{CFL} \cdot h / v_\mathrm{sig}$
2. **Gravity limiter:** $\Delta t_\mathrm{grav} \leq \eta_t \cdot \sqrt{h / |\vec{a}|}$

Final timestep: $\Delta t = \min(\Delta t_\mathrm{hydro}, \Delta t_\mathrm{grav})$

Typical safe values: $C_\mathrm{CFL} \sim 0.2$–$0.4$, $\eta_t \sim 0.1$–$0.2$

### 4.2 Initial Conditions (Oka et al. 2017 N-body Model)

**Exact parameters from Oka et al. Methods section:**

| Parameter | Value |
|-----------|-------|
| BH mass $M_\mathrm{BH}$ | $10^5\,M_\odot$ |
| Cloud mass | $1000\,M_\odot$ |
| Initial position $(X_0, Y_0)$ | $(9.8, -0.65)$ pc |
| Initial velocity $(v_X, v_Y)$ | $(-8.19, 0.4)$ km/s |
| Cloud dispersion $\sigma_r$ | 0.2 pc |
| Velocity dispersion | 1.43 km/s |
| Orbital inclination $i$ | 70° |
| Position angle PA | 41.6° |
| $V_\mathrm{LSR}$ | −120 km/s |
| Best-fit snapshot | $t = 7.2 \times 10^5$ yr |

**Note:** For SPH simulations with resolved cloud ($R_\mathrm{cloud} = 1.13$ pc), pericentre must satisfy $r_\mathrm{peri} \geq 1.5 R_\mathrm{cloud}$ to avoid numerical issues.

### 4.3 Resolution Requirements

#### From Coughlin & Nixon (2021):

Under-resolved simulations systematically **underestimate compression**.

| $\beta$ | Min Particles | Recommended |
|---------|---------------|-------------|
| 1 | $10^4$ | $10^5$ |
| 2 | $10^5$ | $10^6$ |
| 3 | $5 \times 10^5$ | $2 \times 10^6$ |
| 4 | $2 \times 10^6$ | $10^7$ |
| 8 | $2 \times 10^7$ | $10^8$ |
| 16 | $10^8$ | $10^9$ |

**For CAT_OKA ($\beta \approx 3$):**
$$\boxed{N_\mathrm{min} = 10^6, \quad N_\mathrm{rec} = 10^7}$$

#### Physical Constraints:

| Scale | Value | Required $h$ |
|-------|-------|--------------|
| Jeans length (dense) | 0.024 pc | < 0.006 pc |
| Jeans length (post-shock) | 0.016 pc | **< 0.004 pc** |
| In-plane tidal | 0.015 pc | **< 0.004 pc** |
| Vertical tidal | 0.022 pc | < 0.005 pc |

**Critical resolution:** $h_\mathrm{max} = 0.004$ pc

### 4.4 Numerical Safety Checks

**Sink radius vs smoothing length:**
- Rule of thumb: $r_\mathrm{sink} \gtrsim h_\mathrm{min}$
- If $r_\mathrm{sink} \gg h_\mathrm{min}$: sink removes resolved structure
- If $r_\mathrm{sink} \ll h_\mathrm{min}$: sink may not help with timestepping

**Softening vs sink radius:**
- Requirement: $\varepsilon < r_\mathrm{sink}$ (softening inside sink)
- Typical: $\varepsilon \approx r_\mathrm{sink}/5$

**JSON Configuration Example:**
```json
{
  "imbh_parameters": {
    "enabled": true,
    "M_BH": 1e5,
    "BH_initial_position": [0.0, 0.0, 0.0],
    "BH_initial_velocity": [0.0, 0.0, 0.0],
    "softening_epsilon": 0.001,
    "sink_radius": 0.005,
    "enable_sink": true
  }
}
```

### 4.5 Temporal Requirements

| Constraint | Value |
|------------|-------|
| CFL (bulk velocity) | $\Delta t \lesssim 10$ yr |
| Orbital accuracy | $\Delta t < 100$ yr |
| Total simulation time | $\gtrsim 10^6$ yr |
| Timesteps | $\sim 10^5$ |

---

## 5. Post-Processing for Observational Comparison

### 5.1 Chemistry Treatment

**Recommendation:** Equilibrium chemistry with post-processing radiative transfer.

**Justification:**
- CO: $\tau_\mathrm{form} \ll \tau_\mathrm{dyn}$ in dense gas → equilibrium
- HCN: $\tau_\mathrm{form} \sim \tau_\mathrm{dyn}$ → approximate equilibrium
- PV structure is kinematically determined; chemistry affects intensity, not morphology

**Equilibrium abundances:**
- CO: $X(\mathrm{CO}) = 10^{-4}$ (constant in well-shielded regions)
- HCN: $X(\mathrm{HCN}) = 10^{-8} \times (n/10^4)^{0.5}$

### 5.2 P-V Diagram Generation

**Step 1: Transform to observer frame**
```
Rotation: R = R_z(PA=41.6°) · R_x(i=70°)
v_LOS = v · ẑ_obs + V_LSR
```

**Step 2: Compute line-of-sight velocity**
$$v_\mathrm{LOS,i} = \mathbf{v}_i \cdot \hat{z}_\mathrm{obs} + V_\mathrm{LSR}$$

**Step 3: Grid and weight by emissivity**
- CO J=2-1: LTE (densities exceed $n_\mathrm{crit} \sim 10^3$ cm$^{-3}$)
- HCN J=3-2: Non-LTE with RADEX (densities marginal at $n_\mathrm{crit} \sim 3 \times 10^6$ cm$^{-3}$)

**Step 4: Convolve with ALMA beam**
- CO: $1.87'' \times 1.14''$
- HCN: $1.52'' \times 0.60''$

### 5.3 Key Features to Match

| Feature | Observation | Target |
|---------|-------------|--------|
| Velocity range | −105 to −5 km/s | $\sim 100$ km/s width |
| P-V gradient | Parallelogram | Tidal stretching |
| Dense clump offset | 0.2 pc | Track density max |
| Velocity at clump | $V_\mathrm{LSR} \sim -60$ km/s | Main body |

---

## 6. Verification Workflow

### 6.1 Simulation vs Theory Comparisons

| Plot | X-axis | Y-axis | Theory Overlay |
|------|--------|--------|----------------|
| Lagrangian height | $t/t_\mathrm{dyn}$ | $z/z_0$ | $H(\tau)$ solution |
| Central density | $t/t_\mathrm{dyn}$ | $\rho_c/\rho_{c,0}$ | $1/H(\tau)$ |
| Velocity profile | $z$ | $v_z$ | Linear (homologous) |
| β-scaling | $\beta$ | $\rho_\mathrm{max}$ | Power-law fits |

### 6.2 Convergence Test Protocol

1. Run at resolutions $N$, $2N$, $4N$
2. Compare $\rho_c(t)$ evolution
3. Convergence criterion:
   $$\frac{|\rho_\mathrm{max}^{(2N)} - \rho_\mathrm{max}^{(N)}|}{\rho_\mathrm{max}^{(2N)}} < 0.1$$

### 6.3 Shock Detection Methods

1. **Velocity divergence:** $\nabla \cdot \vec{v} < 0$ (strong negative → shock)
2. **Artificial viscosity:** Track $\Pi_{ij}$ (high values → shock heating)
3. **Entropy:** Non-adiabatic increase indicates shock passage

---

## 7. Conclusions

### Required Physics
- Single-fluid hydrodynamics (SPH)
- Self-gravity
- Point-mass external gravity (IMBH)
- Shock-capturing artificial viscosity

### Optional Enhancements
- Radiative cooling (for post-shock structure)
- MHD (for small-scale structure, not primary dynamics)

### NOT Needed
- Two-fluid approximation
- Ambipolar diffusion
- Time-dependent chemistry

### Key Numbers

| Quantity | Value |
|----------|-------|
| MHD importance | $M_A \sim 250$ (negligible) |
| Two-fluid importance | $\tau_\mathrm{AD}/\tau_\mathrm{dyn} \sim 1000$ (negligible) |
| Shock Mach number | $\mathcal{M} \sim 10$–15 (strong shocks) |
| Penetration parameter | $\beta \approx 3$ (shock threshold) |
| Minimum particles | $10^6$ |
| Recommended particles | $10^7$ |

---

## References

1. Oka, T., et al. (2017). "Millimetre-wave Emission from an Intermediate-Mass Black Hole Candidate in the Milky Way." *Nature Astronomy*.

2. Coughlin, E. R., & Nixon, C. J. (2021). "Stars Crushed by Black Holes. II. A Physical Model of Adiabatic Compression and Shock Formation in Tidal Disruption Events." *ApJ*.

3. Inoue, T., & Inutsuka, S. (2008). "Two-Fluid MHD Simulations of Converging HI Flows in the Interstellar Medium." *ApJ*.

4. Koyama, H., & Inutsuka, S. (2000). "Molecular Cloud Formation in Shock-compressed Layers." *ApJ*, 532, 980.

5. van der Tak, F. F. S., et al. (2007). "RADEX: A Computer Program for Fast Non-LTE Analysis." *A&A*, 468, 627.

6. Dullemond, C. P., et al. (2012). "RADMC-3D: A Multi-purpose Radiative Transfer Tool." *ASCL*, 1202.015.

---

## Appendix A: Key Formulas

### A.1 Tidal Radius
$$r_t = \left(\frac{M_\mathrm{cloud}}{3M_\mathrm{BH}}\right)^{1/3} d$$

### A.2 Alfvén Speed
$$v_A = \frac{B}{\sqrt{4\pi\rho}} \approx 0.39\,\mathrm{km/s}$$

### A.3 Jeans Length
$$\lambda_J = c_s \sqrt{\frac{\pi}{G\rho}}$$

### A.4 Ambipolar Diffusion Timescale
$$\tau_\mathrm{AD} = \frac{L^2}{v_A^2 \tau_{ni}}$$

### A.5 Homologous Compression
$$\frac{z}{z_0} = H(\tau), \quad \frac{\rho_c}{\rho_{c,0}} = \frac{1}{H(\tau)}$$

---

## Appendix B: Physical Constants

| Constant | Value |
|----------|-------|
| $G$ | $6.67 \times 10^{-8}$ cm³ g⁻¹ s⁻² |
| $k_B$ | $1.38 \times 10^{-16}$ erg K⁻¹ |
| $m_H$ | $1.67 \times 10^{-24}$ g |
| $M_\odot$ | $2 \times 10^{33}$ g |
| 1 pc | $3.09 \times 10^{18}$ cm |
| 1 km/s | $10^5$ cm/s |
| 1 yr | $3.15 \times 10^7$ s |
