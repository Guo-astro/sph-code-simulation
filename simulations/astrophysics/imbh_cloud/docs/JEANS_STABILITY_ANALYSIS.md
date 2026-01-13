# Jeans Stability Analysis for ISM Cooling Simulations

## Problem Statement

Phase 2.5 ISM cooling simulation experienced timestep collapse from ~0.01 to ~10⁻⁶, with the Newton-Raphson smoothing length solver failing for many particles. Investigation revealed the cloud was undergoing gravitational collapse.

## Root Cause Analysis

### Original Configuration (UNSTABLE)
| Parameter | Value |
|-----------|-------|
| T_cloud | 10 K |
| n_center | 1800 cm⁻³ |
| M_cloud | 40 M_sun |
| R_cloud | 0.75 pc |
| xi_s | 6.0 (near-critical BE) |

### Jeans Stability Check (Original)
- **λ_J** = 0.77 pc at T=10K, n=1800
- **M_J** = 23 M_sun
- **R/λ_J** = 0.97 (marginal!)
- **M/M_J** = 1.75 (**UNSTABLE!**)

### What Went Wrong

1. **The cloud mass exceeded the Jeans mass** (M/M_J = 1.75 > 1)
2. **K&I 2000 cooling reduced temperature** from 10K → 7K (equilibrium at n=1800)
3. **Lower temperature reduced pressure support**, triggering collapse
4. **Free-fall time** t_ff = 1.08 Myr matched observed collapse time (~2 Myr ≈ 2×t_ff)

## Key Physics Insight

For a Bonnor-Ebert sphere, the dimensionless truncation radius xi_s determines stability:

$$\frac{R}{\lambda_J} = \frac{\xi_s \cdot r_0}{c_s \sqrt{\pi/(G\rho_c)}} = \frac{\xi_s}{\sqrt{4\pi}} \approx \frac{\xi_s}{3.54}$$

This means:
| xi_s | R/λ_J | Status |
|------|-------|--------|
| 6.45 (critical) | 1.82 | **ALWAYS UNSTABLE** |
| 6.0 (standard) | 1.70 | **UNSTABLE** |
| 3.54 | 1.00 | Marginal |
| 3.0 | 0.85 | **STABLE** |

**Critical insight**: A critical Bonnor-Ebert sphere (xi_s ≈ 6.45) is **inherently Jeans unstable by construction**!

## Solution: Truncated Bonnor-Ebert Sphere

To achieve true Jeans stability, use a **truncated** BE sphere with xi_s < 3.54.

### New Stable Configuration
| Parameter | Value |
|-----------|-------|
| T_cloud | **20 K** |
| n_center | **3000 cm⁻³** |
| xi_s | **3.0** (truncated) |
| mu | 1.27 |

### Derived Parameters (from BE generator)
| Parameter | Value |
|-----------|-------|
| r_0 (scale) | 0.160 pc |
| R_cloud | 0.48 pc |
| M_cloud | 22.5 M_sun |
| n_edge | 1034 cm⁻³ |
| P_ext/k_B | 20,689 K cm⁻³ |

### Stability Verification

**Jeans Stability:**
- λ_J = 1.0 pc (at T=20K, n=3000)
- M_J = 50 M_sun
- **R/λ_J = 0.48** ✓ (< 0.8)
- **M/M_J = 0.45** ✓ (< 0.8)

**Thermal Stability:**
- T_cloud = 20 K
- T_eq (K&I 2000 at n=3000) ≈ 6 K
- Margin = +14 K above equilibrium
- Cloud will cool but **NOT collapse**

**Stability at Equilibrium (T → 6K):**
- λ_J(6K) ≈ 0.55 pc
- R/λ_J ≈ 0.87 (still < 1, stable)
- M_J(6K) ≈ 8 M_sun (cloud approaches this limit)
- Thermal equilibrium prevents further cooling

## CRITICAL UPDATE: BE Sphere + K&I 2000 Cooling Incompatibility

### The Fundamental Problem

K&I 2000 cooling produces VERY cold equilibrium temperatures at high densities. This makes the Jeans mass extremely small, causing ANY realistic BE sphere to become super-Jeans after cooling.

### Stability Table: BE Sphere vs K&I 2000 Cooling

| n_center (cm⁻³) | T_eq (K) | M_J(T_eq) (M☉) | M_BE (M☉) | M/M_J | Status |
|-----------------|----------|----------------|-----------|-------|--------|
| 100 | 80 | 2192 | 1276 | 0.58 | ✓ STABLE |
| 200 | 45 | 654 | 453 | 0.69 | ✓ STABLE |
| 300 | 28 | 262 | 224 | 0.86 | ⚠ RISKY |
| 500 | 17 | 96 | 112 | 1.16 | ✗ UNSTABLE |
| 750 | 13 | 52 | 75 | 1.42 | ✗ UNSTABLE |
| 1000 | 11 | 35 | 58 | 1.63 | ✗ UNSTABLE |
| 1500 | 9 | 21 | 45 | 2.08 | ✗ UNSTABLE |
| 2000 | 8 | 16 | 39 | 2.48 | ✗ UNSTABLE |
| 3000 | 7 | 10 | 31 | 3.04 | ✗ UNSTABLE |
| 5000 | 4.5 | 4.1 | 24 | 5.89 | ✗ UNSTABLE |

**Key insight**: For n > 300 cm⁻³, NO standard BE sphere can survive K&I 2000 cooling!

### Maximum Safe Density

- **n_max ≈ 200 cm⁻³** for stable cooling
- T_eq = 45 K at this density
- M_cloud ~ 450 M☉, R ~ 3.2 pc
- M/M_J = 0.69 after cooling (STABLE)

### If High Density Required (n ~ 3000 cm⁻³)

Three options exist:

1. **Ultra-small cloud**: M < 8 M☉ (use ξ_s ~ 2.0) - impractical
2. **Start near equilibrium**: T_init ~ 8K - no cooling to test
3. **Disable self-gravity**: G = 0 - test cooling in isolation

## Alternative: K&I 2000 Compatible Sphere

See section below for a thermally self-consistent approach.

## Configuration Files

### Phase 1: Relaxation (Truncated BE)
```
simulations/astrophysics/imbh_cloud/config/presets/phase1_stable.json
```
- Creates truncated BE sphere with xi_s=3.0
- Runs damped relaxation to hydrostatic equilibrium
- Output: snapshot_0021.csv

### Phase 2.5: Cooling Test
```
simulations/astrophysics/imbh_cloud/config/presets/phase2.5_stable_cooling.json
```
- Resumes from Phase 1 relaxed state
- Enables K&I 2000 cooling
- Should run stably without timestep collapse

## Physical Interpretation

The truncated BE sphere (xi_s=3.0) represents a cloud that is:

1. **Pressure-confined** by external medium (not self-gravitating to the edge)
2. **Sub-critical** in the BE stability sense (not on verge of collapse)
3. **Sub-Jeans** in both radius and mass (gravitationally stable)

This is physically reasonable for a dense cloud embedded in the warm ISM, where external pressure truncates the cloud before it reaches the critical BE radius.

## Comparison: Critical vs Truncated BE

| Property | Critical (xi_s=6) | Truncated (xi_s=3) |
|----------|-------------------|---------------------|
| Density contrast | 14:1 | 2.9:1 |
| R/λ_J | ~0.95 | ~0.48 |
| Self-gravity | Dominant | Moderate |
| Stability | Marginal | Stable |
| Physical meaning | Isolated cloud | Pressure-confined |

## Recommendations

1. **For stable ISM simulations**: Use xi_s ≤ 3.0
2. **For collapse studies**: Use xi_s ≥ 6.0 (intentionally unstable)
3. **Always verify**: Check R/λ_J < 0.8 AND M/M_J < 0.8 before running
4. **Account for cooling**: If T will decrease, check stability at T_eq

## K&I 2000 Compatible Sphere Analysis

### Question: Can We Use a Barotropic Sphere Instead of BE?

A Bonnor-Ebert sphere assumes **isothermal** gas (constant T). But what if we use a sphere where T = T_eq(n) according to K&I 2000? This would be a **barotropic** hydrostatic sphere.

The KoyamaInutsukaRelaxation module in the codebase already implements this approach.

### The Physics of Barotropic Spheres

For a barotropic EOS where P = P(ρ), the effective sound speed is:

$$c_{\text{eff}}^2 = \frac{dP}{d\rho} = \frac{k_B T}{m_n} \left(1 + \frac{d \ln T}{d \ln n}\right)$$

This differs from the isothermal sound speed c_s² = k_B T / m_n by the factor (1 + d ln T / d ln n).

### K&I 2000 Cold Branch: The Killer Problem

On the K&I 2000 cold branch (high density), temperature **decreases** with increasing density:

| n (cm⁻³) | T_eq (K) | d ln T / d ln n | c_eff / c_s | M_J_eff / M_J_iso |
|----------|----------|-----------------|-------------|-------------------|
| 100 | 80 | -0.31 | 0.83 | 0.57 |
| 200 | 45 | -1.02 | ~0 | ~0 |
| 300 | 28 | -1.12 | ~0 | ~0 |
| 500 | 17 | -0.92 | 0.29 | 0.02 |
| 1000 | 11 | -0.54 | 0.68 | 0.31 |
| 2000 | 8 | -0.37 | 0.79 | 0.50 |
| 3000 | 7 | -0.53 | 0.69 | 0.33 |
| 5000 | 4.5 | -0.67 | 0.58 | 0.19 |

### Key Insight: K&I 2000 Sphere Makes the Problem WORSE

Since d ln T / d ln n < 0 on the cold branch:
- c_eff² < c_s² (effective sound speed is reduced)
- c_eff ~ 0.7 × c_s typically
- M_J_eff ~ c_eff³ ~ 0.35 × M_J_iso

**The effective Jeans mass is SMALLER than isothermal!**

This means a K&I 2000 barotropic sphere is **MORE prone to collapse** than an equivalent isothermal BE sphere, not less.

### Why This Makes Physical Sense

On the cold branch:
1. As the core compresses, n increases
2. T_eq(n) **decreases** (cooling wins over compression heating)
3. This provides **less** thermal pressure resistance than isothermal
4. The gas is effectively **softer** than isothermal
5. Gravitational collapse is **easier**, not harder

### The Fundamental Conflict

K&I 2000 equilibrium was designed for **shock-compressed layers** (2D), not **self-gravitating spheres** (3D):

| K&I 2000 Context | Spherical Cloud Context |
|------------------|-------------------------|
| 2D compression | 3D compression |
| Thermal equilibrium | Hydrostatic equilibrium |
| No self-gravity in model | Self-gravity essential |
| n ~ 1000 cm⁻³ in dense phase | n ~ 1000-3000 cm⁻³ target |
| Works beautifully | Jeans unstable |

### Conclusion

**A K&I 2000 compatible sphere does NOT solve the Jeans instability problem.**

The only viable options for high-density (n > 300 cm⁻³) simulations are:

1. **Disable self-gravity (G = 0)**: Test cooling physics in isolation
2. **Very small mass (M < 8 M☉)**: Impractical resolution requirements
3. **Accept lower density (n ~ 200 cm⁻³)**: T_eq = 45 K, M ~ 450 M☉

For realistic ISM physics with self-gravity, the maximum viable density is **n ~ 200 cm⁻³**.

## Parameter Survey Collapse Analysis (2026-01-13)

### Observed Behavior

During the IMBH-cloud parameter survey with r_p = 0.4 pc and K&I 2000 cooling enabled:

| Metric | Adiabatic | Cooling |
|--------|-----------|---------|
| Snapshots completed | 68 | 10 (stuck) |
| Max density | 2.3×10⁴ | **3.5×10⁷** |
| Max temperature | 691,821 K | 63.5 K |
| Timestep at end | ~10⁻³ | **10⁻⁸** (collapsed) |

### First Principles Derivation

#### 1. Jeans Instability from Fluid Equations

Starting from the linearized fluid equations with self-gravity:

**Continuity:** ∂ρ/∂t + ∇·(ρv) = 0

**Momentum:** ∂v/∂t = -∇P/ρ - ∇Φ

**Poisson:** ∇²Φ = 4πGρ

Linearizing with ρ = ρ₀ + ρ₁, assuming plane wave solutions ∝ exp[i(k·r - ωt)]:

$$\omega^2 = c_s^2 k^2 - 4\pi G \rho_0$$

Instability occurs when ω² < 0, giving the **Jeans wavenumber**:

$$k_J = \sqrt{\frac{4\pi G \rho_0}{c_s^2}}$$

And the **Jeans mass**:

$$M_J = \frac{\pi^{5/2}}{6} \frac{c_s^3}{G^{3/2} \rho_0^{1/2}}$$

#### 2. Effect of Equation of State on Collapse

For compression (∇·v < 0), the energy equation gives:

$$\frac{du}{dt} = -\frac{P}{\rho} \nabla \cdot v - \Lambda(T, \rho)$$

**Case A: Adiabatic (Λ = 0)**

PdV work heats the gas: T ∝ ρ^(γ-1)

Since c_s² ∝ T:
$$M_J \propto \frac{c_s^3}{\rho^{1/2}} \propto \rho^{\frac{3(\gamma-1)}{2} - \frac{1}{2}} = \rho^{\frac{3\gamma - 4}{2}}$$

For γ = 5/3: **M_J ∝ ρ^(1/6)** → Jeans mass INCREASES with compression → STABLE

**Case B: Isothermal/Efficient Cooling**

Temperature remains constant, cooling removes PdV heat:
$$M_J \propto \rho^{-1/2}$$

**Jeans mass DECREASES with compression → RUNAWAY COLLAPSE**

#### 3. Quantitative Analysis of Observed Collapse

Initial condition (from resume snapshot):
- Core density: ρ₀ = 1318 code units
- Temperature: T ≈ 20 K
- Sound speed: c_s ≈ 0.27 km/s
- Core mass: M_core ≈ 4.0 M☉

**CRITICAL: Unit Conversion**

The simulation uses density_to_n_H = 31.9 cm⁻³ per code unit:
- Code density (max): 1318 → **n_H = 42,044 cm⁻³**
- Code density (median): 5.6 → n_H = 178 cm⁻³

**The core density is 210× higher than the stability limit of 200 cm⁻³!**

This explains why the config setting "n_center = 500 cm⁻³" doesn't match reality:
1. The BE sphere was generated with n_center = 500 cm⁻³
2. But the snapshot was from phase2.5 which already underwent cooling-driven compression
3. The core has compressed to n ~ 42,000 cm⁻³ before the parameter survey even started

Initial Jeans mass:
$$M_J = \frac{c_s^3}{G^{3/2} \sqrt{\rho_0}} = \frac{(0.27)^3}{(0.0043)^{3/2} \sqrt{1318}} \approx 1.9 \, M_\odot$$

**Stability ratio: M_core/M_J = 4.0/1.9 ≈ 2.1 > 1 → UNSTABLE**

#### 4. Free-Fall Time Verification

$$t_{ff} = \sqrt{\frac{3\pi}{32 G \rho}} \approx 0.21 \text{ code time}$$

Observed collapse time: t ≈ 0.18 code time — **matches free-fall prediction!**

#### 5. Timestep Collapse Mechanism

CFL condition: Δt < C × h/c_s

Smoothing length scales as: h ∝ ρ^(-1/3)

Therefore: Δt_max ∝ ρ^(-1/3)

When density increases 10⁴× (from 10³ to 10⁷):
- Timestep decreases by 10^(4/3) ≈ 20×
- But collapse accelerates (shorter t_ff)
- Creates positive feedback → **timestep death spiral**

### Summary Table

| Quantity | Formula | Adiabatic | Cooling |
|----------|---------|-----------|---------|
| T response to compression | - | T ↑↑ (691,000 K) | T = const (63 K) |
| M_J scaling | M_J ∝ T^(3/2) ρ^(-1/2) | ∝ ρ^(+1/6) | ∝ ρ^(-1/2) |
| Stability | - | Self-limiting | Runaway |
| Outcome | - | Runs to completion | Timestep collapse |

### Recommendations for Parameter Survey

1. **Use adiabatic runs** for studying tidal disruption dynamics
2. **If cooling required**, either:
   - Start with lower density cloud (n < 200 cm⁻³)
   - Add temperature floor to prevent runaway cooling
   - Implement sink particles for collapsing regions
3. **Monitor Jeans ratio** M/M_J throughout simulation

## References

- Bonnor, W.B. (1956). "Boyle's Law and gravitational instability"
- Ebert, R. (1955). "Über die Verdichtung von H I-Gebieten"
- Koyama, H. & Inutsuka, S. (2000). "Molecular Cloud Formation in Shock-compressed Layers"
- Truelove et al. (1997). "The Jeans Condition" - Resolution requirements for collapse

---
*Generated: 2026-01-12*
*Updated: 2026-01-13 - Added K&I 2000 compatible sphere analysis*
*Updated: 2026-01-13 - Added parameter survey collapse analysis with first principles derivation*
*Analysis scripts:*
- *simulations/astrophysics/ism_cooling/scripts/design_stable_cloud.py*
- *scripts/find_optimal_density.py*
- *scripts/analyze_ki2000_sphere.py*
