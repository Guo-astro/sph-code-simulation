# IMBH-Molecular Cloud Interaction: Research Setup

## Scientific Objective

Simulate intermediate-mass black hole (IMBH) tidal disruption of molecular clouds with the following parameters:

### Physical Parameters

- **IMBH Mass**: M_BH = 10⁵ M_☉ (100,000 solar masses)
- **Impact Parameters**: b = 3-6 parsecs (parameter scan study)
- **Initial Relative Velocity**: v₀ = 0-20 km/s (velocity effect study)
  - Baseline: v₀ = 10 km/s in x-direction
  - Additional scans: v₀ = 0, 5, 10, 15, 20 km/s
- **Cloud Structure**: 3D Lane-Emden sphere with polytropic index n = 5/3 (γ = 5/3)
- **Thermal Physics**: Koyama & Inutsuka (2000) thermal equilibrium curve
  - Assumption: Thermal timescale << dynamical timescale
  - Cloud remains in thermal equilibrium throughout simulation

### Research Questions

1. **Impact Parameter Study**: How does tidal disruption efficiency vary with impact parameter b?
   - b = 3 pc: Strong disruption (b < r_tidal)
   - b = 4 pc: Moderate disruption
   - b = 5 pc: Weak disruption
   - b = 6 pc: Minimal disruption (reference case)

2. **Velocity Effect Study**: How does initial relative velocity affect:
   - Encounter timescale and disruption morphology
   - Mass accretion rate onto BH
   - Tidal tail formation
   - Shock heating at pericenter passage

## Resolution Analysis

### 1. Jeans Mass Criterion

The Jeans mass for a self-gravitating cloud is:

```
M_J = (π^(5/2) / 6) * (c_s^3 / (G^(3/2) * ρ^(1/2)))
```

For thermal equilibrium in CNM (Cold Neutral Medium):
- n_H ~ 100-1000 cm⁻³ → T_eq ~ 10-50 K (from K&I thermal curve)
- Sound speed: c_s ~ 0.2-0.3 km/s

**Numerical example** (n_H = 100 cm⁻³, T = 50 K):
- ρ = n_H × μ m_H / pc³ = 100 × (1.67×10⁻²⁴ g) / (2.938×10⁵⁷ cm³) ≈ 2.47 M_☉/pc³
- c_s = √(k_B T / μ m_H) = √(1.38×10⁻¹⁶ × 50 / 1.67×10⁻²⁴) ≈ 0.29 km/s
- **M_J ≈ 3.5 M_☉** (typical for molecular clouds)

**Requirement**: To resolve fragmentation and thermal pressure support:
- **N_J ≥ 50-100 particles per Jeans mass**
- For cloud mass M_cloud ~ 10³-10⁴ M_☉:
  - **Minimum N ≥ 5×10⁴ - 10⁶ particles**

### 2. Tidal Disruption Timescale

The tidal (pancake) disruption occurs when tidal force exceeds self-gravity:

```
t_tidal ~ √(R_cloud³ / (G M_BH))
```

**Numerical example** (M_BH = 10⁵ M_☉, R_cloud = 5 pc):
```
t_tidal = √(5³ pc³ / (4.302×10⁻³ pc·(km/s)²·M_☉⁻¹ × 10⁵ M_☉))
        = √(125 / 430.2) [Myr]
        ≈ 0.54 Myr = 5.4×10⁵ years
```

Compare to cloud crossing time:
```
t_cross ~ R_cloud / v_rel
```

**Numerical examples**:
- v_rel = 10 km/s: t_cross = 5 pc / 10 km/s = 0.49 Myr = **4.9×10⁵ years**
- v_rel = 20 km/s: t_cross = 5 pc / 20 km/s = 0.24 Myr = **2.4×10⁵ years**
- v_rel = 50 km/s: t_cross = 5 pc / 50 km/s = 0.10 Myr = **1.0×10⁵ years**

**Timescale ordering**: t_thermal << t_tidal ~ t_cross

### 2a. Velocity Effect on Encounter Timescale

**Pericenter passage time** (when tidal force is strongest):
```
t_peri ~ b / v_rel
```

**Numerical examples** (various b and v_rel combinations):

| b [pc] | v_rel [km/s] | t_peri [Myr] | t_peri [years] | Physical interpretation |
|--------|--------------|--------------|----------------|------------------------|
| 3 | 0 | ∞ | ∞ | Static infall - no pericenter, pure radial |
| 3 | 5 | 0.59 | 5.9×10⁵ | Slow, extended tidal compression |
| 3 | 10 | 0.29 | 2.9×10⁵ | Baseline - moderate encounter |
| 3 | 15 | 0.20 | 2.0×10⁵ | Fast - impulsive tidal shock |
| 3 | 20 | 0.15 | 1.5×10⁵ | Very fast - minimal disruption |
| 4 | 10 | 0.39 | 3.9×10⁵ | Baseline at moderate distance |
| 5 | 10 | 0.49 | 4.9×10⁵ | Weak tidal interaction |
| 6 | 10 | 0.59 | 5.9×10⁵ | Minimal tidal effect (reference) |

**Velocity regimes**:
- **v_rel = 0 km/s** (static case): t_peri → ∞, pure radial infall, maximum disruption
- **v_rel = 5 km/s**: t_peri ~ 0.6 Myr, slow encounter, extended tidal compression
- **v_rel = 10 km/s** (BASELINE): t_peri ~ 0.3 Myr, intermediate regime
- **v_rel = 15 km/s**: t_peri ~ 0.2 Myr, fast encounter, impulsive tidal shock
- **v_rel = 20 km/s**: t_peri ~ 0.15 Myr, very fast, minimal disruption

**Key physics**: Higher v_rel → shorter t_peri → less time for tidal stripping → weaker disruption

### 3. Spatial Resolution Requirements

**Minimum spatial resolution needed**:

a) **Tidal radius** (where BH gravity ~ cloud self-gravity):
   ```
   r_t ~ R_cloud * (M_BH / M_cloud)^(1/3)
   ```
   **Numerical example** (M_cloud = 10⁴ M_☉, M_BH = 10⁵ M_☉, R_cloud = 5 pc):
   ```
   r_t = 5 pc × (10⁵ / 10⁴)^(1/3)
       = 5 pc × 10^(1/3)
       = 5 pc × 2.15
       ≈ 10.8 pc
   ```
   - **Need to resolve h ~ 0.1-1 pc** to capture tidal effects

b) **Hill radius** (tidal truncation):
   ```
   r_H ~ b * (M_cloud / (3 M_BH))^(1/3)
   ```
   **Numerical examples** (M_cloud = 10⁴ M_☉, M_BH = 10⁵ M_☉):
   
   | b [pc] | r_H [pc] | Interpretation |
   |--------|----------|----------------|
   | 3 | 0.99 | Strong tidal truncation |
   | 4 | 1.32 | Moderate truncation |
   | 5 | 1.65 | Weak truncation |
   | 6 | 1.98 | Minimal truncation |
   
   Calculation: r_H = b × (10⁴ / (3×10⁵))^(1/3) = b × 0.3264 pc
   - **Need N_neighbor ~ 50-100 within r_H** for proper resolution

c) **Smoothing length constraint**:
   - h ~ R_cloud / √N for uniform distribution
   
   **Numerical examples** (R_cloud = 5 pc):
   
   | N particles | h [pc] | h/R_cloud | Resolution quality |
   |-------------|--------|-----------|-------------------|
   | 5×10⁴ | 0.022 | 0.0045 | Quick test ✓ |
   | 10⁵ | 0.016 | 0.0032 | Good ✓ |
   | 2×10⁵ | 0.011 | 0.0022 | Standard ✓✓ |
   | 10⁶ | 0.005 | 0.0010 | High resolution ✓✓✓ |

### 4. Recommended Particle Numbers

**Quick Reference Table**:

| Cloud Mass | N_particles | h/R_cloud | N_J | Purpose |
|------------|-------------|-----------|-----|---------|
| 10³ M_☉   | 5×10⁴       | ~0.014    | 50  | Quick test |
| 10⁴ M_☉   | 2×10⁵       | ~0.007    | 100 | Standard resolution |
| 10⁴ M_☉   | 10⁶         | ~0.003    | 500 | High resolution |

**Recommended**: Start with **N = 2×10⁵** for production runs, **N = 5×10⁴** for parameter studies.

---

### 4a. Minimal Acceptable Particle Numbers (For Reviewers)

The following table provides **justifiable minimal particle counts** for different physical scenarios, based on multiple resolution criteria that reviewers will scrutinize:

#### **Resolution Criteria for Minimal N**:
1. **Jeans criterion**: N_J ≥ 50 particles per Jeans mass (fragmenttion resolution)
2. **Neighbor criterion**: N_neighbor ≥ 50 within smoothing kernel (SPH accuracy)
3. **Hill sphere**: ≥50 particles within r_H at pericenter (tidal truncation)
4. **Spatial resolution**: h ≤ 0.1 r_H (resolve tidal features)

---

#### **A. Cloud Mass Variation Study** (Fixed: b = 4 pc, v = 10 km/s, R ∝ M^(1/3))

| M_cloud | R_cloud | <n_H> | M_J | r_H | N_min | N_recommended | Justification |
|---------|---------|-------|-----|-----|-------|---------------|--------------|
| 10² M_☉ | 2.3 pc | 100 cm⁻³ | 3.5 M_☉ | 0.68 pc | 1.4×10³ | 1×10⁴ | M/M_J = 29 → 29×50 = 1450; h(10⁴)=0.023pc < r_H ✓ |
| 5×10² M_☉ | 3.7 pc | 100 cm⁻³ | 3.5 M_☉ | 0.91 pc | 7.1×10³ | 3×10⁴ | M/M_J = 143 → need 7150; increased for h < 0.1×r_H |
| 10³ M_☉ | 4.6 pc | 100 cm⁻³ | 3.5 M_☉ | 1.08 pc | 1.4×10⁴ | 5×10⁴ | **Quick test baseline**; h(5×10⁴)=0.021pc ✓ |
| 5×10³ M_☉ | 7.4 pc | 100 cm⁻³ | 3.5 M_☉ | 1.46 pc | 7.1×10⁴ | 1.5×10⁵ | M/M_J = 1429 → larger cloud needs more resolution |
| 10⁴ M_☉ | 9.3 pc | 100 cm⁻³ | 3.5 M_☉ | 1.73 pc | 1.4×10⁵ | 2×10⁵ | **Production baseline**; h(2×10⁵)=0.021pc < 0.1×r_H ✓ |
| 5×10⁴ M_☉ | 15 pc | 100 cm⁻³ | 3.5 M_☉ | 2.33 pc | 7.1×10⁵ | 1×10⁶ | High-mass cloud, publication quality |

**Scaling**: N_min ≈ (M_cloud / M_J) × 50, rounded up to nearest convenient number

**Calculation example** (M = 10⁴ M_☉):
```
M/M_J = 10⁴ / 3.5 ≈ 2857 Jeans masses
N_min = 2857 × 50 = 142,850 ≈ 1.4×10⁵
Recommended = 2×10⁵ (40% safety margin)
h = R/√N = 9.3pc / √(2×10⁵) ≈ 0.021 pc < 0.1 × r_H(1.73pc) ✓
```

---

#### **B. Impact Parameter Variation** (Fixed: M = 10⁴ M_☉, R = 5 pc, v = 10 km/s)

| b [pc] | r_H [pc] | h_max | N_min | N_recommended | Reviewer argument |
|--------|----------|-------|-------|---------------|-------------------|
| 2.0 | 0.66 | 0.066 | 5.8×10⁵ | 8×10⁵ | **Very close**: h(5.8×10⁵)=0.0066pc; 50 ptcls in r_H requires extreme resolution |
| 3.0 | 0.99 | 0.099 | 2.6×10⁵ | 4×10⁵ | **Close encounter**: h ≤ 0.1pc to resolve Hill sphere; strong tidal stripping |
| 4.0 | 1.32 | 0.132 | 1.4×10⁵ | 2×10⁵ | **Baseline**: h(2×10⁵)=0.011pc < 0.1×r_H ✓; 50 neighbors guaranteed |
| 5.0 | 1.65 | 0.165 | 9.2×10⁴ | 1.5×10⁵ | Weaker tidal effect; slightly relaxed but still resolve r_H |
| 6.0 | 1.98 | 0.198 | 6.4×10⁴ | 1×10⁵ | **Distant**: minimal disruption; lower N acceptable for reference |
| 8.0 | 2.64 | 0.264 | 3.6×10⁴ | 5×10⁴ | Control case: very weak tidal effect |

**Calculation** (h_max = 0.1 × r_H):
```
For b = 4 pc: r_H = 1.32 pc → h_max = 0.132 pc
N_min = (R_cloud / h_max)² = (5 / 0.132)² ≈ 1.4×10⁵
Plus Jeans: M/M_J = 2857, need 142,850
Maximum: N_min = max(1.4×10⁵, 1.4×10⁵) ≈ 1.4×10⁵
Recommended: 2×10⁵ (add safety margin)
```

**Reviewer note**: Close encounters (b < 4 pc) require **higher resolution** because:
1. Smaller r_H → need smaller h → need more particles
2. Stronger tidal gradients → need better force resolution
3. Higher compression → more Jeans masses form → need more N_J

---

#### **C. Velocity Variation** (Fixed: M = 10⁴ M_☉, R = 5 pc, b = 4 pc)

| v [km/s] | t_peri [Myr] | Δt_output | N_snapshots | N_particles | Justification |
|----------|--------------|-----------|-------------|-------------|---------------|
| 0 (static) | ∞ | 0.02 | 100 | 2×10⁵ | **Quasi-static infall**: slow evolution, standard resolution adequate |
| 5 | 0.78 | 0.01 | 200 | 3×10⁵ | Slow encounter: t_peri ~ 1.5×t_tidal → need finer time sampling + spatial |
| 10 | 0.39 | 0.005 | 400 | 2×10⁵ | **Baseline**: moderate timescale, standard resolution |
| 15 | 0.26 | 0.003 | 600 | 2.5×10⁵ | Fast encounter: impulsive shock → need good Riemann solver + slightly more N |
| 20 | 0.20 | 0.002 | 1000 | 3×10⁵ | **Very fast shock**: ΔT_shock ~ 160K → resolve shock heating → higher N |
| 30 | 0.13 | 0.001 | 1300 | 4×10⁵ | Extreme case: shock-dominated, publication-level resolution |

**Velocity-dependent N**: Higher v requires:
1. **More particles** to resolve shock fronts (Δx_shock ~ few h)
2. **Finer timesteps** (CFL condition: Δt ∝ h / v_shock)
3. **More snapshots** to capture rapid pericenter passage

**Calculation** (v = 20 km/s shock):
```
Shock width: λ_shock ~ 5-10 h (numerical diffusion)
v_shock ~ v_rel = 20 km/s
To resolve shock: need h < λ_shock / 5
Post-shock compression: n → 4n (strong shock) → 4× more Jeans masses
N_min = 1.4×10⁵ × √2 ≈ 2×10⁵ (increased by ~√2 for shock)
Recommended: 3×10⁵ (50% safety for shock physics)
```

---

#### **D. Cloud Radius Variation** (Fixed: M = 10⁴ M_☉, b = 4 pc, v = 10 km/s)

| R_cloud | <ρ> | <n_H> | M_J | r_H/R | N_min | N_recommended | Reviewer concern |
|---------|-----|-------|-----|-------|-------|---------------|------------------|
| 3 pc | 88 M_☉/pc³ | 3564 cm⁻³ | 1.3 M_☉ | 0.44 | 3.8×10⁵ | 5×10⁵ | **Dense compact**: M/M_J=7692 → many Jeans masses |
| 4 pc | 37 M_☉/pc³ | 1500 cm⁻³ | 1.9 M_☉ | 0.33 | 2.6×10⁵ | 3.5×10⁵ | Moderate density, resolve r_H well |
| 5 pc | 19 M_☉/pc³ | 775 cm⁻³ | 2.7 M_☉ | 0.26 | 1.9×10⁵ | 2.5×10⁵ | **Baseline** |
| 7 pc | 6.9 M_☉/pc³ | 280 cm⁻³ | 4.5 M_☉ | 0.19 | 1.1×10⁵ | 1.5×10⁵ | Lower density, fewer Jeans masses |
| 10 pc | 2.4 M_☉/pc³ | 97 cm⁻³ | 7.6 M_☉ | 0.13 | 6.6×10⁴ | 1×10⁵ | **Diffuse**: M/M_J=1316, relaxed requirement |

**Key insight**: Smaller R_cloud at fixed M:
- Higher <ρ> → lower M_J → **more Jeans masses** → need more particles
- r_H/R smaller → tidal effect covers larger fraction → need better resolution

**Calculation** (R = 3 pc, M = 10⁴ M_☉):
```
<ρ> = M / (4π/3 R³) = 10⁴ / 113 ≈ 88 M_☉/pc³ → n_H ≈ 3564 cm⁻³
M_J(3564 cm⁻³) ≈ 1.3 M_☉ (higher density → smaller Jeans mass)
M/M_J = 10⁴ / 1.3 ≈ 7692 Jeans masses
N_min = 7692 × 50 ≈ 3.8×10⁵
```

---

#### **E. Combined Worst-Case Scenarios**

| Scenario | M | R | b | v | N_min | N_recommended | Why this is challenging |
|----------|---|---|---|---|-------|---------------|------------------------|
| **Dense close fast** | 10⁴ M_☉ | 3 pc | 3 pc | 20 km/s | 8×10⁵ | 1×10⁶ | Small M_J + small r_H + shock → maximum resolution |
| **Massive close slow** | 5×10⁴ M_☉ | 15 pc | 3 pc | 5 km/s | 9×10⁵ | 1.2×10⁶ | Many M_J + extended tidal tails |
| **Standard publication** | 10⁴ M_☉ | 5 pc | 4 pc | 10 km/s | 1.4×10⁵ | 2×10⁵ | **Baseline for papers** |
| **Quick exploratory** | 10³ M_☉ | 4.6 pc | 5 pc | 10 km/s | 1.4×10⁴ | 5×10⁴ | Parameter scan, proof of concept |
| **Control (weak)** | 10⁴ M_☉ | 7 pc | 6 pc | 10 km/s | 8×10⁴ | 1×10⁵ | Reference case, minimal disruption |

**Worst-case calculation** (Dense close fast):
```
M = 10⁴ M_☉, R = 3 pc → <n_H> = 3564 cm⁻³ → M_J = 1.3 M_☉
N_Jeans = (10⁴ / 1.3) × 50 = 3.8×10⁵

b = 3 pc → r_H = 0.99 pc
h_max = 0.1 × r_H = 0.099 pc
N_spatial = (3 pc / 0.099 pc)² = 9.2×10⁵

v = 20 km/s → shock resolution: N_shock ≈ 1.5 × N_baseline = 2.1×10⁵

N_min = max(3.8×10⁵, 9.2×10⁵, 2.1×10⁵) ≈ 9.2×10⁵
Recommended: 1×10⁶ (round up, add 10% safety)
```

---

#### **F. Summary: Minimal N vs Recommended N**

| Use Case | N_min (bare minimum) | N_recommended | Safety margin |
|----------|---------------------|---------------|---------------|
| Quick parameter scan | 5×10⁴ | 5×10⁴ | None (accept lower quality) |
| Exploratory run | 1×10⁵ | 1.5×10⁵ | +50% |
| **Standard publication** | 1.4×10⁵ | **2×10⁵** | +43% |
| High-quality paper | 2×10⁵ | 4×10⁵ | +100% |
| Flagship simulation | 5×10⁵ | 1×10⁶ | +100% |

**Reviewer checklist** (what they will scrutinize):
1. ✓ **Jeans resolution**: N_J ≥ 50 particles per M_J
2. ✓ **Neighbor count**: N_neighbor ≥ 50 within smoothing kernel
3. ✓ **Tidal feature**: h ≤ 0.1 r_H (resolve Hill sphere)
4. ✓ **Convergence test**: Compare N and 2N results (< 10% difference)
5. ✓ **Energy conservation**: ΔE/E₀ < 1% over simulation
6. ✓ **Shock resolution**: At least 5-10 particles across shock front (if v > 15 km/s)

**Recommended approach for papers**:
- Main results: N = 2×10⁵ (baseline) + N = 4×10⁵ (convergence test)
- Parameter scan: N = 1×10⁵ (acceptable) → upgrade to 2×10⁵ for final published figures
- Supplementary: N = 5×10⁴ (proof of concept only, clearly labeled)

### 5. Thermal Equilibrium Constraint

From Koyama & Inutsuka (2000) cooling timescale:

```
t_cool ~ u / |du/dt| ~ u / |(u_eq - u) / τ_relax|
```

**Numerical example** (τ_relax = 0.05 Myr, t_tidal = 0.54 Myr):
```
τ_relax / t_tidal = 0.05 / 0.54 ≈ 0.09 = 9%
```
- **Thermal equilibrium maintained**: T → T_eq(n) within ~10% of t_tidal
- **Cooling is ~10× faster than tidal disruption** → quasi-equilibrium valid

**Specific internal energy** (T = 50 K, γ = 5/3):
```
u = k_B T / ((γ-1) μ m_H)
  = (1.38×10⁻¹⁶ erg/K × 50 K) / (0.667 × 1.67×10⁻²⁴ g)
  ≈ 6.2×10⁹ erg/g
  ≈ 6.2×10¹² cm²/s² (code units: ~0.062 in [pc²/Myr²])
```

**Implementation**: Use `KoyamaInutsukaCooling` with relaxation time:
```cpp
real tau_relax = 0.05;  // Myr (5% of tidal timescale ~0.5 Myr)
real du_dt = cooling.cooling_rate(n_H, T_current, tau_relax);
// For T >> T_eq: du/dt ≈ -(u - u_eq) / tau_relax ≈ -u / 0.05 Myr
```

## Physics Modules Required

### 1. External Force: Point-Mass Black Hole

**Header**: `include/external_forces/point_mass_bh.hpp`

```cpp
namespace sph {
namespace external_forces {

class PointMassBlackHole {
    vec_t m_position;      // BH position [pc]
    real m_mass;           // BH mass [M_☉]
    real m_softening;      // Gravitational softening [pc]
    
public:
    vec_t acceleration(const vec_t& r_particle) const;
    real potential(const vec_t& r_particle) const;
};

}
}
```

**Force**:
```
F_BH = -G M_BH (r - r_BH) / (|r - r_BH|² + ε²)^(3/2)
```

where ε = softening length ~ 0.01 R_cloud (prevents singularity).

**Numerical example** (M_BH = 10⁵ M_☉, r = 3 pc from BH, ε = 0.05 pc):
```
|F_BH| = G M_BH / r²
       = (4.302×10⁻³ pc·(km/s)²·M_☉⁻¹) × 10⁵ M_☉ / (3 pc)²
       = 430.2 / 9 [km²/s²/pc]
       ≈ 47.8 km²/s²/pc
       ≈ 47.8 pc/Myr² (code units: ~47.8)

Acceleration: a_BH = F_BH / m_particle
```

**Gravitational potential energy** (cloud at b = 3 pc):
```
E_pot = -G M_BH M_cloud / b
      = -(4.302×10⁻³) × 10⁵ × 10⁴ / 3 [M_☉·(km/s)²]
      ≈ -1.43×10⁸ M_☉·(km/s)²
      ≈ -2.85×10⁵⁴ erg (enormous binding energy!)
```

### 2. Thermal Physics: K&I Cooling

**Existing module**: `include/thermal/koyama_inutsuka_cooling.hpp`

- Already implemented ✓
- Provides T_eq(n_H), P_eq(n_H)
- Cooling rate: du/dt = (u_eq - u) / τ_relax

### 3. Initial Conditions: Lane-Emden + Thermal Equilibrium

**Based on**: `src/sample/lane_emden.cpp`

Modifications needed:
1. Set initial T = T_eq(n) from K&I curve (not polytro

pic P = K ρ^γ)
2. Add IMBH at offset position
3. Optionally add relative velocity

### 4. SPH Method: GDISPH (Recommended)

**Why GDISPH**:
- Best for density gradients (cloud disruption creates steep gradients)
- Hybrid Godunov + pressure-energy formulation
- Built-in Riemann solver handles shocks from tidal compression
- Already tested in `sample/sedov/`, `sample/khi/`

**Configuration**:
```json
{
  "numerical": {
    "sph_type": "gdisph",
    "kernel": "wendland",
    "neighbor_number": 50,
    "use_gravity": true
  },
  "artificial_viscosity": {
    "alpha": 1.0,
    "use_balsara_switch": true,
    "use_time_dependent_av": true
  }
}
```

## Code Units and Conversion

### Recommended Unit System

**Length unit**: 1 code unit = 1 parsec (pc)
**Mass unit**: 1 code unit = 1 M_☉
**Time unit**: Derived from G

```
G = 4.302 × 10⁻³ pc (km/s)² M_☉⁻¹
t_unit = √(L³ / (G M)) = √(pc³ / (G M_☉))
       ≈ 0.978 Myr
```

**Velocity unit**: v_unit = L/t_unit ≈ 1.02 km/s

### Physical Parameters in Code Units

| Quantity | Physical Value | Code Units |
|----------|----------------|------------|
| M_BH | 10⁵ M_☉ | 10⁵ |
| M_cloud | 10⁴ M_☉ | 10⁴ |
| R_cloud | 5 pc | 5.0 |
| b (impact) | 3-6 pc | 3.0-6.0 |
| v_rel (baseline) | 10 km/s | ~9.8 |
| v_rel (slow) | 5 km/s | ~4.9 |
| v_rel (fast) | 15 km/s | ~14.7 |
| v_rel (very fast) | 20 km/s | ~19.6 |
| T_CNM | 50 K | — |
| n_H | 100 cm⁻³ | — |
| t_end | 2 Myr | ~2.0 |

**Velocity conversion**: v[code units] = v[km/s] / 1.02

### Initial Conditions Setup

**Cloud placement**: Center at origin (0, 0, 0)
**BH placement**: Position (b, 0, 0) with velocity (-v_x, 0, 0)
  - Moves toward cloud along x-axis
  - Pericenter passage at x ≈ 0 (closest approach)
**Simulation time**: 2 Myr (sufficient for multiple crossing times)

### Density Conversion

From code density ρ [M_☉/pc³] to number density n_H [cm⁻³]:

```
n_H = ρ * (M_☉/pc³) / (μ m_H)
    = ρ * 40.5  [for μ = 1, neutral H]
```

where:
- μ = 1 (neutral hydrogen)
- m_H = 1.67 × 10⁻²⁴ g
- 1 pc³ = 2.938 × 10⁵⁷ cm³

**Derivation**:
```
n_H = ρ [M_☉/pc³] × (1.989×10³³ g/M_☉) / (2.938×10⁵⁷ cm³/pc³) / (1.67×10⁻²⁴ g)
    = ρ × (1.989×10³³ / 2.938×10⁵⁷ / 1.67×10⁻²⁴)
    = ρ × 40.5 cm⁻³
```

**Numerical examples**:

| ρ [M_☉/pc³] | n_H [cm⁻³] | Physical regime |
|-------------|------------|-----------------|
| 0.1 | 4.05 | Warm neutral medium (WNM) |
| 1.0 | 40.5 | Cold neutral medium (CNM) - low density |
| 2.47 | 100 | **CNM - typical molecular cloud** |
| 10 | 405 | Dense molecular cloud |
| 24.7 | 1000 | Very dense - star forming |

**Cloud mass to density** (M_cloud = 10⁴ M_☉, R_cloud = 5 pc):
```
<ρ> = M_cloud / (4π/3 R_cloud³)
    = 10⁴ M_☉ / (4π/3 × 125 pc³)
    ≈ 19.1 M_☉/pc³
    ≈ 775 cm⁻³ (central density higher!)
```

## Implementation Checklist

### Phase 1: External Force Module
- [ ] `include/external_forces/point_mass_bh.hpp`
- [ ] `src/external_forces/point_mass_bh.cpp`
- [ ] Add to CMakeLists.txt
- [ ] Unit tests

### Phase 2: Initial Conditions
- [ ] `src/sample/imbh_cloud.cpp`
  - [ ] Lane-Emden sphere with n=5/3
  - [ ] K&I thermal equilibrium T_eq(n)
  - [ ] IMBH placement at distance b
  - [ ] Optional relative velocity
  
### Phase 3: Configuration Presets

**Impact Parameter Study**:
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b3pc_v10_gdisph.json`
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b4pc_v10_gdisph.json`
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b5pc_v10_gdisph.json`
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b6pc_v10_gdisph.json`

**Velocity Study** (fixed b = 4 pc):
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b4pc_v0_gdisph.json` (static infall)
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b4pc_v5_gdisph.json`
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b4pc_v10_gdisph.json` (baseline)
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b4pc_v15_gdisph.json`
- [ ] `sample/imbh_cloud/config/presets/imbh_cloud_b4pc_v20_gdisph.json`

### Phase 4: Visualization & Analysis
- [ ] `sample/imbh_cloud/scripts/visualize_disruption.py`
  - [ ] Density evolution
  - [ ] Tidal deformation (axis ratios)
  - [ ] Mass loss rate
  - [ ] Temperature vs density phase diagram
  
- [ ] `sample/imbh_cloud/scripts/compare_impact_parameters.py`
  - [ ] Multi-panel comparison for b = 3, 4, 5, 6 pc
  - [ ] Bound vs unbound mass evolution
  - [ ] Tidal tail morphology differences
  
- [ ] `sample/imbh_cloud/scripts/compare_velocities.py`
  - [ ] Multi-panel comparison for v = 0, 5, 10, 15, 20 km/s
  - [ ] Encounter timescale analysis
  - [ ] Peak compression vs velocity
  
### Phase 5: Makefile
- [ ] `sample/imbh_cloud/Makefile.imbh_cloud`
  - [ ] Single run targets (specific b, v combinations)
  - [ ] Impact parameter scan: `imbh_scan_impact` (4 runs)
  - [ ] Velocity scan: `imbh_scan_velocity` (5 runs)
  - [ ] Full matrix scan: `imbh_scan_full` (20 runs)
  - [ ] Comparison visualization targets

## Expected Physics

### Tidal Disruption Sequence

1. **Initial approach** (t = 0 - 0.3 t_tidal):
   - Cloud feels increasing tidal gradient
   - Elongation along radial direction
   - Compression perpendicular to orbital plane (pancake)
   - **Velocity dependence**: Higher v → faster approach → less pre-compression

2. **Maximum compression** (t ~ t_peri = b/v_rel):
   - Peak density reached at pericenter passage
   - Shocks from rapid compression (stronger for high v_rel)
   - Potential triggered star formation (future work)
   - **Impact parameter effect**: Smaller b → stronger tidal force → more compression
   - **Velocity effect**: Higher v → impulsive shock → heating vs cooling competition

3. **Disruption** (t > t_tidal):
   - Leading/trailing tidal tails form
   - Bound vs unbound material separation
   - Self-gravity determines remnant mass
   - **b-dependence**: Close encounters (b ~ 3 pc) → 50-90% mass loss
   - **b-dependence**: Distant encounters (b ~ 6 pc) → <10% mass loss
   - **v-dependence**: Slow encounters → gradual tidal stripping
   - **v-dependence**: Fast encounters → impulsive disruption with shocks

4. **Thermal response**:
   - Compression → n increases → T_eq decreases (CNM branch)
   - Cooling maintains thermal equilibrium (for slow encounters)
   - **Shock heating** (for fast encounters): May temporarily exceed T_eq
   - Cooling time vs compression time determines thermal state

### Diagnostics to Track

1. **Morphology**:
   - Axis ratios (a:b:c) from moment of inertia tensor
   - Elongation parameter: e = 1 - c/a
   - **Parameter study**: Plot e(t) for different b and v
   
   **Expected values**:
   - t = 0: e ≈ 0 (spherical cloud, a ≈ b ≈ c)
   - t = t_peri: e ~ 0.5-0.8 (strong tidal elongation, "cigar" shape)
   - t >> t_tidal: e → varies (tidal tails, disrupted structure)

2. **Energetics**:
   - Kinetic energy: E_kin = Σ(½ m v²)
   - Thermal energy: E_th = Σ(m u)
   - Gravitational potential: E_pot (cloud self-gravity + BH)
   - Total energy conservation: ΔE/E₀ < 1%
   - **Parameter study**: E_kin vs time shows velocity damping
   
   **Numerical estimates** (M_cloud = 10⁴ M_☉, v_rel = 10 km/s, T = 50 K):
   ```
   E_kin,initial ≈ ½ M_cloud v_rel² 
                ≈ 0.5 × 10⁴ M_☉ × (10 km/s)²
                ≈ 5×10⁵ M_☉·(km/s)²
                ≈ 10⁵⁴ erg
   
   E_th,initial ≈ M_cloud × u
               ≈ 10⁴ M_☉ × 6.2×10⁹ erg/g × 1.989×10³³ g/M_☉
               ≈ 1.2×10⁴⁷ erg (negligible compared to E_kin!)
   
   E_pot,BH ≈ -G M_BH M_cloud / b
            ≈ -1.43×10⁸ M_☉·(km/s)² (for b=3pc)
            ≈ -2.85×10⁵⁴ erg (dominates total energy!)
   ```

3. **Mass Budget**:
   - Bound mass: M_bound(t) = Σm for particles with E_total < 0
   - Unbound mass: M_unbound(t) = Σm for particles with E_total > 0
   - Accretion rate: dM_BH/dt (particles within R_acc ~ 0.01 pc)
   - **Key metric**: f_disrupted = (M₀ - M_bound(t_final)) / M₀
   
   **Expected disruption fractions**:
   
   | b [pc] | v [km/s] | f_disrupted | M_lost [M_☉] | Physical interpretation |
   |--------|----------|-------------|--------------|------------------------|
   | 3 | 10 | ~0.70 | 7000 | Strong disruption |
   | 4 | 10 | ~0.40 | 4000 | Moderate disruption |
   | 5 | 10 | ~0.20 | 2000 | Weak disruption |
   | 6 | 10 | ~0.10 | 1000 | Minimal disruption |
   | 4 | 5 | ~0.50 | 5000 | Slower → more stripping |
   | 4 | 20 | ~0.25 | 2500 | Faster → less stripping |

4. **Thermal State**:
   - n-T diagram (should follow K&I curve for slow encounters)
   - Departures from equilibrium: ΔT = T - T_eq(n)
   - Cooling time: τ_cool ~ u / |du/dt|
   - **Parameter study**: Fast encounters (high v) show shock heating (T > T_eq)
   
   **Expected temperature evolution** (initial T = 50 K):
   ```
   Compression at pericenter: n increases 10× → ρ → 10× higher
   K&I curve (CNM branch): T_eq ∝ n^(-0.7) (roughly)
   Expected: T_eq drops to ~50 K / 10^0.7 ≈ 10 K (colder!)
   
   Shock heating (v=20 km/s): ΔT_shock ~ (γ-1) v² / (2 k_B / μ m_H)
                                        ~ 0.667 × (20 km/s)² / (2×8.3×10⁷)
                                        ~ 160 K (hot shock!)
   
   Result: Competition between compression cooling vs shock heating
   ```

## References

1. Koyama & Inutsuka (2000), ApJ 532, 980: Thermal equilibrium curves
2. Guillochon & Loeb (2015), ApJ 811, 20: Tidal disruption events
3. Burkert & Naab (2013), MNRAS 434, 36: Cloud-BH interactions
4. Stone et al. (1998), ApJS 114, 345: Numerical methods for tidal disruption
5. Rees (1988), Nature 333, 523: Tidal disruption by supermassive black holes
6. Hills (1975), Nature 254, 295: Tidal capture and disruption

## Appendix: Derivation of Physical Criteria

This appendix provides detailed derivations of all key physical scales and timescales used in the resolution analysis.

---

### A.1 Jeans Mass and Jeans Length

**Physical origin**: Self-gravitating clouds are unstable to collapse when gravity overcomes thermal pressure support.

**Jeans criterion**: For a perturbation of wavelength λ, gravitational collapse occurs when:
```
Gravitational potential energy > Thermal kinetic energy
G M² / R > N k_B T
```

**Derivation from sound wave dispersion**:

Consider a density perturbation: ρ = ρ₀ + δρ e^(i(k·r - ωt))

The dispersion relation for sound waves in self-gravitating medium:
```
ω² = c_s² k² - 4πG ρ₀
```

**Critical wavenumber** (ω² = 0, marginal stability):
```
k_J = √(4πG ρ₀) / c_s
```

**Jeans length**:
```
λ_J = 2π / k_J = c_s √(π / (G ρ₀))
```

**Numerical example** (n_H = 100 cm⁻³, T = 50 K):
```
ρ₀ = 100 cm⁻³ × 1.67×10⁻²⁴ g = 1.67×10⁻²² g/cm³
c_s = √(k_B T / μ m_H) = 0.29 km/s = 2.9×10⁵ cm/s
G = 6.67×10⁻⁸ cm³/(g·s²)

λ_J = 2.9×10⁵ cm/s × √(π / (6.67×10⁻⁸ × 1.67×10⁻²²))
    = 2.9×10⁵ × √(2.81×10¹³)
    ≈ 1.54×10¹⁸ cm
    ≈ 0.50 pc
```

**Jeans mass** (mass within Jeans length sphere):
```
M_J = ρ₀ × (4π/3) × (λ_J/2)³
    = (π/6)^(1/2) × (c_s³ / (G^(3/2) ρ₀^(1/2)))
    = (π^(5/2) / 6) × (c_s³ / (G^(3/2) ρ₀^(1/2)))
```

**Numerical example**:
```
M_J = (π^(5/2) / 6) × (2.9×10⁵)³ / ((6.67×10⁻⁸)^(3/2) × (1.67×10⁻²²)^(1/2))
    ≈ 6.9×10³³ g
    ≈ 3.5 M_☉
```

**SPH resolution requirement**: To resolve fragmentation, need **N_J ≥ 50** particles per M_J.

---

### A.2 Tidal Disruption: Pancaking and Elongation

**Physical setup**: Cloud of mass M_cloud, radius R_cloud, at distance r from BH of mass M_BH.

#### **A.2.1 Tidal Force Components**

The tidal force arises from the **gradient** of the BH gravitational field across the cloud.

**BH acceleration at cloud center** (distance r along x-axis):
```
a_center = -G M_BH / r² x̂
```

**Acceleration at cloud edge** (at position r + δr):
```
a_edge = -G M_BH / (r + δr)² x̂
       ≈ -G M_BH / r² × (1 - 2δr/r) x̂  [Taylor expansion]
```

**Tidal acceleration** (differential force):
```
a_tidal = a_edge - a_center
        = 2 G M_BH δr / r³ x̂
```

**General tidal tensor** (for arbitrary displacement δr = (δx, δy, δz)):
```
T_ij = -G M_BH / r³ × [3 r̂_i r̂_j - δ_ij]
```

where r̂ is the unit vector from BH to cloud center.

**For BH along x-axis**, the tidal tensor is:
```
        ⎡ +2    0    0 ⎤
T = -GM/r³ ⎢  0   -1    0 ⎥
        ⎣  0    0   -1 ⎦
```

#### **A.2.2 Pancaking (Vertical Compression)**

**Perpendicular direction** (y and z): Compression force
```
F_⊥ = m × T_yy × δy = +G M_BH m δy / r³
```

This is a **restoring force** (positive → compression toward orbital plane).

**Equation of motion** (simple harmonic oscillator):
```
d²δy/dt² = (G M_BH / r³) δy
```

**Vertical oscillation frequency**:
```
ω_⊥ = √(G M_BH / r³)
```

**Numerical example** (M_BH = 10⁵ M_☉, r = 3 pc at pericenter):
```
ω_⊥ = √(4.302×10⁻³ × 10⁵ / 3³)  [in units of Myr⁻¹]
    = √(430.2 / 27)
    ≈ 4.0 Myr⁻¹

T_compress = 2π / ω_⊥ ≈ 1.6 Myr
```

#### **A.2.3 Elongation (Radial Stretching)**

**Radial direction** (x): Stretching force
```
F_∥ = m × T_xx × δx = -2 G M_BH m δx / r³
```

This is **anti-restoring** (negative → exponential growth of tidal tails).

**Equation of motion**:
```
d²δx/dt² = -2 (G M_BH / r³) δx
```

**Exponential growth rate**:
```
γ_∥ = √(2 G M_BH / r³) = √2 × ω_⊥
```

**Tidal tail growth**:
```
δx(t) = δx₀ × exp(γ_∥ t)
```

**Numerical example** (same parameters):
```
γ_∥ = √2 × 4.0 ≈ 5.7 Myr⁻¹
e-folding time: τ_stretch = 1/γ_∥ ≈ 0.18 Myr
```

**Physical interpretation**: Material stretches exponentially along the BH direction, forming **leading and trailing tidal tails**.

#### **A.2.4 Tidal Disruption Timescale**

Cloud is **tidally disrupted** when tidal force exceeds self-gravity:
```
G M_BH R_cloud / r³ > G M_cloud / R_cloud²
```

**Critical radius** (tidal radius):
```
r_t = R_cloud × (M_BH / M_cloud)^(1/3)
```

**Disruption timescale** (time for cloud to cross tidal radius):
```
t_tidal ~ r_t / v_orb
        ~ r_t / √(G M_BH / r_t)
        ~ √(r_t³ / (G M_BH))
        ~ √(R_cloud³ × (M_BH / M_cloud) / (G M_BH))
        ~ √(R_cloud³ / (G M_cloud))
```

**Numerical example** (M_cloud = 10⁴ M_☉, R_cloud = 5 pc):
```
t_tidal = √(125 pc³ / (4.302×10⁻³ × 10⁴))
        = √(125 / 43.02)
        ≈ 0.54 Myr
```

**Note**: This is independent of M_BH! Disruption time depends only on cloud properties.

---

### A.3 Tidal Radius vs Hill Radius

These are **different concepts** often confused:

#### **A.3.1 Tidal Radius** (r_t)

**Definition**: Radius at which BH tidal force equals cloud self-gravity.

**Derivation**: Balance tidal force and self-gravity at cloud surface:
```
G M_BH R_cloud / r_t³ = G M_cloud / R_cloud²
```

**Result**:
```
r_t = R_cloud × (M_BH / M_cloud)^(1/3)
```

**Numerical example** (M_cloud = 10⁴ M_☉, M_BH = 10⁵ M_☉, R_cloud = 5 pc):
```
r_t = 5 pc × (10⁵ / 10⁴)^(1/3)
    = 5 pc × 2.15
    ≈ 10.8 pc
```

**Physical meaning**: Cloud passing within r_t will be **tidally disrupted**.

#### **A.3.2 Hill Radius** (r_H or Roche Lobe)

**Definition**: Maximum distance from cloud where cloud's gravity dominates over tidal force.

**Derivation**: Consider test particle at distance r_H from cloud center. It feels:
1. Cloud gravity: F_cloud = G M_cloud / r_H²
2. Tidal force from BH: F_tidal = 2 G M_BH r_H / b³ (where b = impact parameter)

**Balance point** (Hill sphere boundary):
```
G M_cloud / r_H² = 2 G M_BH r_H / b³
```

**Solve for r_H**:
```
r_H³ = M_cloud b³ / (2 M_BH)
r_H = b × (M_cloud / (2 M_BH))^(1/3)
```

**More accurate formula** (including factor of 3 from orbit geometry):
```
r_H = b × (M_cloud / (3 M_BH))^(1/3)
```

**Numerical example** (M_cloud = 10⁴ M_☉, M_BH = 10⁵ M_☉, b = 4 pc):
```
r_H = 4 pc × (10⁴ / (3×10⁵))^(1/3)
    = 4 pc × (0.0333)^(1/3)
    = 4 pc × 0.322
    ≈ 1.29 pc
```

**Physical meaning**: Material beyond r_H will be **tidally stripped** from the cloud.

**Comparison**:
- **r_t**: Where disruption begins (intrinsic to cloud, scales with R_cloud)
- **r_H**: Extent of bound material (depends on impact parameter b)

---

### A.4 Pericenter Passage and Orbital Timescales

#### **A.4.1 Pericenter Time**

For a **parabolic encounter** (typical for cloud-BH interactions):

**Initial conditions**: Cloud at distance b (impact parameter), velocity v_∞ at infinity.

**Energy conservation**:
```
½ v_∞² = ½ v_peri² - G M_BH / r_peri
```

For **hyperbolic orbit** (v_∞ > 0):
```
r_peri ≈ b  (for weak encounters)
```

**Pericenter passage time** (time spent within distance ~ b):
```
t_peri ~ b / v_rel
```

**Numerical examples**:
```
b = 3 pc, v = 10 km/s:
t_peri = 3 pc / 10 km/s
       = 3 × 3.086×10¹⁸ cm / (10 × 10⁵ cm/s)
       = 9.26×10¹² s
       ≈ 0.29 Myr

b = 6 pc, v = 20 km/s:
t_peri = 6 pc / 20 km/s ≈ 0.30 Myr
```

**Physical interpretation**: 
- **Slow encounters** (v small): Long t_peri → gradual tidal stripping
- **Fast encounters** (v large): Short t_peri → impulsive shock

---

### A.5 Smoothing Length and Neighbor Number

#### **A.5.1 SPH Kernel Support**

The SPH smoothing kernel W(r, h) has compact support within radius **2h** (for cubic spline) or **3h** (for Wendland C4).

**Neighbor number**:
```
N_neighbor = ∫ n(r') W(|r - r'|, h) dV
           ≈ n₀ × V_kernel
```

For **Wendland C4 kernel** (support radius = 2h):
```
V_kernel = (4π/3) × (2h)³ = 32π h³ / 3
```

**Average neighbor count**:
```
N_neighbor = n₀ × 32π h³ / 3
```

For **uniform distribution** with N total particles in volume V:
```
n₀ = N / V
```

**Smoothing length** (from neighbor number requirement):
```
N_neighbor = (N / V) × 32π h³ / 3
```

For **spherical cloud** (V = 4πR³/3):
```
N_neighbor = N × (h/R)³ × 8
```

**Solve for h**:
```
h = R × (N_neighbor / (8N))^(1/3)
```

**For N_neighbor = 50** (SPH standard):
```
h = R × (50 / (8N))^(1/3)
  = R × (6.25 / N)^(1/3)
  ≈ R / (N^(1/3) × 0.54)
  ≈ 1.84 R / N^(1/3)
```

**Simplified approximation**:
```
h ≈ R / √N  (slightly overestimates)
```

**Numerical example** (R = 5 pc, N = 2×10⁵):
```
h = 5 pc / √(2×10⁵)
  = 5 / 447
  ≈ 0.011 pc
```

#### **A.5.2 Resolution Requirement**

To resolve a physical scale λ, need:
```
h < λ / 5  (at least 5 smoothing lengths across feature)
```

For **Hill radius** (λ = r_H):
```
h_max = r_H / 5 ≈ 0.2 r_H
```

**Conservative criterion** (used in this work):
```
h ≤ 0.1 r_H
```

This ensures **≥10 smoothing lengths** across the Hill sphere → well-resolved tidal truncation.

---

### A.6 Shock Heating and Compression

#### **A.6.1 Rankine-Hugoniot Shock Conditions**

For a **strong shock** in ideal gas (γ = 5/3):

**Density jump**:
```
ρ₂ / ρ₁ = ((γ+1) M²) / ((γ-1) M² + 2)
```

For **Mach number M → ∞** (strong shock limit):
```
ρ₂ / ρ₁ → (γ+1) / (γ-1) = 4  (for γ = 5/3)
```

**Temperature jump**:
```
T₂ / T₁ = [2γ M² - (γ-1)] × [(γ-1) M² + 2] / [(γ+1)² M²]
```

For **strong shock**:
```
T₂ / T₁ ≈ 2(γ-1) M² / (γ+1)² ≈ 0.4 M²
```

**Numerical example** (v_shock = 20 km/s, c_s,1 = 0.3 km/s):
```
M = v_shock / c_s,1 = 20 / 0.3 ≈ 67

ρ₂ / ρ₁ ≈ 4  (maximum compression)

T₂ / T₁ ≈ 0.4 × 67² ≈ 1800

For initial T₁ = 50 K:
T₂ ≈ 90,000 K  (but cooling brings it down quickly!)
```

#### **A.6.2 Post-Shock Temperature**

**Kinetic energy → thermal energy**:
```
½ v² = (γ / (γ-1)) × (k_B ΔT / μ m_H)
```

**Temperature increase**:
```
ΔT = (γ-1) μ m_H v² / (2 k_B)
```

**Numerical example** (v = 20 km/s, μ = 1):
```
ΔT = 0.667 × 1.67×10⁻²⁴ g × (20×10⁵ cm/s)² / (2 × 1.38×10⁻¹⁶ erg/K)
   ≈ 160 K
```

**But**: Cooling time in CNM is short! Koyama-Inutsuka equilibrium curve wins for **slow encounters**.

---

### A.7 Thermal Equilibrium Timescale

#### **A.7.1 Cooling Function**

From Koyama & Inutsuka (2000), cooling rate per unit volume:
```
Λ(n, T) = n² × L(T)
```

where L(T) is the cooling function [erg·cm³/s].

**Cooling time**:
```
t_cool = thermal energy / cooling rate
       = (3/2) n k_B T / (n² L(T))
       = (3/2) k_B T / (n L(T))
```

**Typical values** (CNM, n = 100 cm⁻³, T = 50 K):
```
L(50 K) ≈ 10⁻²⁶ erg·cm³/s  [from K&I curve]
t_cool ≈ (3/2) × 1.38×10⁻¹⁶ × 50 / (100 × 10⁻²⁶)
       ≈ 10¹⁰ s
       ≈ 0.03 Myr
```

**Compare to tidal timescale**:
```
t_cool / t_tidal ≈ 0.03 / 0.54 ≈ 0.06 = 6%
```

**Conclusion**: Cooling is ~17× faster than tidal disruption → **thermal equilibrium maintained**.

---

### A.8 Summary of Scaling Relations

| Quantity | Scaling | Numerical example |
|----------|---------|-------------------|
| **Jeans mass** | M_J ∝ ρ^(-1/2) T^(3/2) | 3.5 M_☉ (n=100 cm⁻³, T=50K) |
| **Jeans length** | λ_J ∝ ρ^(-1/2) T^(1/2) | 0.5 pc |
| **Tidal radius** | r_t ∝ R (M_BH/M)^(1/3) | 10.8 pc (M=10⁴M_☉, M_BH=10⁵M_☉) |
| **Hill radius** | r_H ∝ b (M/M_BH)^(1/3) | 1.3 pc (b=4pc) |
| **Tidal time** | t_tidal ∝ √(R³/(GM)) | 0.54 Myr (M=10⁴M_☉, R=5pc) |
| **Pericenter time** | t_peri = b / v | 0.29 Myr (b=3pc, v=10km/s) |
| **Cooling time** | t_cool ∝ T / (n L(T)) | 0.03 Myr (n=100 cm⁻³, T=50K) |
| **Smoothing length** | h ∝ R / √N | 0.011 pc (R=5pc, N=2×10⁵) |
| **Shock temperature** | ΔT ∝ v² | 160 K (v=20 km/s) |

**Timescale hierarchy**:
```
t_cool << t_tidal ~ t_peri << t_Hubble
0.03 Myr << 0.5 Myr ~ 0.5 Myr << 10⁴ Myr
```

This justifies the **thermal equilibrium assumption**: clouds cool much faster than they are tidally disrupted.

---

## Next Steps

1. Review existing gravity force implementation in `src/gravity_force.cpp`
2. Implement external point-mass force module
3. Adapt Lane-Emden setup to include K&I thermal equilibrium
4. Create test case with N = 5×10⁴ particles
5. Validate conservation and thermal equilibrium before production runs
