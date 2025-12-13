# Resolution Analysis for BNS Cocoon Shock Breakout Simulations

Based on Gutiérrez et al. (2024) arXiv:2408.15973v3

## Physical Scales to Resolve

### 1. Ejecta Structure

From NR simulations (Radice+2018), the ejecta has:
- **v_min** ~ 0.05c → r_min = v_min × t₀
- **v_max** ~ 0.6–0.8c → r_max = v_max × t₀
- **Radial extent**: Δr = r_max - r_min ≈ 0.55–0.75 × c × t₀

At reference time t₀ = 1 s (typical for GRB 170817A analysis):
- r_min ≈ 0.05c × 1s ≈ 1.5 × 10⁹ cm
- r_max ≈ 0.7c × 1s ≈ 2.1 × 10¹⁰ cm

### 2. Cocoon/Jet Structure

From paper:
- **Jet opening angle**: θ_j,0 ~ 8–16° ≈ 0.14–0.28 rad
- **Cocoon width**: θ_c ~ 20–40° at breakout
- **Jet head size**: r_h ~ c × t / Γ_h² (relativistic beaming)

The cocoon is formed by jet-ejecta interaction. Key scales:
- **Cocoon thickness**: δr_c ~ r_c / 10 at breakout
- **Shock width**: δ_shock ~ few × smoothing length

### 3. Density Contrast

From paper Table 1 (ejecta profiles):
- Peak density (inner ejecta): ρ_peak ~ 10⁻⁸ g/cm³
- Outer ejecta density: ρ_outer ~ 10⁻¹² g/cm³
- Density contrast: ~10⁴

### 4. Extended Tail (if used)

For extended tail profile:
- ρ(v > v_break) ∝ (Γβ)^(-α) with α ~ 5–10
- v_break ~ 0.5c
- Tail extends to v ~ 0.9c in some models

---

## Resolution Requirements

### Minimum Resolution Criteria

1. **Radial resolution**: Need ≥ 5–10 particles across shock width
   - Shock width ~ 3h (smoothing length)
   - Ejecta radial extent ~ 50–100 × h

2. **Angular resolution**: Need ≥ 3–5 particles across jet angle
   - For θ_j = 15°, need Δθ ≤ 3–5°
   - n_angular ≥ 180° / 5° = 36

3. **Density contrast**: Need enough dynamic range
   - With 10⁴ contrast, need mass resolution ~ 10⁻⁵ × M_ej

### Particle Count Estimates

#### 2D Axisymmetric (r, θ plane)

```
N_particles ≈ n_radial × n_angular × 2 (both hemispheres)
```

| Resolution | n_radial | n_angular | N_total | Purpose |
|-----------|----------|-----------|---------|---------|
| **Test** | 50 | 40 | ~4,000 | Quick tests, debugging |
| **Low** | 100 | 80 | ~16,000 | Proof of concept |
| **Medium** | 200 | 160 | ~64,000 | Publication quality (2D) |
| **High** | 400 | 320 | ~256,000 | Convergence study |
| **Production** | 800 | 640 | ~1,000,000 | Full resolution |

#### 3D Full Simulation

```
N_particles ≈ (4/3)π × n_r³ × fill_fraction
```

For 3D with (r, θ, φ):

| Resolution | n_r | n_θ | n_φ | N_total | Notes |
|-----------|-----|-----|-----|---------|-------|
| **Minimal** | 50 | 30 | 60 | ~90,000 | Barely adequate |
| **Low** | 100 | 60 | 120 | ~720,000 | Basic structure |
| **Medium** | 200 | 120 | 240 | ~5.8M | Publication |
| **High** | 400 | 240 | 480 | ~46M | High fidelity |

---

## Critical Phenomena to Resolve

### 1. Shock Breakout Timescale

From paper Eq. (11)–(15):
```
t_bo ≈ r_bo / (v_s × Γ_c²)
```

Need temporal resolution:
- Δt < t_bo / 100 for accurate breakout detection
- Typically Δt ~ 10⁻³–10⁻² seconds

SPH: CFL condition gives Δt ~ h / c_s, so:
- Need h small enough that CFL timestep resolves breakout

### 2. Optical Depth Unity Surface (τ = 1)

From paper Eq. (6):
```
τ = κ × ∫ρ dr
```

With κ ~ 0.2 cm²/g (Thomson):
- Need to resolve where τ crosses unity
- Requires ≥ 5 particles in the τ ~ 1 layer

### 3. Cocoon Lorentz Factor at Breakout

From paper:
```
Γ_c,bo ~ 2–5 for successful cocoon breakout
```

Need to resolve:
- Internal energy distribution
- Pressure gradient driving expansion
- Requires accurate EOS (γ = 4/3)

---

## Recommended Configurations

### For GRB 170817A-like Events

Based on paper parameters:
- t_GRB = 1.73 s (observed delay)
- E_iso ~ 5 × 10⁴⁶ erg
- θ_j ~ 15°

**Minimum viable (2D)**: 
```json
{
  "n_radial": 150,
  "n_angular": 120,
  "N_total": ~36,000
}
```

**Production (2D)**:
```json
{
  "n_radial": 400,
  "n_angular": 300,
  "N_total": ~240,000
}
```

### Memory and Runtime Estimates

| N_particles | Memory (2D) | Time (100 steps) | Time (full) |
|-------------|------------|------------------|-------------|
| 4,000 | ~50 MB | ~10 s | ~10 min |
| 16,000 | ~200 MB | ~1 min | ~1 hr |
| 64,000 | ~800 MB | ~10 min | ~10 hr |
| 256,000 | ~3 GB | ~1 hr | ~4 days |

---

## Physics Parameters Summary (Paper Values)

### Ejecta (NR-derived)

| Model | v_max (c) | M_ej,iso (M☉) | Profile |
|-------|-----------|---------------|---------|
| DD2_M135-135 | 0.73 | 0.0032 | Extended |
| DD2_M180-108 | 0.68 | 0.0023 | Extended |
| BLh_M1146-1635 | 0.78 | 0.0076 | Extended |
| SFHo_M135-135 | 0.73 | 0.0080 | Extended |
| SLy_M145-125 | 0.58 | 0.0104 | Extended |

### Engine (magnetar/SGRB)

| Parameter | Range | Fiducial |
|-----------|-------|----------|
| B (G) | 5×10¹⁴ – 10¹⁶ | 10¹⁵ |
| t_del (s) | 0.1 – 2 | 0.5 |
| θ_j,0 | 8° – 16° | 12° |

### Cocoon at Breakout

| Parameter | Value |
|-----------|-------|
| Γ_c,bo | 1 – 5 |
| T_obs (keV) | ~50 × Γ_bo |
| t_bo (s) | 0.5 – 2 |

---

## Code Units Conversion

For SRGSPH simulations, typical normalization:
- Length: r₀ = c × t₀ = 3 × 10¹⁰ cm (for t₀ = 1 s)
- Mass: M₀ = M☉ = 2 × 10³³ g
- Time: t₀ = 1 s
- Velocity: c = 1 (natural units)
- Density: ρ₀ = M₀ / r₀³ = 7.4 × 10⁻²¹ g/cm³

Physical to code:
```
M_ej = 0.01 (code) = 0.01 M☉
v_max = 0.7 (code) = 0.7c
ρ_ejecta ~ 10⁻⁸ g/cm³ = 10¹³ ρ₀
```

---

## References

1. Gutiérrez et al. (2024) arXiv:2408.15973v3 - Cocoon shock breakout
2. Radice et al. (2018) ApJL 869, L35 - NR ejecta profiles
3. Shibata & Hotokezaka (2019) ARNPS 69, 41 - GW170817 review
4. Kasliwal et al. (2017) Science 358, 1559 - AT2017gfo observations
