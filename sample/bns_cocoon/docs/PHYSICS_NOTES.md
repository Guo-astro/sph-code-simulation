# Physics Notes: Cocoon Shock Breakout

Detailed physics derivations for cocoon shock breakout simulations.

## 1. Ejecta Structure from GR Simulations

### 1.1 Homologous Expansion

After the merger (t > few ms), the dynamical ejecta expand approximately homologously:

```
r(t) = v × t
```

At reference time t₀ (typically 1 s), a mass element with asymptotic velocity v is located at:

```
r = v × t₀ = v × c × t₀/c = β × c × t₀
```

### 1.2 Velocity-to-Density Conversion

Given the differential mass distribution dM/dv from GR simulations:

```
ρ(r) × 4πr² dr = dM/dv × dv
```

Since dr = t₀ dv (homologous expansion):

```
ρ(r) = (dM/dv) / (4πr² t₀)
```

In code units (length L, time L/c):

```
r_code = v × c × t₀ / L
ρ_code = (dM/dv) / (4π r²_code × t₀ × c)
```

### 1.3 Analytic Velocity Distributions

**Hotokezaka et al. (2013, 2015):**

Kinetic energy distribution:
```
E(≥β) ∝ β^(-0.5)
```

This implies:
```
dE/dβ ∝ β^(-1.5)
```

For non-relativistic ejecta (E ≈ ½Mv²):
```
dM/dβ ∝ β^(-2.5)  [approximately]
```

We use:
```
dM/dv ∝ v^(-2) × exp(-v/v_max)
```

with v_max ≈ 0.4-0.5 c and average velocity ⟨v⟩ ≈ 0.2 c.

**Bauswein et al. (2013):**

Angular concentration near equator:
```
ρ(θ) ∝ exp[-(θ - π/2)² / (2σ²)]
```

with σ ≈ 30° (ejecta within ~60° of orbital plane).

## 2. Shock Dynamics

### 2.1 Shock Jump Conditions (Relativistic)

For a strong relativistic shock:

**Compression ratio:**
```
ρ₂/ρ₁ = (γ̂ + 1) / (γ̂ - 1) × Γ_sh + 1
```

where γ̂ is the adiabatic index and Γ_sh is the shock Lorentz factor.

For ultra-relativistic shock (Γ_sh >> 1):
```
ρ₂/ρ₁ ≈ 4Γ_sh  [for γ̂ = 4/3]
```

**Post-shock pressure:**
```
P₂ ≈ (2/3) ρ₁ c² Γ_sh²
```

**Internal energy density:**
```
e₂ ≈ ρ₁ c² Γ_sh²
```

### 2.2 Shock Propagation

The shock velocity in the lab frame:

```
v_sh = dr_sh/dt
Γ_sh = 1/√(1 - v_sh²/c²)
```

For a blast wave in a power-law density profile ρ ∝ r^(-k):

```
Γ_sh ∝ r^(-(3-k)/2)    [Sedov-Taylor analog]
```

## 3. Optical Depth and Breakout

### 3.1 Optical Depth

```
τ(r) = ∫_r^∞ κ ρ(r') dr'
```

where κ is the opacity (typically 0.1-1 cm²/g for electron scattering in merger ejecta).

### 3.2 Breakout Condition

**Non-relativistic shock:**
```
τ(r_bo) = 1
```

**Relativistic shock:**
```
τ(r_bo) = c/v_sh
```

At breakout, photons can escape on the shock crossing timescale.

### 3.3 Breakout Radius

For power-law density profile:
```
r_bo ≈ (κ M_ej / 4π)^(1/2)
```

For typical parameters:
```
r_bo ~ 10^{11} cm × (M_ej/0.01 M_☉)^(1/2) × (κ/0.2 cm²g⁻¹)^(1/2)
```

## 4. Breakout Observables

### 4.1 Internal Energy at Breakout

Energy in the shocked shell:
```
E_int = ∫_shell ρ u dV ≈ ⟨ρ u⟩ V_shell
```

For strong shock:
```
u ≈ Γ_sh c² / (γ̂ - 1)
```

### 4.2 Comoving Temperature

Radiation energy density in the shocked region:
```
e_rad ≈ E_int / V_shell
```

Assuming thermal equilibrium:
```
e_rad = a T_co⁴
T_co = (e_rad / a)^(1/4)
```

where a = 4σ/c is the radiation constant.

### 4.3 Observed Temperature

Doppler boosting for emission toward observer:
```
T_obs = δ T_co / (1 + z)
```

where δ ≈ Γ_bo (for head-on emission) is the Doppler factor.

For typical breakout:
```
T_obs ~ 10-100 keV
```

### 4.4 Flash Timescales

**Curvature timescale** (angular spreading):
```
t_curv = R_bo / (2c Γ_bo²)
```

**Shell crossing time:**
```
t_shell = ΔR / c ≈ R_bo / (c Γ_bo)
```

**Flash duration:**
```
t_flash = max(t_curv, t_shell)
```

For relativistic breakout (Γ_bo >> 1):
```
t_flash ≈ t_curv ∝ Γ_bo^(-2)
```

### 4.5 Radiated Energy and Luminosity

Radiated energy:
```
E_rad = f_rad × E_int
```

where f_rad ≈ 0.1 is the radiation efficiency.

Isotropic equivalent:
```
E_iso = (4π / ΔΩ) × E_rad
```

where ΔΩ is the solid angle of emission.

Peak luminosity:
```
L_peak = E_rad / t_flash
L_iso = (4π / ΔΩ) × L_peak
```

## 5. Application to GRB 170817A

### 5.1 Observed Properties

| Parameter | Value | Reference |
|-----------|-------|-----------|
| E_iso | (4.6 ± 1.5) × 10⁴⁶ erg | Goldstein+2017 |
| E_peak | 185 ± 62 keV | Goldstein+2017 |
| T₉₀ | 2.0 ± 0.5 s | Abbott+2017 |
| Distance | 40 Mpc | Abbott+2017 |
| Viewing angle | ~20° | Mooley+2018 |

### 5.2 Inferred Cocoon Properties

From cocoon breakout models:

```
E_cocoon ~ 10^{49-50} erg
Γ_bo ~ 2-5
θ_open ~ 20-40°
```

### 5.3 Ejecta Constraints from Kilonova

| Parameter | Range | Reference |
|-----------|-------|-----------|
| M_ej (total) | 0.03-0.05 M_☉ | Villar+2017 |
| M_ej (blue) | 0.01-0.02 M_☉ | Cowperthwaite+2017 |
| v_ej | 0.1-0.3 c | Mooley+2018 |
| v_fast | up to 0.6-0.8 c | Hotokezaka+2018 |

## 6. Numerical Implementation Notes

### 6.1 SRGSPH Requirements

- Gaussian kernel (mandatory for GSPH)
- Relativistic energy equation
- Baryon number tracking (nu field)
- Proper velocity handling (v/c < 1)

### 6.2 Resolution Requirements

For accurate shock capturing:
- Resolve shock width: Δx < R_bo / Γ_bo
- Smoothing length: h ~ 30 neighbors
- Time step: Δt < CFL × Δx / c

### 6.3 Unit System

Recommended relativistic units:
- Length: L = 10^10 cm (or c × t₀)
- Time: L/c
- Velocity: c
- Density: arbitrary (set by M_ej)

## References

1. Nakar, E. & Piran, T. (2017), ApJ 834:28 - Cocoon emission
2. Gottlieb, O. et al. (2018), MNRAS 479:588 - Jet-cocoon structure
3. Hotokezaka, K. et al. (2013), PRD 87:024001 - Dynamical ejecta
4. Radice, D. et al. (2018), ApJ 869:130 - Ejecta dataset
5. Nakar, E. & Sari, R. (2012), ApJ 747:88 - Shock breakout theory
