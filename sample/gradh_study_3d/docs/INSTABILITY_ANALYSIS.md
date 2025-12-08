# GSPH Grad-H Instability: Complete Analysis

## Executive Summary

GSPH (Godunov SPH) without the grad-h correction (Ω factor) exhibits a **secular numerical drift** that leads to core collapse in self-gravitating systems like the Lane-Emden sphere. This document presents a data-driven derivation of the instability mechanism.

## Key Findings

| Quantity | Value |
|----------|-------|
| Measured drift rate Γ | 0.020 code units⁻¹ |
| Free-fall frequency ω_ff | 4.25 code units⁻¹ |
| Ratio Γ/ω_ff | 0.005 |
| Collapse timescale | ~50 code units ≈ 110 × t_ff |
| Effective force coefficient C_eff | 0.1 - 0.7 (decreasing) |

## 1. The Instability Mechanism

### 1.1 Why Standard SPH is Stable (C = 1)

In Standard SPH, the pressure force on particle i is:

$$\mathbf{F}_i^P = -\sum_j m_j \left( \frac{P_i}{\rho_i^2} \nabla_i W_{ij} + \frac{P_j}{\rho_j^2} \nabla_j W_{ij} \right)$$

Each term uses the **LOCAL** pressure. Under uniform compression with δρ:
- P responds: δP/P = γ·δρ/ρ
- h changes: h ∝ ρ^(-1/3), so |∇W| ∝ ρ^(4/3)
- Net force response: δF/F = (γ - 2/3)·δρ/ρ = δρ/ρ for γ = 5/3

**Result: Force coefficient C = 1, matching gravity's response.**

### 1.2 Why GSPH Without Ω is Unstable (C < 1)

In GSPH, the pressure force uses the **shared Riemann solver pressure** p*:

$$\mathbf{F}_i^P = -\sum_j m_j \, p^*_{ij} \left( \frac{1}{\rho_i^2} \nabla_i W_{ij} + \frac{1}{\rho_j^2} \nabla_j W_{ij} \right)$$

The Riemann pressure **AVERAGES** left and right states:

$$p^* \approx \frac{c_R P_L + c_L P_R}{c_L + c_R}$$

### 1.3 The Non-Uniform Perturbation Effect

From simulation data at t=50:

| Region | δρ/ρ₀ | P | Comment |
|--------|-------|---|---------|
| Center (r<0.1) | +0.95 | 1.30 | Strong compression |
| Edge (r~0.5) | -0.16 | 0.17 | Slight expansion |

The Riemann pressure for a center-edge pair:
- P_center = 1.30
- P_edge = 0.17
- **p* = 0.68** (average!)

The GSPH force on the center particle uses p* = 0.68 instead of P_center = 1.30.

**The force is reduced by factor p*/P_center = 0.52!**

### 1.4 The Effective Force Coefficient

Measured from simulation:

```
snap   ρ/ρ₀    F/F₀    C_eff = (δF/F)/(δρ/ρ)
--------------------------------------------
  10   1.32   1.19    0.58
  30   1.57   1.09    0.16
  50   2.02   1.13    0.13
  70   2.65   1.12    0.07
```

**C_eff starts around 0.6 and decreases toward 0.1 as compression increases.**

## 2. The Secular Drift

### 2.1 Why Not Exponential Instability?

Simple theory predicts: Γ ≈ √(1-C_eff) × ω_ff ≈ 2.8

But measured: Γ ≈ 0.02

**Resolution: The system starts in equilibrium!**

The Lane-Emden sphere is initially in hydrostatic balance. The weak force coefficient means:
1. Equilibrium is **approximate** but not exact
2. Small imbalance causes slow **drift**, not rapid collapse
3. Growth is **secular** (polynomial) not exponential

### 2.2 The Drift Rate

The force-to-gravity ratio ε = F_P/ρ relative to initial:

```
snap   ρ/ρ₀    ε = (F_P/ρ)/(F₀/ρ₀) - 1
------------------------------------------
   0   1.00    0.00
  25   1.51   -0.26
  50   2.02   -0.44
  75   2.84   -0.61
```

ε becomes increasingly negative → pressure force increasingly inadequate.

**Drift rate: d ln(ρ)/dt ≈ 0.02, giving τ_drift ≈ 50 code units ≈ 110 × t_ff**

## 3. How Ω Fixes the Instability

### 3.1 The Ω Factor

The grad-h correction Ω is defined as:

$$\Omega = \frac{1}{1 + \frac{h}{D\rho}\frac{d\rho}{dh}}$$

where dρ/dh is computed from the kernel derivative sum.

### 3.2 Force with Ω

$$\mathbf{F}_i^P = -\sum_j m_j \, p^*_{ij} \left( \frac{\Omega_i}{\rho_i^2} \nabla_i W_{ij} + \frac{\Omega_j}{\rho_j^2} \nabla_j W_{ij} \right)$$

At compressed regions:
- h decreases, ρ increases
- dρ/dh < 0 → Ω > 1
- **Force is BOOSTED by factor Ω > 1**

From data: Center Ω ≈ 1.04-1.13 when compression occurs.

### 3.3 Variational Consistency

More fundamentally, Ω ensures:
1. **Discrete energy conservation**: dE/dt = 0 exactly
2. **Exact equilibrium**: Hydrostatic balance is a discrete fixed point
3. **No secular drift**: System oscillates around equilibrium

## 4. Stability Comparison

### With Ω (grad-h enabled):
- Core density: ρ₀ → 1.66ρ₀ (stable oscillation)
- Compression ratio: 1.15x
- No runaway collapse

### Without Ω (grad-h disabled):
- Core density: ρ₀ → 22ρ₀ (collapse!)
- Compression ratio: 22x
- Secular drift leads to runaway

## 5. Theoretical Framework

### 5.1 Stability Criterion

For a self-gravitating system, stability requires:

$$\frac{\delta F_P}{F_P} \geq \frac{\delta F_g}{F_g}$$

With gravity giving C_g = 1, we need C_P ≥ 1.

### 5.2 GSPH Force Response

The effective coefficient in GSPH without Ω:

$$C_{\text{eff}} = R \cdot \gamma - \frac{2}{3}$$

where R = (δp*/p*) / (δP/P) is the Riemann pressure response reduction.

For R ≈ 0.77 (from data): C_eff ≈ 0.62 < 1 → **UNSTABLE**

### 5.3 Drift Rate

The secular drift rate:

$$\Gamma_{\text{drift}} \approx 0.02 \ll \omega_{ff} \approx 4.2$$

This slow rate explains the quasi-static nature of the collapse.

## 6. Conclusions

1. **Mechanism**: GSPH without Ω has weak pressure response (C < 1) due to averaged Riemann pressure

2. **Nature**: Secular drift, not dynamical instability; timescale ~100 × t_ff

3. **Fix**: The Ω factor restores C = 1 through variational consistency and discrete energy conservation

4. **Practical implication**: Always use grad-h correction for self-gravitating GSPH simulations

## References

- Springel & Hernquist (2002): Grad-h formulation for GSPH
- Hopkins (2013): GIZMO and modern SPH formulations
- Price (2012): SPH review - smoothing length and grad-h terms
