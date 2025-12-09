# GSPH Grad-h Instability: Comprehensive Analysis

## Executive Summary

This document provides a comprehensive analysis of why GSPH (Godunov SPH) without the grad-h correction experiences catastrophic gravitational collapse in stratified self-gravitating systems, while Standard SPH (SSPH) remains stable.

**Key Finding**: The instability is a **Variational Inconsistency Instability** (also called "Pressure Deficit Feedback Instability") - a purely numerical artifact arising from incomplete variational derivative treatment in the SPH Lagrangian formulation.

**Practical Recommendation**: **Always use the grad-h correction (useGradH = true) for self-gravitating systems with density gradients.**

---

## 1. Theoretical Foundation

### 1.1 SPH Density Estimation

The SPH density at particle i is:

$$\rho_i = \sum_j m_j W(|\mathbf{r}_i - \mathbf{r}_j|, h_i)$$

### 1.2 Adaptive Smoothing Length

The smoothing length is determined self-consistently:

$$h_i = \eta \left(\frac{m_i}{\rho_i}\right)^{1/D}$$

This creates an implicit relation where $\rho_i$ depends on $h_i$, and $h_i$ depends on $\rho_i$.

---

## 2. The Grad-h Correction Factor Ω

### 2.1 Derivation

From the SPH Lagrangian and chain rule with h(ρ) dependence:

$$\boxed{\Omega_k = 1 + \frac{h_k}{D\rho_k}\sum_j m_j \frac{\partial W_{kj}}{\partial h_k}}$$

### 2.2 Correct Pressure Force

The variationally consistent pressure force is:

$$\mathbf{a}_i = -\sum_j m_j \left[\frac{P_i}{\Omega_i \rho_i^2}\nabla W_{ij}(h_i) + \frac{P_j}{\Omega_j \rho_j^2}\nabla W_{ij}(h_j)\right]$$

### 2.3 Ω Values in Stratified Media

| Location | ∇ρ | Ω Value |
|----------|-----|---------|
| Core (uniform) | ≈ 0 | → 1.0 |
| Near surface (steep gradient) | large | → 1.1-1.3 |

---

## 3. Pressure Force Error

### 3.1 Error Without Grad-h

The fractional error is:

$$\epsilon = 1 - \frac{1}{\Omega} = \frac{\Omega - 1}{\Omega}$$

For typical stratified equilibria: **ε ≈ 5% - 13% UNDERESTIMATE**

### 3.2 Net Force Imbalance

At hydrostatic equilibrium: $\mathbf{a}_P^{\text{true}} + \mathbf{a}_g = 0$

Without grad-h:
$$\mathbf{a}_{\text{net}} = \epsilon \times \mathbf{a}_g$$

Since gravity points inward and ε > 0: **NET INWARD ACCELERATION**

---

## 4. The Instability Mechanism

### 4.1 Positive Feedback Loop

```
1. Initial state: stratified equilibrium (∇ρ ≠ 0)
         ↓
2. Pressure underestimate: ε ∝ |∇ρ| > 0
         ↓
3. Net inward force: a_net = ε × a_g
         ↓
4. Material flows inward → Density increases
         ↓
5. Gradient steepens → Error grows
         ↓
6. [POSITIVE FEEDBACK → RUNAWAY COLLAPSE]
```

### 4.2 Observed Collapse Dynamics

The collapse follows the **delayed singularity** form:

$$\rho(t) = \rho_{bg} + \frac{A}{(t_c - t)^2}$$

| Time | ρ_max | Phase |
|------|-------|-------|
| 0.0 - 7.4 | 1.4 - 2.1 | **Quasi-stable** (oscillating) |
| 7.4 - 8.4 | 2.0 - 4.2 | **Acceleration** |
| 8.4 - 8.8 | 4.2 - 18.4 | **Runaway** |
| 8.8 - 9.0 | 18.4 - 4521 | **Singularity approach** |

**Key insight**: The collapse is **SUDDEN, not gradual** - system appears stable for ~7 time units, then collapses in ~1.5 time units.

---

## 5. Why SSPH Survives

### 5.1 SSPH Pressure Formulation

SSPH uses symmetric pressure averaging:
$$\mathbf{a}_i^{\text{SSPH}} = -\sum_j m_j \frac{P_i + P_j}{2\rho_i\rho_j}\nabla W_{ij}$$

### 5.2 Error Cancellation

The averaging cancels systematic biases - if both $P_i$ and nearby $P_j$ are underestimated similarly, the force **difference** is preserved.

### 5.3 GSPH Vulnerability

GSPH's single Riemann pressure $p^*$ directly propagates density estimation errors without compensating averaging.

---

## 6. Simulation Results Summary

| Method | Grad-h | Final ρ_max | Status |
|--------|--------|-------------|--------|
| GSPH | ON | 2.03 | **STABLE** |
| GSPH | OFF | 4521 | **COLLAPSED** |
| SSPH | ON | 2.00 | **STABLE** |
| SSPH | OFF | 2.18 | **STABLE** |

---

## 7. Classification

**Name**: Pressure Deficit Induced Singularity (PDIS) / Variational Inconsistency Instability

**Characteristics**:
- ✓ NUMERICAL (not physical) - disappears with grad-h correction
- ✓ DELAYED (long quasi-stable phase before collapse)
- ✓ SINGULAR (approaches $\rho \to \infty$ at finite time $t_c$)
- ✓ SCHEME-DEPENDENT (affects GSPH but not SSPH)

**What it is NOT**:
- ✗ Jeans instability (physical, occurs above $M_J$)
- ✗ Tensile instability (creates particle clumping patterns)
- ✗ Pairing instability (related to kernel properties)

---

## 8. The Cure

### Always use the grad-h correction:

$$\mathbf{a}_i = -\sum_j m_j \left[\frac{P_i}{\Omega_i\rho_i^2}\nabla W_{ij}(h_i) + \frac{P_j}{\Omega_j\rho_j^2}\nabla W_{ij}(h_j)\right]$$

**Benefits**:
1. ✓ Exact energy and momentum conservation
2. ✓ Correct pressure forces in stratified media
3. ✓ Variational consistency of the SPH Lagrangian
4. ✓ Stable hydrostatic equilibria
5. ✓ Minimal computational overhead (<5%)

---

## References

1. Springel & Hernquist (2002) - Cosmological SPH with variable h
2. Hopkins (2013) - GIZMO: A new class of accurate hydrodynamics methods
3. Rosswog (2009) - Astrophysical SPH for relativistic systems
4. Price & Monaghan (2007) - On the variable smoothing length in SPH

---

*Document consolidated from first-principles analysis, December 2024*
