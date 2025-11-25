# Special Relativistic Godunov SPH: Complete Mathematical Formulation

This document provides a comprehensive compilation of all formulas for Special Relativistic Godunov Smoothed Particle Hydrodynamics (SRGSPH), including detailed Riemann solver formulations. **All formulations follow Kitajima et al. (SRGSPH) as the primary source**.

**Reference Papers:**
- **[SRGSPH]**: Kitajima, K., Inutsuka, S., & Seno, I. (2025). Special Relativistic Smoothed Particle Hydrodynamics Based on Riemann Solver.
- **[Pons]**: Pons, J.A., Martí, J.M., & Müller, E. (2000). The exact solution of the Riemann problem with non-zero tangential velocities in relativistic hydrodynamics.

**Important Note on Units and Notation:**
- SRGSPH keeps factors of $c$ explicit in theoretical formulations, setting $c=1$ only for numerical calculations
- Pons sets $c=1$ throughout and uses $h$ (specific enthalpy) instead of $H$ (enthalpy per baryon)
- This document follows SRGSPH's approach: explicit $c$ in theory, numerical calculations with $c=1$

---

## Table of Contents

1. [Notation and Conventions](#notation-and-conventions)
2. [Lagrangian Derivation of Equations of Motion](#lagrangian-derivation-of-equations-of-motion)
   - [Non-Relativistic Lagrangian Formulation](#non-relativistic-lagrangian-formulation)
   - [Relativistic Lagrangian Formulation](#relativistic-lagrangian-formulation)
3. [Theory of Relativistic Hydrodynamics](#theory-of-relativistic-hydrodynamics)
4. [Basic Equations (SRGSPH Formulation)](#basic-equations-srgsph-formulation)
5. [Conservative Formulation](#conservative-formulation)
6. [SRGSPH Formulation with Fixed Smoothing Length](#srgsph-formulation-with-fixed-smoothing-length)
7. [SRGSPH Formulation with Variable Smoothing Length](#srgsph-formulation-with-variable-smoothing-length)
8. [Riemann Problem Theory](#riemann-problem-theory)
9. [Riemann Problem: Rarefaction Waves](#riemann-problem-rarefaction-waves)
10. [Riemann Problem: Shock Waves](#riemann-problem-shock-waves)
11. [Complete Riemann Solver Algorithm](#complete-riemann-solver-algorithm)
12. [Primitive Variable Recovery](#primitive-variable-recovery)
13. [Time Integration](#time-integration)
14. [Numerical Implementation (c=1)](#numerical-implementation-c1)

---

## Notation and Conventions

### Units and Speed of Light
- **Theory**: Speed of light $c$ kept explicit
- **Numerical calculations**: $c = 1$ **[SRGSPH, numerical implementation]**
- Minkowski metric: $\eta^{\mu\nu} = \text{diag}(-1, 1, 1, 1)$

### SRGSPH Symbols (Kitajima et al.)

**Spacetime coordinates:**
- $t$: time
- $\mathbf{x} = (x, y, z)$ or $\bm{x}$: spatial position vector
- $x^\mu = (t, x, y, z)$: spacetime coordinates ($\mu = 0, 1, 2, 3$)

**Primitive variables (physical quantities):**
- $n$: baryon number density in rest frame **[SRGSPH Eq. 114]**
- $N$: baryon number density in lab frame **[SRGSPH Eq. 91]**
- $P$: pressure **[SRGSPH Eq. 92]**
- $\bm{v} = (v^x, v^y, v^z)$: three-velocity **[SRGSPH]**
- $v^n$: normal velocity component
- $v^t = \sqrt{(v^y)^2 + (v^z)^2}$: tangential velocity magnitude
- $u$: thermal energy per baryon **[SRGSPH Eq. 120]**
- $\gamma_c$: ratio of specific heats (adiabatic index) **[SRGSPH Eq. 119]**
- $c_s$: sound speed

**Relativistic quantities:**
- $\gamma$: Lorentz factor **[SRGSPH Eq. 111]**
- $H$: enthalpy per baryon **[SRGSPH Eq. 112]**
- $\mathbf{S}$ or $\bm{S}$: relativistic canonical momentum per baryon **[SRGSPH Eq. 104]**
- $e$: canonical energy per baryon in lab frame **[SRGSPH Eq. 107]**
- $u^\mu$: four-velocity

**Conserved variables (grid-based formulation):**
- $D$: conserved rest-mass density
- $S^i$: conserved momentum density
- $\tau$: conserved energy density

**SPH-specific quantities:**
- $\nu$: baryon number in an SPH particle **[SRGSPH Eq. 134]**
- $W(\bm{x}, h)$: kernel function **[SRGSPH Eq. 138]**
- $h$ or $h(\bm{x})$: smoothing length
- $V_p$ or $V_{\rm p}$: particle volume **[SRGSPH Eq. 221]**
- $\langle f_i \rangle$: convolved quantity for particle $i$ **[SRGSPH Eq. 148]**
- $d$: number of spatial dimensions
- $\eta$: smoothing length parameter **[SRGSPH Eq. 231]**
- $C_{\rm smooth}$: smoothing coefficient **[SRGSPH Eq. 233]**

**Riemann problem notation:**
- $L, R$: left and right initial states
- $L_*, R_*$: left and right intermediate states
- $P_*$: pressure in intermediate state
- $v^x_*$: normal velocity in intermediate state
- $a, b$: states ahead and behind a wave
- $\xi$: self-similarity variable ($x/t$)
- $j$: invariant mass flux across shock
- $V_s$: shock velocity
- $\gamma_s$: shock Lorentz factor

### Relationship to Pons et al. Notation

| Quantity | SRGSPH (this doc) | Pons et al. |
|----------|-------------------|-------------|
| Lorentz factor | $\gamma$ | $W$ |
| Enthalpy per baryon | $H$ (with explicit $c^2$) | $h$ (with $c=1$) |
| Pressure | $P$ | $p$ |
| Adiabatic index | $\gamma_c$ | $\gamma$ |
| Baryon density (rest) | $n$ | $\rho/m_b$ or absorbed into $\rho$ |
| Thermal energy/baryon | $u$ | $\epsilon \rho/n$ or absorbed |

**Key difference**: SRGSPH uses enthalpy **per baryon** $H$ with baryon number density $n$, while Pons uses specific enthalpy $h$ with mass density $\rho$.

When $c=1$ and using natural units: $H \approx h$ numerically, but the formulations differ in how they're derived.

---

## Lagrangian Derivation of Equations of Motion

This section provides the fundamental Lagrangian derivations for both non-relativistic and special relativistic hydrodynamics, showing how the equations of motion used in the simulation arise from first principles.

### Non-Relativistic Lagrangian Formulation

#### Action Principle

The non-relativistic hydrodynamics equations can be derived from a variational principle applied to the action:

```math
S = \int dt \int d^3x \, \mathcal{L}(\rho, \bm{v}, \nabla\rho, \nabla\bm{v})
```

where $\mathcal{L}$ is the Lagrangian density.

#### Lagrangian Density for Ideal Fluid

For a non-relativistic ideal fluid, the Lagrangian density is:

```math
\mathcal{L} = \rho \left[ \frac{1}{2} \bm{v}^2 - u(\rho, s) \right]
```

where:
- $\rho$ is the mass density
- $\bm{v}$ is the velocity field
- $u(\rho, s)$ is the specific internal energy
- $s$ is the specific entropy (constant in isentropic flow)

#### Euler-Lagrange Equations

Using the Lagrangian formulation with material (Lagrangian) coordinates $\bm{a}$ and spatial (Eulerian) coordinates $\bm{x}(\bm{a}, t)$, the equations of motion are:

**Continuity Equation**:
```math
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \bm{v}) = 0
```

**Euler Equation (Momentum Conservation)**:
```math
\rho \frac{d\bm{v}}{dt} = -\nabla P
```

where:
```math
\frac{d}{dt} = \frac{\partial}{\partial t} + \bm{v} \cdot \nabla
```

is the material derivative, and the pressure is:
```math
P = \rho^2 \frac{\partial u}{\partial \rho}\bigg|_s
```

**Energy Equation**:
```math
\frac{de}{dt} = -P \nabla \cdot \bm{v}
```

where $e = \rho u$ is the internal energy density.

#### Derivation of Momentum Equation

Starting from the Lagrangian $\mathcal{L}$, the momentum equation follows from:

```math
\frac{d}{dt}\left(\frac{\delta \mathcal{L}}{\delta \bm{v}}\right) = \frac{\delta \mathcal{L}}{\delta \bm{x}}
```

Computing the variational derivatives:
```math
\frac{\delta \mathcal{L}}{\delta \bm{v}} = \rho \bm{v}
```

```math
\frac{\delta \mathcal{L}}{\delta \bm{x}} = -\nabla \left( \rho \frac{\partial u}{\partial \rho} \right) = -\nabla P
```

This yields:
```math
\frac{d(\rho \bm{v})}{dt} = -\nabla P
```

Using the continuity equation to expand the material derivative:
```math
\rho \frac{d\bm{v}}{dt} + \bm{v} \frac{d\rho}{dt} = -\nabla P
```

Since $\frac{d\rho}{dt} = -\rho \nabla \cdot \bm{v}$, we recover:
```math
\rho \frac{d\bm{v}}{dt} = -\nabla P
```

#### Canonical Momentum (Non-Relativistic)

The canonical momentum density is:
```math
\bm{\pi} = \frac{\delta \mathcal{L}}{\delta \bm{v}} = \rho \bm{v}
```

The momentum per unit mass is simply $\bm{v}$.

### Relativistic Lagrangian Formulation

#### Four-Dimensional Action Principle

In special relativity, the action is constructed to be Lorentz invariant:

```math
S = \int d^4x \, \mathcal{L}
```

where $d^4x = dt \, d^3x$ and $\mathcal{L}$ is the Lagrangian density (scalar under Lorentz transformations).

#### Lagrangian Density for Perfect Fluid

For a relativistic perfect fluid, the Lagrangian density is:

```math
\mathcal{L} = -\rho c^2 \sqrt{1 - \frac{\bm{v}^2}{c^2}} - \rho u(\rho, s) \sqrt{1 - \frac{\bm{v}^2}{c^2}}
```

This can be rewritten using the Lorentz factor $\gamma = (1 - \bm{v}^2/c^2)^{-1/2}$:

```math
\mathcal{L} = -\frac{\rho}{\gamma} \left( c^2 + u \right)
```

where $\rho/\gamma = n$ is the rest-frame (proper) baryon density.

#### Alternative Form Using Baryon Number Density

Using the rest-frame baryon density $n$ and enthalpy per baryon $H$:

```math
\mathcal{L} = -n H c^2 / \gamma = -n c^2 (1 + u/c^2 + P/(nc^2)) / \gamma
```

where:
```math
H = 1 + \frac{u}{c^2} + \frac{P}{nc^2}
```

is the enthalpy per baryon **[SRGSPH Eq. 112]**.

#### Energy-Momentum Tensor

The energy-momentum tensor for a perfect fluid is:

```math
T^{\mu\nu} = (e + P) u^\mu u^\nu + P \eta^{\mu\nu}
```

where:
- $e$ is the energy density (including rest mass)
- $P$ is the pressure
- $u^\mu = \gamma(c, \bm{v})$ is the four-velocity
- $\eta^{\mu\nu} = \text{diag}(-1, 1, 1, 1)$ is the Minkowski metric

#### Conservation Laws from Energy-Momentum Tensor

The equations of motion follow from energy-momentum conservation:

```math
\partial_\mu T^{\mu\nu} = 0
```

**Spatial components ($\nu = i$)**: Momentum conservation
```math
\partial_t T^{0i} + \partial_j T^{ji} = 0
```

**Time component ($\nu = 0$)**: Energy conservation
```math
\partial_t T^{00} + \partial_i T^{i0} = 0
```

#### Canonical Momentum (Relativistic)

The canonical momentum density is derived from the Lagrangian:

```math
\bm{\pi} = \frac{\delta \mathcal{L}}{\delta \bm{v}} = \gamma^2 n H \bm{v}
```

The **canonical momentum per baryon** is **[SRGSPH Eq. 104]**:

```math
\bm{S} = \frac{\bm{\pi}}{N} = \gamma H \bm{v}
```

where $N = \gamma n$ is the lab-frame baryon density.

#### Derivation of Relativistic Euler Equation

Starting from the spatial components of $\partial_\mu T^{\mu i} = 0$:

```math
\partial_t [(\rho h \gamma^2) v^i] + \partial_j [(\rho h \gamma^2) v^i v^j + P \delta^{ij}] = 0
```

where $\rho = m_b n$ is the mass density and $h = H$ when $c=1$.

Expanding and using the continuity equation, this becomes:

```math
\rho h \gamma^2 \frac{dv^i}{dt} + v^i \frac{d(\rho h \gamma^2)}{dt} = -\partial^i P
```

#### SRGSPH Form of Momentum Equation

Using the lab-frame baryon density $N = \gamma n$ and canonical momentum per baryon $\bm{S} = \gamma H \bm{v}$, the momentum equation becomes **[SRGSPH Eq. 92]**:

```math
\frac{d\bm{S}}{dt} = -\frac{1}{N} \nabla P
```

This is the **equation of motion** used in SRGSPH.

#### Derivation of Canonical Energy Equation

The canonical energy per baryon in the lab frame is defined as **[SRGSPH Eq. 107]**:

```math
e = \gamma H - \frac{P}{Nc^2}
```

The time evolution follows from energy conservation **[SRGSPH Eq. 94]**:

```math
\frac{de}{dt} = -\frac{1}{N} \nabla \cdot (P\bm{v})
```

#### Physical Interpretation

**Non-relativistic limit**: As $\bm{v}/c \to 0$:
- $\gamma \to 1$
- $H \to 1 + u/c^2 + P/(nc^2) \approx 1$ (for non-relativistic temperatures)
- $\bm{S} = \gamma H \bm{v} \approx \bm{v}$ (canonical momentum → velocity)
- $e = \gamma H - P/(Nc^2) \approx 1 + u/c^2$ (rest mass + thermal energy)

**Relativistic regime**: When $\bm{v} \sim c$ or $u \sim nc^2$:
- Lorentz factor $\gamma \gg 1$ couples all velocity components
- Enthalpy $H \gg 1$ amplifies pressure effects
- Canonical momentum $\bm{S} = \gamma H \bm{v}$ differs significantly from $\bm{v}$
- Energy $e$ includes kinetic and internal contributions weighted by $\gamma H$

#### Comparison of Lagrangian Approaches

| Aspect | Non-Relativistic | Relativistic (SRGSPH) |
|--------|------------------|----------------------|
| Action | $S = \int dt \int d^3x \, \mathcal{L}$ | $S = \int d^4x \, \mathcal{L}$ |
| Lagrangian density | $\mathcal{L} = \rho[\frac{1}{2}\bm{v}^2 - u]$ | $\mathcal{L} = -nH/\gamma$ (with $c=1$) |
| Canonical momentum | $\bm{\pi} = \rho \bm{v}$ | $\bm{\pi} = \gamma^2 n H \bm{v}$ |
| Momentum per particle | $\bm{v}$ | $\bm{S} = \gamma H \bm{v}$ |
| EoM form | $\frac{d\bm{v}}{dt} = -\frac{1}{\rho}\nabla P$ | $\frac{d\bm{S}}{dt} = -\frac{1}{N}\nabla P$ |
| Energy equation | $\frac{de}{dt} = -P\nabla \cdot \bm{v}$ | $\frac{de}{dt} = -\frac{1}{N}\nabla \cdot (P\bm{v})$ |
| Symmetry | Galilean invariance | Lorentz invariance |

#### Connection to Simulation Implementation

The simulation implements these Lagrangian-derived equations in discretized SPH form:

**Non-relativistic SPH** (classical GSPH, DISPH):
- Evolves $\bm{v}$ directly
- Pressure gradient computed from kernel-weighted sums
- Energy equation for thermal evolution

**Relativistic SPH** (SRGSPH):
- Evolves canonical variables $\bm{S}$ and $e$
- Riemann solver computes $P_{ij}^*$ and $\bm{v}_{ij}^*$ at interfaces
- Recovers primitive variables $(\bm{v}, n, P, u)$ from $(\bm{S}, e)$ using quartic solver
- All formulations preserve Lorentz invariance of the underlying Lagrangian

---

## Theory of Relativistic Hydrodynamics

### Why Relativistic Hydrodynamics is Different

**[Based on Pons §1 and SRGSPH §1]**

In classical hydrodynamics, tangential velocities remain constant across waves. In relativistic hydrodynamics, this fundamentally changes:

1. **Lorentz factor coupling**: All velocity components couple through:
   ```math
   \gamma = \frac{1}{\sqrt{1 - \bm{v}^2/c^2}} = \frac{1}{\sqrt{1 - v^2/c^2}}
   ```
   where $v^2 = (v^x)^2 + (v^y)^2 + (v^z)^2$

2. **Enthalpy-momentum coupling**: Momentum per baryon is:
   ```math
   \bm{S} = \gamma H \bm{v}
   ```
   Even for slow flows ($\gamma \approx 1$), if $H > 1$ (thermodynamically relativistic), tangential velocities matter.

3. **No universal rest frame**: For Riemann problems with arbitrary tangential velocities, no single frame exists where all tangential components vanish.

### Conservation Laws

The fundamental conservation laws:
- **Rest-mass (baryon number) conservation**: $N = \gamma n$ conserved per fluid element
- **Energy-momentum conservation**: Encoded in energy-momentum tensor

For a perfect fluid:
```math
T^{\mu\nu} = (energy\_density + P) u^\mu u^\nu + P \eta^{\mu\nu}
```

where the energy-momentum includes rest mass energy.

---

## Basic Equations (SRGSPH Formulation)

### Lagrangian Derivative

**[SRGSPH Eq. 100]**:
```math
\frac{d}{dt} \equiv \partial_t + \bm{v} \cdot \nabla
```

### Fundamental Equations

**Equation of Continuity** **[SRGSPH Eq. 90]**:
```math
\frac{dN}{dt} = -N \nabla \cdot \bm{v}
```

**Equation of Motion** **[SRGSPH Eq. 92]**:
```math
\frac{d\bm{S}}{dt} = -\frac{1}{N} \nabla P
```

**Equation of Energy** **[SRGSPH Eq. 94]**:
```math
\frac{de}{dt} = -\frac{1}{N} \nabla \cdot (P\bm{v})
```

where $N$, $\bm{S}$, and $e$ are baryon number density, relativistic canonical momentum per baryon, and canonical energy per baryon in the lab frame **[SRGSPH, discussion after Eq. 96]**.

### Relativistic Canonical Quantities

**Canonical momentum per baryon** **[SRGSPH Eq. 104]**:
```math
\bm{S} = \gamma H \bm{v}
```

**Canonical energy per baryon** **[SRGSPH Eq. 107]**:
```math
e = \gamma H - \frac{P}{Nc^2}
```

**Lorentz factor** **[SRGSPH Eq. 111]**:
```math
\gamma = \frac{1}{\sqrt{1 - \bm{v}^2/c^2}}
```

**Enthalpy per baryon** **[SRGSPH Eq. 112]**:
```math
H = 1 + \frac{u}{c^2} + \frac{P}{nc^2}
```

where $u$ is thermal energy per baryon and $n$ is baryon number density in rest frame.

**Lab-rest frame relation** **[SRGSPH Eq. 114]**:
```math
N = \gamma n
```

### Equation of State

For an ideal gas **[SRGSPH Eq. 119]**:
```math
P = (\gamma_c - 1) n u
```

where $\gamma_c$ is the adiabatic index (ratio of specific heats).

### Sound Speed

For an ideal gas **[SRGSPH, after Eq. 394]**:
```math
c_s^2 = \frac{(\gamma_c - 1)(H - 1) c^2}{H}
```

**With $c=1$**:
```math
c_s = \sqrt{\frac{(\gamma_c - 1)(H - 1)}{H}}
```

---

## Conservative Formulation

For grid-based methods (context for understanding Riemann problems), the conserved variables are:

**Rest-mass density** (with $c=1$):
```math
D = \rho \gamma = m_b n \gamma = m_b N
```

**Momentum density** (with $c=1$, using specific enthalpy $h = H$ when $c=1$):
```math
S^i = \rho h \gamma^2 v^i
```

**Energy density** (with $c=1$):
```math
\tau = \rho h \gamma^2 - P - D
```

These satisfy conservation form:
```math
\mathbf{U}_{,t} + \mathbf{F}^{(i)}_{,i} = 0
```

---

## SRGSPH Formulation with Fixed Smoothing Length

### Kernel Function

**[SRGSPH Eq. 138]**: Gaussian kernel in $d$ dimensions:
```math
W(\bm{x}, h) = \left[ \frac{1}{h\sqrt{\pi}} \right]^d \exp\left[ -\frac{\bm{x}^2}{h^2} \right]
```

### Number Density Field

**[SRGSPH Eq. 134]**:
```math
N(\bm{x}) = \nu \sum_j W(\bm{x} - \bm{x}_j, h)
```

where $\nu$ is the (constant) baryon number per SPH particle.

### Convolution

**[SRGSPH Eqs. 146-148]**:
```math
\langle N_i \rangle = N(\bm{x}_i)
```
```math
\langle f_i \rangle = \int f(\bm{x}) W(\bm{x} - \bm{x}_i, h) d\bm{x}
```

**Accuracy** **[SRGSPH Eq. 154]**: Error is $\mathcal{O}(h^2)$:
```math
\langle f_i \rangle = f(\bm{x}_i) + \frac{h^2}{4} f''(\bm{x}_i) + \mathcal{O}(h^4)
```

### Equation of Motion (Fixed h)

**[SRGSPH Eqs. 162-166]**, starting from Lagrangian form:
```math
\langle \dot{\bm{S}}_i \rangle = -\sum_j \int \frac{\nu P(\bm{x})}{N^2(\bm{x})} [\nabla_i - \nabla_j] W(\bm{x} - \bm{x}_i, h) W(\bm{x} - \bm{x}_j, h) d\bm{x}
```

### Equation of Energy (Fixed h)

**[SRGSPH Eqs. 176-181]**:
```math
\langle \dot{e}_i \rangle = -\sum_j \int \frac{\nu P(\bm{x}) \bm{v}(\bm{x})}{N^2(\bm{x})} \cdot [\nabla_i - \nabla_j] W(\bm{x} - \bm{x}_i, h) W(\bm{x} - \bm{x}_j, h) d\bm{x}
```

### Evaluation of Convolution Integrals

**[SRGSPH Eqs. 188-191]**, using product of Gaussians:
```math
\int \frac{f(\bm{x})}{N^2(\bm{x})} [\nabla_i - \nabla_j] W(\bm{x} - \bm{x}_i, h) W(\bm{x} - \bm{x}_j, h) d\bm{x}
```
```math
= f_{ij} V_{ij}^2 [\nabla_i W(\bm{x}_i - \bm{x}_j, \sqrt{2}h) - \nabla_j W(\bm{x}_i - \bm{x}_j, \sqrt{2}h)]
```

where:

**Effective volume factor** **[SRGSPH Eq. 194]**:
```math
V_{ij}^2(h) = \int \frac{1}{N^2(\bm{x})} \left( \frac{\sqrt{2}}{h\sqrt{\pi}} \right)^d \exp\left[ -\frac{2}{h^2} \left( \bm{x} - \frac{\bm{x}_i + \bm{x}_j}{2} \right)^2 \right] d\bm{x}
```

**Weighted average** **[SRGSPH Eq. 196]**:
```math
f_{ij}(h) = \frac{1}{V_{ij}^2} \int \frac{f(\bm{x})}{N^2(\bm{x})} \left( \frac{\sqrt{2}}{h\sqrt{\pi}} \right)^d \exp\left[ -\frac{2}{h^2} \left( \bm{x} - \frac{\bm{x}_i + \bm{x}_j}{2} \right)^2 \right] d\bm{x}
```

### Approximations

**[SRGSPH, discussion after Eq. 200]**:
- $V_{ij}^2$ approximated by interpolation: $V_{ij}^2 \approx V_{ij,\text{interp}}^2$
- $f_{ij}$ replaced by Riemann solution: $f_{ij} \approx f_{ij}^*$

### Final SRGSPH Equations (Fixed h)

**Equation of Motion** **[SRGSPH Eq. 209]**:
```math
\langle \dot{\bm{S}}_i \rangle = -\sum_j \nu P_{ij}^* V_{ij,\text{interp}}^2 [\nabla_i W(\bm{x}_i - \bm{x}_j, \sqrt{2}h) - \nabla_j W(\bm{x}_i - \bm{x}_j, \sqrt{2}h)]
```

**Equation of Energy** **[SRGSPH Eq. 212]**:
```math
\langle \dot{e}_i \rangle = -\sum_j \nu P_{ij}^* \bm{v}_{ij}^* \cdot V_{ij,\text{interp}}^2 [\nabla_i W(\bm{x}_i - \bm{x}_j, \sqrt{2}h) - \nabla_j W(\bm{x}_i - \bm{x}_j, \sqrt{2}h)]
```

**Conservation**: Anti-symmetry in $(i,j)$ ensures total momentum and energy conservation **[SRGSPH Eq. 214]**.

---

## SRGSPH Formulation with Variable Smoothing Length

### Particle Volume Definition

**[SRGSPH Eq. 221]**:
```math
V_{\rm p}(\bm{x}) = \left[ \sum_j W(\bm{x} - \bm{x}_j, h(\bm{x})) \right]^{-1}
```

### Smoothing Length Definition

**[SRGSPH Eq. 231]**:
```math
h(\bm{x}) = \eta [V_{\rm p}^*(\bm{x})]^{1/d}
```

where **[SRGSPH Eq. 233]**:
```math
V_{\rm p}^*(\bm{x}) = \left[ \sum_j W(\bm{x} - \bm{x}_j, C_{\rm smooth} h(\bm{x})) \right]^{-1}
```

**Typical values** **[SRGSPH, after Eq. 236]**: $\eta = 1.0$, $C_{\rm smooth} = 2.0$

**Iterative solution** **[SRGSPH, discussion after Eq. 238]**: Iterate Eqs. (231) and (233) until convergence.

### Number Density (Volume-Based)

**[SRGSPH Eq. 243]**:
```math
N(\bm{x}) = \nu(\bm{x}) V_{\rm p}^{-1}(\bm{x})
```

where $\nu(\bm{x})$ varies spatially.

### Final SRGSPH Equations (Variable h)

**Equation of Motion** **[SRGSPH Eq. 371]**:
```math
\langle \nu_i \dot{\bm{S}}_i \rangle = -\sum_j P_{ij}^* V_{ij,\text{interp}}^2 [\nabla_i W(\bm{x}_i - \bm{x}_j, \sqrt{2}h_i) - \nabla_j W(\bm{x}_i - \bm{x}_j, \sqrt{2}h_j)]
```

**Equation of Energy** **[SRGSPH Eq. 373]**:
```math
\langle \nu_i \dot{e}_i \rangle = -\sum_j P_{ij}^* \bm{v}_{ij}^* \cdot V_{ij,\text{interp}}^2 [\nabla_i W(\bm{x}_i - \bm{x}_j, \sqrt{2}h_i) - \nabla_j W(\bm{x}_i - \bm{x}_j, \sqrt{2}h_j)]
```

where averaging ensures symmetry: $V_{ij}^2 = (V_{ij}^2(h_i) + V_{ij}^2(h_j))/2$ **[SRGSPH Eq. 365]**.

---

## Riemann Problem Theory

### Problem Setup

Initial left and right states at $t=0$, $x=0$:
- Left: $(P_L, n_L, v_L^x, v_L^y, v_L^z, u_L)$
- Right: $(P_R, n_R, v_R^x, v_R^y, v_R^z, u_R)$

Compute derived quantities using SRGSPH formulations:
- $\gamma_L, \gamma_R$ from velocities [Eq. 111]
- $H_L, H_R$ from EOS [Eqs. 112, 119]

### Wave Structure

**[Pons §5]**: Three waves emerge:
```math
L \; \mathcal{W}_\leftarrow \; L_* \; \mathcal{C} \; R_* \; \mathcal{W}_\rightarrow \; R
```

- $\mathcal{W}$: shock or rarefaction
- $\mathcal{C}$: contact discontinuity

**Contact discontinuity**:
- Continuous: $P$, $v^x$
- Discontinuous: $n$, $v^{y,z}$, $u$

### Self-Similarity

All waves satisfy $\xi = x/t = \text{constant}$, reducing PDEs to ODEs.

### Riemann Invariants (Key Difference from Newtonian)

**[Pons Eqs. 3.8-3.9, adapted to SRGSPH]**:

Across rarefactions and shocks, with $c=1$:
```math
H \gamma v^y = \text{constant}
```
```math
H \gamma v^z = \text{constant}
```

**Physical meaning**: The "relativistic tangential momentum per baryon" is conserved, but actual tangential velocity changes due to enthalpy and Lorentz factor variations.

---

## Riemann Problem: Rarefaction Waves

### Characteristic Speeds

**[Pons Eq. 3.6]**, with $c=1$:
```math
\xi = \frac{v^x(1 - c_s^2) \pm c_s \sqrt{(1 - v^2)[1 - v^2 c_s^2 - (v^x)^2(1 - c_s^2)]}}{1 - v^2 c_s^2}
```

### Reduced System

**[Pons Eqs. 3.7-3.9, with SRGSPH variables, $c=1$]**:

**ODE:**
```math
\rho_m h \gamma^2 (v^x - \xi) dv^x + (1 - \xi v^x) dP = 0
```

where $\rho_m = m_b n$ is mass density, and $h = H$ when $c=1$.

**Riemann invariants:**
```math
H \gamma v^y = \text{constant}
```
```math
H \gamma v^z = \text{constant}
```

### ODE in Standard Form

**[Pons Eq. 3.10, with $c=1$]**:
```math
\frac{dv^x}{dP} = \pm \frac{1}{\rho_m h \gamma^2 c_s} \frac{1}{\sqrt{1 + g(\xi_\pm, v^x, v^t)}}
```

where:
```math
g(\xi_\pm, v^x, v^t) = \frac{(v^t)^2 (\xi_\pm^2 - 1)}{(1 - \xi_\pm v^x)^2}
```

### Post-Rarefaction Tangential Velocity

**[Pons Eq. 3.11, with SRGSPH $H$, $c=1$]**:
```math
v_b^t = H_a \gamma_a v_a^t \left\{ \frac{1 - (v_b^x)^2}{H_b^2 + (H_a \gamma_a v_a^t)^2} \right\}^{1/2}
```

---

## Riemann Problem: Shock Waves

### Rankine-Hugoniot Jump Conditions

**[Pons Eqs. 4.1-4.2]**: Mass and energy-momentum conservation across shocks.

### Invariant Mass Flux

**[Pons Eq. 4.7, with $c=1$]**:
```math
j = \gamma_s D_a (V_s - v_a^x) = \gamma_s D_b (V_s - v_b^x)
```

where $D = m_b N = m_b \gamma n$ and:
```math
\gamma_s = \frac{1}{\sqrt{1 - V_s^2}}
```

### Tangential Velocity Invariant

**[Pons Eqs. 4.10-4.11, with SRGSPH $H$, $c=1$]**:
```math
H \gamma v^y = \text{constant across shock}
```
```math
H \gamma v^z = \text{constant across shock}
```

### Taub Adiabat

**[Pons Eq. 4.15, with SRGSPH $H$]**:
```math
[H^2] = \left( \frac{H_b}{n_b} + \frac{H_a}{n_a} \right) [P]
```

For ideal gas, this becomes a quadratic for $H_b$ **[Pons Eq. 4.16, adapted]**.

### Mass Flux from Pressure

**[Pons Eq. 4.17]**:
```math
j^2 = \frac{-[P]}{[H/(m_b n)]}
```

---

## Complete Riemann Solver Algorithm

### Solution Procedure

**[Pons §5, SRGSPH §2.7.2]**:

1. **Setup**: Convert primitive variables to SRGSPH canonical form
   - Compute $\gamma_L, \gamma_R, H_L, H_R$ using SRGSPH Eqs. 111, 112, 119

2. **Define wave functions**: $v_{L*}^x(P)$, $v_{R*}^x(P)$ using rarefaction/shock relations

3. **Solve for $P_*$**: Root of $v_{L*}^x(P) - v_{R*}^x(P) = 0$

4. **Compute intermediate states**: $L_*$, $R_*$ using appropriate wave relations

5. **Extract interface state** **[SRGSPH Eq. 423]**:
   - If $v_*^x > 0$: use $L_*$
   - If $v_*^x < 0$: use $R_*$
   - Return $(P_{ij}^*, \bm{v}_{ij}^*)$

---

## Primitive Variable Recovery

**[SRGSPH §2.8]**:

After time integration, recover primitives from $(\langle \bm{S}_i \rangle, \langle e_i \rangle)$.

### Quartic Equation for Lorentz Factor

**[SRGSPH, after Eq. 409]**: With $X = \gamma_c/(\gamma_c - 1)$:
```math
(\langle \gamma_i \rangle^2 - 1)(X \langle e_i \rangle \langle \gamma_i \rangle - 1)^2 - \langle \bm{S}_i \rangle^2 (X \langle \gamma_i \rangle^2 - 1)^2 = 0
```

Solve numerically for $\langle \gamma_i \rangle$.

### Velocity Recovery

**[SRGSPH Eq. 416]**:
```math
\langle \bm{v}_i \rangle = \frac{X \langle \gamma_i \rangle^2 - 1}{\langle \gamma_i \rangle (X \langle e_i \rangle \langle \gamma_i \rangle - 1)} \langle \bm{S}_i \rangle
```

### Other Primitive Variables

From $\langle \gamma_i \rangle$ and known $\langle N_i \rangle$:
```math
\langle n_i \rangle = \frac{\langle N_i \rangle}{\langle \gamma_i \rangle}
```

```math
\langle H_i \rangle = \frac{X \langle e_i \rangle \langle \gamma_i \rangle - 1}{\langle \gamma_i \rangle^2}
```

```math
\langle P_i \rangle = \langle n_i \rangle (\langle H_i \rangle - 1) \frac{\gamma_c - 1}{\gamma_c} c^2
```

With $c=1$: $\langle P_i \rangle = \langle n_i \rangle (\langle H_i \rangle - 1) (\gamma_c - 1)/\gamma_c$

---

## Time Integration

### Euler Method

**[SRGSPH Eqs. 425-428]**:
```math
\langle \nu_i \bm{S}_i \rangle^{n+1} = \langle \nu_i \bm{S}_i \rangle^n + \langle \nu_i \dot{\bm{S}}_i \rangle \Delta t
```
```math
\langle \nu_i e_i \rangle^{n+1} = \langle \nu_i e_i \rangle^n + \langle \nu_i \dot{e}_i \rangle \Delta t
```
```math
\bm{x}_i^{n+1} = \bm{x}_i^n + \langle \bm{v}_i \rangle \Delta t
```

### CFL Condition

**[SRGSPH Eqs. 432-433]**:
```math
\Delta t = C_{\text{CFL}} \min_i \left[ \frac{h_i}{c_{s,i}} \right]
```

Typical: $C_{\text{CFL}} = 0.3$

### MUSCL Reconstruction

**[SRGSPH §2.7.3]**: For second-order accuracy, reconstruct states at interfaces using gradients with shock detection and limiting.

---

## Numerical Implementation (c=1)

### Setting c=1

For numerical calculations **[SRGSPH, numerical sections]**, set $c = 1$:

- $\gamma = 1/\sqrt{1 - v^2}$
- $H = 1 + u + P/n$
- $\bm{S} = \gamma H \bm{v}$
- $e = \gamma H - P/N$
- $c_s = \sqrt{(\gamma_c - 1)(H - 1)/H}$

### Correspondence with Pons Notation

When $c=1$ and using mass density $\rho = m_b n$:
- SRGSPH $H$ ↔ Pons $h$ (numerically equal)
- SRGSPH $n$ ↔ Pons $\rho/m_b$ (or $\rho$ with $m_b=1$)
- SRGSPH $u$ ↔ Pons $\epsilon \rho/n$ (or $\epsilon$ with $m_b=1$)

This allows direct use of Pons's Riemann solver formulas with SRGSPH variables.

---

## References

1. **[SRGSPH]**: Kitajima, K., Inutsuka, S., & Seno, I. (2025). Special Relativistic Smoothed Particle Hydrodynamics Based on Riemann Solver. Journal of Computational Physics (accepted).

2. **[Pons]**: Pons, J.A., Martí, J.M., & Müller, E. (2000). The exact solution of the Riemann problem with non-zero tangential velocities in relativistic hydrodynamics. Journal of Fluid Mechanics, 422, 125-139.

3. Inutsuka, S. (2002). Reformulation of Smoothed Particle Hydrodynamics with Riemann Solver. Journal of Computational Physics, 179, 238-267.

---

*This document uses SRGSPH (Kitajima et al.) formulations throughout, with explicit $c$ in theory and $c=1$ for numerics. Riemann solver details from Pons et al. are adapted to SRGSPH notation.*
