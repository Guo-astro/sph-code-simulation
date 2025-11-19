# Special Relativistic Godunov SPH: Complete Mathematical Formulation

This document provides a comprehensive compilation of all formulas for Special Relativistic Godunov Smoothed Particle Hydrodynamics (SRGSPH), including detailed Riemann solver formulations. **All symbols follow the Kitajima et al. (SRGSPH) convention** as the single source of truth, with Pons et al. formulas translated to use consistent notation.

**Reference Papers:**
- **[SRGSPH]**: Kitajima, K., Inutsuka, S., & Seno, I. (2025). Special Relativistic Smoothed Particle Hydrodynamics Based on Riemann Solver.
- **[Pons]**: Pons, J.A., Martí, J.M., & Müller, E. (2000). The exact solution of the Riemann problem with non-zero tangential velocities in relativistic hydrodynamics.
  - **Note**: Pons et al. use $W$ for Lorentz factor (here: $\gamma$), $p$ for pressure (here: $P$), and $\gamma$ for adiabatic index (here: $\gamma_c$)

---

## Table of Contents

1. [Notation and Conventions](#notation-and-conventions)
2. [Theory of Relativistic Hydrodynamics](#theory-of-relativistic-hydrodynamics)
3. [Basic Equations of Special Relativistic Hydrodynamics](#basic-equations-of-special-relativistic-hydrodynamics)
4. [Conservative Formulation](#conservative-formulation)
5. [SRGSPH Formulation with Fixed Smoothing Length](#srgsph-formulation-with-fixed-smoothing-length)
6. [SRGSPH Formulation with Variable Smoothing Length](#srgsph-formulation-with-variable-smoothing-length)
7. [Theory of Relativistic Riemann Problems](#theory-of-relativistic-riemann-problems)
8. [Riemann Problem: Rarefaction Waves](#riemann-problem-rarefaction-waves)
9. [Riemann Problem: Shock Waves](#riemann-problem-shock-waves)
10. [Complete Riemann Solver Algorithm](#complete-riemann-solver-algorithm)
11. [Primitive Variable Recovery](#primitive-variable-recovery)
12. [Time Integration](#time-integration)

---
# Special Relativistic Godunov SPH: Complete Mathematical Formulation

This document provides a comprehensive compilation of all formulas for Special Relativistic Godunov Smoothed Particle Hydrodynamics (SRGSPH), including detailed Riemann solver formulations. **All symbols follow the Kitajima et al. (SRGSPH) convention** as the single source of truth, with Pons et al. formulas translated to use consistent notation.

**Reference Papers:**
- **[SRGSPH]**: Kitajima, K., Inutsuka, S., & Seno, I. (2025). Special Relativistic Smoothed Particle Hydrodynamics Based on Riemann Solver.
- **[Pons]**: Pons, J.A., Martí, J.M., & Müller, E. (2000). The exact solution of the Riemann problem with non-zero tangential velocities in relativistic hydrodynamics.
  - **Note**: Pons et al. use $W$ for Lorentz factor (here: $\gamma$), $p$ for pressure (here: $P$), and $\gamma$ for adiabatic index (here: $\gamma_c$)

---

## Table of Contents

1. [Notation and Conventions](#notation-and-conventions)
2. [Theory of Relativistic Hydrodynamics](#theory-of-relativistic-hydrodynamics)
3. [Basic Equations of Special Relativistic Hydrodynamics](#basic-equations-of-special-relativistic-hydrodynamics)
4. [Conservative Formulation](#conservative-formulation)
5. [SRGSPH Formulation with Fixed Smoothing Length](#srgsph-formulation-with-fixed-smoothing-length)
6. [SRGSPH Formulation with Variable Smoothing Length](#srgsph-formulation-with-variable-smoothing-length)
7. [Theory of Relativistic Riemann Problems](#theory-of-relativistic-riemann-problems)
8. [Riemann Problem: Rarefaction Waves](#riemann-problem-rarefaction-waves)
9. [Riemann Problem: Shock Waves](#riemann-problem-shock-waves)
10. [Complete Riemann Solver Algorithm](#complete-riemann-solver-algorithm)
11. [Primitive Variable Recovery](#primitive-variable-recovery)
12. [Time Integration](#time-integration)

---

## Notation and Conventions

### Units
- Speed of light: $c = 1$ (for numerical calculations, set explicitly in coordinate system choice)
- Minkowski metric: $\eta^{\mu\nu} = \text{diag}(-1, 1, 1, 1)$ **[Pons Eq. after 2.3]**

### Symbols (SRGSPH Convention)

**Spacetime coordinates:**
- $t$: time
- $\mathbf{x} = (x, y, z)$: spatial position vector
- $x^\mu = (t, x, y, z)$: spacetime coordinates ($\mu = 0, 1, 2, 3$) **[Pons §2]**

**Primitive variables (physical quantities):**
- $n$: baryon number density in rest frame **[SRGSPH]**
- $N$: baryon number density in lab frame **[SRGSPH Eq. 91]**
- $\rho$: proper rest-mass density **[Pons Eq. 2.1]**
- $P$: pressure **[SRGSPH Eq. 92]** (Pons uses lowercase $p$)
- $\mathbf{v} = (v^x, v^y, v^z)$: three-velocity **[Pons Eq. 2.5]**
- $v^n$: normal velocity component (perpendicular to discontinuity)
- $v^t = \sqrt{(v^y)^2 + (v^z)^2}$: tangential velocity magnitude **[Pons §3]**
- $u$: thermal energy per baryon **[SRGSPH]**
- $\epsilon$: specific internal energy **[Pons Eq. 2.2]**
- $\gamma_c$: ratio of specific heats (adiabatic index) **[SRGSPH Eq. 119]** (Pons uses $\gamma$)
- $c_s$: sound speed **[Pons Eq. 2.9]**

**Relativistic quantities:**
- $\gamma$: Lorentz factor **[SRGSPH Eq. 111]** (Pons uses $W$)
- $H$: enthalpy per baryon **[SRGSPH Eq. 112]**
- $h$: specific enthalpy **[Pons Eq. 2.2]**
- $\mathbf{S}$: relativistic canonical momentum per baryon **[SRGSPH Eq. 104]**
- $e$: canonical energy per baryon in lab frame **[SRGSPH Eq. 107]**
- $u^\mu$: four-velocity **[Pons Eq. 2.1]**

**Conserved variables:**
- $D$: rest-mass density (conserved) **[Pons Eq. 2.10]**
- $S^i$: momentum density **[Pons Eq. 2.10]**
- $\tau$: energy density **[Pons Eq. 2.10]**

**SPH-specific quantities:**
- $\nu$: baryon number in an SPH particle **[SRGSPH]**
- $W(\mathbf{x}, h)$: kernel function **[SRGSPH Eq. 138]** (note: different from Lorentz factor)
- $h$: smoothing length (context distinguishes from enthalpy)
- $V_p$: particle volume **[SRGSPH Eq. 221]**
- $\langle f_i \rangle$: convolved quantity for particle $i$ **[SRGSPH Eq. 148]**
- $d$: number of spatial dimensions
- $\eta$: smoothing length parameter **[SRGSPH Eq. 231]**
- $C_{\text{smooth}}$: smoothing coefficient **[SRGSPH Eq. 233]**

**Riemann problem notation:**
- $L, R$: left and right initial states **[Pons §1]**
- $L_*, R_*$: left and right intermediate states **[Pons §5]**
- $P_*$: pressure in intermediate state **[Pons Eq. 1.1]**
- $v^x_*$: normal velocity in intermediate state **[Pons Eq. 1.2]**
- $a, b$: states ahead and behind a wave
- $\xi$: self-similarity variable ($x/t$) **[Pons Eq. 3.1]**
- $j$: invariant mass flux across shock **[Pons Eq. 4.7]**
- $V_s$: shock velocity **[Pons Eq. 4.4]**
- $\gamma_s$: shock Lorentz factor **[Pons Eq. 4.5]** (Pons uses $W_s$)

---

## Theory of Relativistic Hydrodynamics

### Why Relativistic Hydrodynamics is Different

**[Based on Pons §1 and SRGSPH §1]**

In classical (Newtonian) hydrodynamics, the decay of an initial discontinuity (Riemann problem) has a key simplifying feature: **tangential velocities do not affect the solution**. The tangential velocity component remains constant across both shock waves and rarefaction waves, allowing the problem to be solved in a frame where tangential velocities vanish.

**In relativistic hydrodynamics, this is fundamentally different:**

1. **Lorentz factor coupling**: All velocity components (normal and tangential) are coupled through the Lorentz factor:
   ```math
   \gamma = \frac{1}{\sqrt{1 - v^2}} = \frac{1}{\sqrt{1 - (v^x)^2 - (v^y)^2 - (v^z)^2}}
   ```

2. **No universal rest frame**: For a Riemann problem with arbitrary tangential velocities on both sides, there is **no single reference frame** in which the tangential velocity vanishes everywhere. Unlike isolated shocks or rarefactions (where one can always choose a comoving frame), the Riemann problem requires the full treatment.

3. **Enthalpy coupling**: The specific enthalpy $h$ appears in the momentum density:
   ```math
   S^i = \rho h \gamma^2 v^i
   ```
   This means even for slow flows ($\gamma \approx 1$), if the fluid is thermodynamically relativistic ($h > 1$, i.e., internal energy or pressure comparable to rest mass), tangential velocities significantly affect the dynamics.

### Physical Interpretation

**Conservation Laws:**

The fundamental conservation laws in relativistic hydrodynamics are:
- **Rest-mass conservation**: Baryon number is conserved
- **Energy-momentum conservation**: The energy-momentum tensor is conserved

These are encoded in:
```math
J^\mu_{,\mu} = 0 \quad \text{(rest-mass)}
```
```math
T^{\mu\nu}_{,\mu} = 0 \quad \text{(energy-momentum)}
```

Unlike Newtonian mechanics, energy and momentum are unified in the energy-momentum tensor $T^{\mu\nu}$, and their conservation is geometrically linked.

**Why Tangential Velocities Matter:**

1. **In rarefactions**: Fluid elements experience continuous changes. The Riemann invariant:
   ```math
   h \gamma v^t = \text{constant along fluid lines}
   ```
   means that as enthalpy $h$ changes (due to pressure/density changes), the tangential velocity must adjust to maintain the invariant.

2. **Across shocks**: The Rankine-Hugoniot jump conditions preserve:
   ```math
   h \gamma v^{y,z} = \text{constant across shock}
   ```
   The jump in pressure and density causes jumps in $h$ and $\gamma$, requiring the tangential velocity to jump as well.

3. **Limiting behavior**: When $v^t \to \sqrt{1 - (v^x)^2}$ (maximum allowed tangential velocity for given normal velocity), the characteristic speeds of waves approach the local flow velocity, and wave propagation becomes severely restricted.

---

## Basic Equations of Special Relativistic Hydrodynamics

### Four-Velocity and Normalization

The four-velocity is **[Pons Eq. 2.4, translated]**:
```math
u^\mu = \gamma(1, v^x, v^y, v^z)
```

with normalization condition **[Pons Eq. 2.3]**:
```math
u^\mu u_\mu = -1
```

The Lorentz factor is **[SRGSPH Eq. 111]**, **[Pons Eq. 2.5, $W \to \gamma$]**:
```math
\gamma = \frac{1}{\sqrt{1 - v^2/c^2}} = \frac{1}{\sqrt{1 - v^2}}
```

where the total velocity squared is **[Pons Eq. 2.6]**:
```math
v^2 = (v^x)^2 + (v^y)^2 + (v^z)^2
```

### Lagrangian Form

The Lagrangian derivative is **[SRGSPH Eq. 100]**:
```math
\frac{d}{dt} \equiv \partial_t + \mathbf{v} \cdot \nabla
```

**Equation of Continuity** **[SRGSPH Eq. 90]**:
```math
\frac{dN}{dt} = -N \nabla \cdot \mathbf{v}
```

**Equation of Motion** **[SRGSPH Eq. 92]**:
```math
\frac{d\mathbf{S}}{dt} = -\frac{1}{N} \nabla P
```

**Equation of Energy** **[SRGSPH Eq. 94]**:
```math
\frac{de}{dt} = -\frac{1}{N} \nabla \cdot (P\mathbf{v})
```

### Relativistic Canonical Quantities

The relativistic canonical momentum per baryon **[SRGSPH Eq. 104]**:
```math
\mathbf{S} = \gamma H \mathbf{v}
```

The canonical energy per baryon in the lab frame **[SRGSPH Eq. 107]**:
```math
e = \gamma H - \frac{P}{Nc^2}
```

The enthalpy per baryon **[SRGSPH Eq. 112]**:
```math
H = 1 + \frac{u}{c^2} + \frac{P}{nc^2}
```

The specific enthalpy **[Pons Eq. 2.2]**:
```math
h = 1 + \epsilon + \frac{P}{\rho}
```

Relation between lab frame and rest frame densities **[SRGSPH Eq. 114]**:
```math
N = \gamma n
```

### Equation of State

For an ideal gas **[SRGSPH Eq. 119]**, **[Pons Eq. 2.8, $\gamma \to \gamma_c$]**:
```math
P = (\gamma_c - 1) n u
```

or equivalently:
```math
P = (\gamma_c - 1) \rho \epsilon
```

The sound speed is defined by **[Pons Eq. 2.9]**:
```math
h c_s^2 = \left. \frac{\partial P}{\partial \rho} \right|_s
```

For an ideal gas, this gives **[SRGSPH, after Eq. 394]**:
```math
c_s = \sqrt{\frac{(\gamma_c - 1)(H - 1)}{H}}
```

---

## Conservative Formulation

### Conservation Laws

The equations in conservative form **[Pons Eq. 2.7]**:
```math
\mathbf{U}_{,t} + \mathbf{F}^{(i)}_{,i} = 0
```

where the vector of conserved variables is **[Pons Eq. 2.10]**:
```math
\mathbf{U} = (D, S^1, S^2, S^3, \tau)^T
```

and the flux vector is **[Pons Eq. 2.11]**:
```math
\mathbf{F}^{(i)} = (D v^i, S^1 v^i + P \delta^{1i}, S^2 v^i + P \delta^{2i}, S^3 v^i + P \delta^{3i}, S^i - D v^i)^T
```

### Conserved Variables

The conserved variables in terms of primitive variables **[Pons Eq. 2.10, translated]**:

**Rest-mass density:**
```math
D = \rho \gamma
```

**Momentum density:**
```math
S^i = \rho h \gamma^2 v^i
```

**Energy density:**
```math
\tau = \rho h \gamma^2 - P - D
```

### Density Current and Energy-Momentum Tensor

The density current **[Pons Eq. 2.1]**:
```math
J^\mu = \rho u^\mu
```

The energy-momentum tensor for a perfect fluid **[Pons Eq. 2.2]**:
```math
T^{\mu\nu} = \rho h u^\mu u^\nu + P \eta^{\mu\nu}
```

Conservation equations **[Pons Eqs. 2.4, 2.5]**:
```math
J^\mu_{,\mu} = 0
```

```math
T^{\mu\nu}_{,\mu} = 0
```

---

## SRGSPH Formulation with Fixed Smoothing Length

### Kernel Function

The Gaussian kernel (in $d$ dimensions) **[SRGSPH Eq. 138]**:
```math
W(\mathbf{x}, h) = \left[ \frac{1}{h\sqrt{\pi}} \right]^d \exp\left[ -\frac{\mathbf{x}^2}{h^2} \right]
```

**Note**: $W(\mathbf{x}, h)$ is the SPH kernel function, distinct from the Lorentz factor $\gamma$.

### Number Density Field

**[SRGSPH Eq. 134]**:
```math
N(\mathbf{x}) \equiv \nu \sum_j W(\mathbf{x} - \mathbf{x}_j, h)
```

where $\nu$ is the baryon number in an SPH particle (constant in this section).

### Convolution for Physical Quantities

For particle $i$ **[SRGSPH Eq. 146]**:
```math
\langle N_i \rangle \equiv N(\mathbf{x}_i)
```

**[SRGSPH Eq. 148]**:
```math
\langle f_i \rangle \equiv \int f(\mathbf{x}) W(\mathbf{x} - \mathbf{x}_i, h) d\mathbf{x}
```

### Accuracy of Convolution

Taylor expansion around $\mathbf{x}_i$ **[SRGSPH Eq. 154]**:
```math
\langle f_i \rangle = f(\mathbf{x}_i) + \frac{h^2}{4} f''(\mathbf{x}_i) + \mathcal{O}(h^4)
```

The error is $\mathcal{O}(h^2)$.

### Derivative of Kernel

Important property **[SRGSPH Eq. 171]**:
```math
\nabla W(\mathbf{x} - \mathbf{x}_j, h) = -\nabla_j W(\mathbf{x} - \mathbf{x}_j, h)
```

### Equation of Motion (Fixed h)

Starting from the Lagrangian form **[SRGSPH Eq. 162]**:
```math
\langle \dot{\mathbf{S}}_i \rangle \equiv \int \frac{d\mathbf{S}(\mathbf{x})}{dt} W(\mathbf{x} - \mathbf{x}_i, h) d\mathbf{x}
```

After manipulation through **[SRGSPH Eqs. 163-166]**:
```math
\langle \dot{\mathbf{S}}_i \rangle = -\sum_j \int \frac{\nu P(\mathbf{x})}{N^2(\mathbf{x})} [\nabla_i - \nabla_j] W(\mathbf{x} - \mathbf{x}_i, h) W(\mathbf{x} - \mathbf{x}_j, h) d\mathbf{x}
```

### Equation of Energy (Fixed h)

**[SRGSPH Eqs. 176-181]**:
```math
\langle \dot{e}_i \rangle = -\sum_j \int \frac{\nu P(\mathbf{x}) \mathbf{v}(\mathbf{x})}{N^2(\mathbf{x})} \cdot [\nabla_i - \nabla_j] W(\mathbf{x} - \mathbf{x}_i, h) W(\mathbf{x} - \mathbf{x}_j, h) d\mathbf{x}
```

### Evaluation of Convolution Integrals

The key transformation uses the product of two Gaussians **[SRGSPH §2.2.3]**:

**[SRGSPH Eqs. 188-191]**:
```math
\int \frac{f(\mathbf{x})}{N^2(\mathbf{x})} [\nabla_i - \nabla_j] W(\mathbf{x} - \mathbf{x}_i, h) W(\mathbf{x} - \mathbf{x}_j, h) d\mathbf{x}
```
```math
= f_{ij} V_{ij}^2 [\nabla_i W(\mathbf{x}_i - \mathbf{x}_j, \sqrt{2}h) - \nabla_j W(\mathbf{x}_i - \mathbf{x}_j, \sqrt{2}h)]
```

where:

**Effective Volume Factor** **[SRGSPH Eq. 194]**:
```math
V_{ij}^2(h) \equiv \int \frac{1}{N^2(\mathbf{x})} \left( \frac{\sqrt{2}}{h\sqrt{\pi}} \right)^d \exp\left[ -\frac{2}{h^2} \left( \mathbf{x} - \frac{\mathbf{x}_i + \mathbf{x}_j}{2} \right)^2 \right] d\mathbf{x}
```

**Weighted Average** **[SRGSPH Eq. 196]**:
```math
f_{ij}(h) \equiv \frac{1}{V_{ij}^2} \int \frac{f(\mathbf{x})}{N^2(\mathbf{x})} \left( \frac{\sqrt{2}}{h\sqrt{\pi}} \right)^d \exp\left[ -\frac{2}{h^2} \left( \mathbf{x} - \frac{\mathbf{x}_i + \mathbf{x}_j}{2} \right)^2 \right] d\mathbf{x}
```

### Approximations

**[SRGSPH, discussion after Eq. 200]:**

The volume factor $V_{ij}^2$ cannot be solved analytically due to the SPH particle distribution, so it is approximated by interpolation:
```math
V_{ij}^2 \approx V_{ij,\text{interp}}^2
```

**[SRGSPH, discussion before Eq. 209]:**

The weighted average is replaced by the Riemann solution to accurately capture shocks:
```math
f_{ij} \approx f_{ij}^* \quad \text{(from Riemann solver)}
```

### Final SRGSPH Equations (Fixed h)

**Equation of Motion** **[SRGSPH Eq. 209]**:
```math
\langle \dot{\mathbf{S}}_i \rangle = -\sum_j \nu P_{ij}^* V_{ij,\text{interp}}^2 [\nabla_i W(\mathbf{x}_i - \mathbf{x}_j, \sqrt{2}h) - \nabla_j W(\mathbf{x}_i - \mathbf{x}_j, \sqrt{2}h)]
```

**Equation of Energy** **[SRGSPH Eq. 212]**:
```math
\langle \dot{e}_i \rangle = -\sum_j \nu P_{ij}^* \mathbf{v}_{ij}^* \cdot V_{ij,\text{interp}}^2 [\nabla_i W(\mathbf{x}_i - \mathbf{x}_j, \sqrt{2}h) - \nabla_j W(\mathbf{x}_i - \mathbf{x}_j, \sqrt{2}h)]
```

where $P_{ij}^*$ and $\mathbf{v}_{ij}^*$ are obtained from the Riemann solver **[SRGSPH Eq. 423]**.

**Conservation property** **[SRGSPH, discussion after Eq. 214]**: These equations satisfy anti-symmetry in indices $i$ and $j$, ensuring conservation of total momentum and energy.

---

## SRGSPH Formulation with Variable Smoothing Length

### Why Volume-Based Approach?

**[SRGSPH §2.4, §2.5]**

Traditional SPH defines number density first. With variable smoothing lengths, this can be done via:
- **Gather**: $N(\mathbf{x}) = \nu \sum_j W(\mathbf{x} - \mathbf{x}_j, h(\mathbf{x}))$
- **Scatter**: $N(\mathbf{x}) = \nu \sum_j W(\mathbf{x} - \mathbf{x}_j, h_j)$

Both have issues:
- Gather introduces $\nabla h$ terms requiring smoothing
- Scatter causes overshooting at discontinuities

With non-uniform baryon numbers, the "standard" approach:
```math
N_{\text{standard}}(\mathbf{x}) = \sum_j \nu_j W(\mathbf{x} - \mathbf{x}_j, h(\mathbf{x}))
```
causes $h$ to become sharply varying at density jumps, increasing numerical errors and reducing time steps.

**Solution**: Define particle volume first, ensuring smooth $h$ variations.

### Particle Volume Definition

**[SRGSPH Eq. 221]**:
```math
V_p(\mathbf{x}) \equiv \left[ \sum_j W(\mathbf{x} - \mathbf{x}_j, h(\mathbf{x})) \right]^{-1}
```

### Smoothing Length Definition

The smoothing length $h(\mathbf{x})$ is defined over the entire spatial domain **[SRGSPH Eq. 231]**:
```math
h(\mathbf{x}) \equiv \eta [V_p^*(\mathbf{x})]^{1/d}
```

where **[SRGSPH Eq. 233]**:
```math
V_p^*(\mathbf{x}) \equiv \left[ \sum_j W(\mathbf{x} - \mathbf{x}_j, C_{\text{smooth}} h(\mathbf{x})) \right]^{-1}
```

**Typical values** **[SRGSPH, after Eq. 236]**: $\eta = 1.0$, $C_{\text{smooth}} = 2.0$

**Purpose of $C_{\text{smooth}} > 1$**: This "pre-smooths" the volume field, reducing $\nabla h$ so it can be neglected in the equations of motion/energy.

### Iterative Solution

**[SRGSPH, discussion after Eq. 238]**:

The equations for $h(\mathbf{x})$ and $V_p^*(\mathbf{x})$ are solved iteratively until convergence:
- First time step: iterate Eqs. (231) and (233) until convergence
- Subsequent steps: use previous $h(\mathbf{x})$ as initial guess for $V_p^*(\mathbf{x})$

### Number Density with Variable Baryon Number

The volume-based approach (preferred) **[SRGSPH Eq. 243]**:
```math
N(\mathbf{x}) = \nu(\mathbf{x}) V_p^{-1}(\mathbf{x})
```

where $\nu(\mathbf{x})$ is the spatially-varying baryon number.

The standard approach (for comparison) **[SRGSPH Eq. 285]**:
```math
N_{\text{standard}}(\mathbf{x}) \equiv \sum_j \nu_j W(\mathbf{x} - \mathbf{x}_j, h(\mathbf{x}))
```

### Convolution with Variable h

For particle $i$ **[SRGSPH Eq. 249]**:
```math
\langle N_i \rangle \equiv N(\mathbf{x}_i)
```

**[SRGSPH Eq. 252]**:
```math
\langle \nu_i f_i \rangle \equiv \int \nu(\mathbf{x}) f(\mathbf{x}) W(\mathbf{x} - \mathbf{x}_i, h(\mathbf{x})) d\mathbf{x}
```

### Kernel Gradient with Variable h

**[SRGSPH Eq. 332]**:
```math
\nabla W(\mathbf{x} - \mathbf{x}_i, h(\mathbf{x})) = -\left[ \nabla_i - \nabla h \left( \frac{h}{2} \nabla_i^2 + \frac{1-d}{h} \right) \right] W(\mathbf{x} - \mathbf{x}_i, h(\mathbf{x}))
```

### Equation of Motion (Variable h)

**[SRGSPH Eqs. 323-328]**:
```math
\langle \nu_i \dot{\mathbf{S}}_i \rangle = -\sum_j \int P(\mathbf{x}) V_p^2(\mathbf{x}) \left[ (\nabla_i - \nabla_j) - \frac{h}{2} \nabla h (\nabla_i^2 - \nabla_j^2) \right] W(\mathbf{x} - \mathbf{x}_i, h(\mathbf{x})) W(\mathbf{x} - \mathbf{x}_j, h(\mathbf{x})) d\mathbf{x}
```

### Equation of Energy (Variable h)

**[SRGSPH Eqs. 338-344]**:
```math
\langle \nu_i \dot{e}_i \rangle = -\sum_j \int P(\mathbf{x}) \mathbf{v}(\mathbf{x}) \cdot V_p^2(\mathbf{x}) \left[ (\nabla_i - \nabla_j) - \frac{h}{2} \nabla h (\nabla_i^2 - \nabla_j^2) \right] W(\mathbf{x} - \mathbf{x}_i, h(\mathbf{x})) W(\mathbf{x} - \mathbf{x}_j, h(\mathbf{x})) d\mathbf{x}
```

### Approximations for Variable h

**[SRGSPH §2.6.1, Eqs. 352-353]**

Assuming smoothing length varies smoothly (thanks to $C_{\text{smooth}} > 1$):
```math
\nabla h \approx 0
```

For preserving symmetry between particles $i$ and $j$:
```math
(\nabla_i - \nabla_j) \int W(\mathbf{x}_i - \mathbf{x}_j, \sqrt{2}h(\mathbf{x})) d\mathbf{x} \approx \nabla_i W(\mathbf{x}_i - \mathbf{x}_j, \sqrt{2}h_i) - \nabla_j W(\mathbf{x}_i - \mathbf{x}_j, \sqrt{2}h_j)
```

### Averaged Volume Factor

**[SRGSPH Eq. 365]**:
```math
V_{ij}^2 \equiv \frac{1}{2} \left( V_{ij}^2(h_i) + V_{ij}^2(h_j) \right)
```

This symmetric averaging ensures action-reaction principle.

### Final SRGSPH Equations (Variable h)

**Equation of Motion** **[SRGSPH Eq. 371]**:
```math
\langle \nu_i \dot{\mathbf{S}}_i \rangle = -\sum_j P_{ij}^* V_{ij,\text{interp}}^2 [\nabla_i W(\mathbf{x}_i - \mathbf{x}_j, \sqrt{2}h_i) - \nabla_j W(\mathbf{x}_i - \mathbf{x}_j, \sqrt{2}h_j)]
```

**Equation of Energy** **[SRGSPH Eq. 373]**:
```math
\langle \nu_i \dot{e}_i \rangle = -\sum_j P_{ij}^* \mathbf{v}_{ij}^* \cdot V_{ij,\text{interp}}^2 [\nabla_i W(\mathbf{x}_i - \mathbf{x}_j, \sqrt{2}h_i) - \nabla_j W(\mathbf{x}_i - \mathbf{x}_j, \sqrt{2}h_j)]
```

**Conservation property** **[SRGSPH, discussion after Eq. 378]**: Anti-symmetry is preserved even with different smoothing lengths $h_i \neq h_j$, ensuring conservation.

---

## Theory of Relativistic Riemann Problems

### The Riemann Problem Structure

**[Pons §1, §5]**

The Riemann problem is the initial value problem of a discontinuous piecewise constant initial state separating two uniform regions $L$ (left) and $R$ (right). At $t = 0$, a discontinuity exists at $x = 0$.

**Classical (Newtonian) case:**
- Tangential velocities remain constant across all waves
- The problem can be solved using only normal velocity components
- Solution depends only on 4 quantities: $P_L, \rho_L, P_R, \rho_R$ (with $v_L^x, v_R^x$)

**Relativistic case with tangential velocities:**
- Tangential velocities **change** across rarefactions and shocks
- Cannot transform to a frame where all tangential velocities vanish
- Solution depends on **all 6 velocity components**: $v_L^{x,y,z}, v_R^{x,y,z}$
- Thermodynamic coupling through enthalpy affects even slow flows

### Wave Structure and Self-Similarity

**[Pons §5]**

The solution consists of three elementary waves:
```math
I \rightarrow L \; \mathcal{W}_\leftarrow \; L_* \; \mathcal{C} \; R_* \; \mathcal{W}_\rightarrow \; R
```

where:
- $\mathcal{W}_\leftarrow, \mathcal{W}_\rightarrow$: Left- and right-going waves (shock or rarefaction)
- $\mathcal{C}$: Contact discontinuity (moves with the fluid)
- $L_*, R_*$: Intermediate constant states

**Self-similarity**: All waves propagate as $\xi = x/t$ = constant. This allows reduction from PDE to ODE for rarefactions and algebraic jump conditions for shocks.

### Contact Discontinuity Properties

**[Pons §5]**

Across the contact discontinuity:
- **Continuous**: Pressure $P$, normal velocity $v^x$
- **Discontinuous**: Density $\rho$, tangential velocities $v^{y,z}$, internal energy $\epsilon$

This is fundamentally different from Newtonian case where tangential velocities are continuous.

### Key Difference: Riemann Invariants

**[Pons §3, Eqs. 3.9-3.10, translated $W \to \gamma$]**

In rarefaction waves and across shocks:
```math
h \gamma v^y = \text{constant}
```
```math
h \gamma v^z = \text{constant}
```

This means:
- Direction of tangential velocity is preserved: $v^y/v^z = \text{const}$
- Magnitude changes due to changes in $h$ and $\gamma$

**Physical interpretation**: As the fluid expands/compresses, the "relativistic tangential momentum" $h \gamma v^t$ is conserved, but the actual tangential velocity adjusts.

---

## Riemann Problem: Rarefaction Waves

### Self-Similarity and Physical Picture

**[Pons §3]**

Rarefaction waves are **simple waves** where fluid elements smoothly transition from one state to another. In the self-similar solution, each fluid element moves at constant $\xi = x/t$.

**Physical process:**
- Fluid elements continuously expand (rarefaction) or compress
- Pressure and density monotonically decrease/increase
- Entropy remains constant along fluid lines (isentropic process)
- All changes occur smoothly (no discontinuities)

### Self-Similarity Variable

**[Pons, before Eq. 3.1]**:

Rarefaction waves depend only on:
```math
\xi = \frac{x}{t}
```

### Self-Similarity Equations

From the continuity equation **[Pons Eq. 3.1, translated $W \to \gamma$]**:
```math
(v^x - \xi) \frac{d\rho}{d\xi} + \{\rho \gamma^2 v^x (v^x - \xi) + \rho\} \frac{dv^x}{d\xi} + \rho \gamma^2 v^y (v^x - \xi) \frac{dv^y}{d\xi} + \rho \gamma^2 v^z (v^x - \xi) \frac{dv^z}{d\xi} = 0
```

From $x$-momentum conservation **[Pons Eq. 3.2, translated]**:
```math
\rho h \gamma^2 (v^x - \xi) \frac{dv^x}{d\xi} + (1 - v^x \xi) \frac{dP}{d\xi} = 0
```

From $y$-momentum conservation **[Pons Eq. 3.3, translated]**:
```math
\rho h \gamma^2 (v^x - \xi) \frac{dv^y}{d\xi} - v^y \xi \frac{dP}{d\xi} = 0
```

From $z$-momentum conservation **[Pons Eq. 3.4, translated]**:
```math
\rho h \gamma^2 (v^x - \xi) \frac{dv^z}{d\xi} - v^z \xi \frac{dP}{d\xi} = 0
```

From isentropic condition (entropy conservation along fluid lines) **[Pons Eq. 3.5]**:
```math
\frac{dP}{d\xi} = h c_s^2 \frac{d\rho}{d\xi} = \rho \frac{dh}{d\xi}
```

**Note on Eqs. (3.3-3.4)** **[Pons, discussion after Eq. 3.4]**: If tangential velocities vanish in the initial state, they remain zero throughout the rarefaction. This recovers the Newtonian behavior.

### Characteristic Speeds (Eigenvalues)

**[Pons Eq. 3.6]**

Non-trivial similarity solutions exist only when:
```math
\xi = \frac{v^x(1 - c_s^2) \pm c_s \sqrt{(1 - v^2)[1 - v^2 c_s^2 - (v^x)^2(1 - c_s^2)]}}{1 - v^2 c_s^2}
```

The $+$ sign corresponds to right-going waves ($\mathcal{R}_\rightarrow$), the $-$ sign to left-going waves ($\mathcal{R}_\leftarrow$).

**Physical meaning** **[Pons, discussion after Eq. 3.6]**: These are the maximum and minimum eigenvalues of the Jacobian matrix $\partial \mathbf{F}^{(x)}/\partial \mathbf{U}$. They represent the characteristic speeds of wave propagation.

**Effect of tangential velocity**: The presence of $v^2 = (v^x)^2 + (v^y)^2 + (v^z)^2$ in the eigenvalues shows that tangential velocities affect wave propagation speed.

### Reduced System for Rarefaction Waves

**[Pons Eqs. 3.7-3.9, translated]**

The system reduces to one ODE and two algebraic constraints:

**ODE:**
```math
\rho h \gamma^2 (v^x - \xi) dv^x + (1 - \xi v^x) dP = 0
```

**Algebraic constraints (Riemann invariants):**
```math
h \gamma v^y = \text{constant}
```

```math
h \gamma v^z = \text{constant}
```

**Physical interpretation** **[Pons, discussion after Eq. 3.9]**:
- The tangential "relativistic momentum per baryon" is conserved
- Direction $v^y/v^z$ is constant across the rarefaction
- Only magnitude adjusts due to enthalpy and Lorentz factor changes

**Limit to classical case** **[Pons, discussion after Eq. 3.9]**: In the Newtonian limit ($v \ll 1$), $\gamma = 1$, but Eqs. (3.8-3.9) do **not** reduce to $v^{y,z} = \text{const}$ because enthalpy $h$ still varies. Thus, even for slow flows, thermodynamically relativistic situations ($h > 1$) require the full treatment.

### ODE in Standard Form

**[Pons Eq. 3.10, translated]**:
```math
\frac{dv^x}{dP} = \pm \frac{1}{\rho h \gamma^2 c_s} \frac{1}{\sqrt{1 + g(\xi_\pm, v^x, v^t)}}
```

where **[Pons, after Eq. 3.10]**:
```math
g(\xi_\pm, v^x, v^t) = \frac{(v^t)^2 (\xi_\pm^2 - 1)}{(1 - \xi_\pm v^x)^2}
```

and $v^t = \sqrt{(v^y)^2 + (v^z)^2}$ is the tangential velocity magnitude.

**When $v^t = 0$**: The function $g = 0$ and the equation reduces to **[Pons, after Eq. 3.10, translated]**:
```math
\gamma^2 dv^x = \pm \frac{c_s}{\gamma_c P} dP = \pm \frac{c_s}{\rho} d\rho
```
recovering the purely normal flow result.

### Tangential Velocity in Post-Rarefaction State

**[Pons Eq. 3.11, translated]**

Given the constraint $h_a \gamma_a v_a^t = h_b \gamma_b v_b^t$ and solving for $v_b^t$:
```math
v_b^t = h_a \gamma_a v_a^t \left\{ \frac{1 - (v_b^x)^2}{h_b^2 + (h_a \gamma_a v_a^t)^2} \right\}^{1/2}
```

**Derivation**: From $v_b^t = (h_a \gamma_a v_a^t)/(h_b \gamma_b)$ and $\gamma_b^2 = 1/(1 - (v_b^x)^2 - (v_b^t)^2)$, solve for $v_b^t$.

### Rarefaction Function

**[Pons, discussion after Eq. 3.11]**

The integration of the ODE from state $a$ (ahead of wave) to state $b$ (behind wave) gives:
```math
v_b^x = \mathcal{R}_{\rightleftharpoons}^a(P_b)
```

This function connects the post-rarefaction normal velocity to the pressure, given the pre-rarefaction state $a$.

**Integration details**:
- Use Eq. (3.10) with initial condition $v^x(P_a) = v_a^x$
- At each pressure step, compute $\xi_\pm$ from Eq. (3.6)
- Update $v^t$ using Eq. (3.11)
- Update thermodynamic variables using EOS
- Integrate to desired $P_b$

**Physical behavior** **[Pons, discussion of Figure 1]**:
- Rarefactions moving toward state $a$: $P_b < P_a$
- Rarefactions moving away from state $a$: $P_b > P_a$
- Tangential velocity restricts achievable normal velocities
- As $v^t \to \sqrt{1 - (v^x)^2}$, the range of $v^x$ in the rarefaction shrinks

---

## Riemann Problem: Shock Waves

### Physical Picture of Shocks

**[Pons §4, intro]**

Shock waves are **discontinuities** where fluid quantities jump instantaneously. Unlike rarefactions:
- Changes occur across a surface (zero thickness)
- Entropy increases (non-isentropic)
- Governed by jump conditions (Rankine-Hugoniot), not differential equations
- Mass, momentum, and energy fluxes are continuous

**Shock types**:
- **Compressive**: Density and pressure increase ($\rho_b > \rho_a$, $P_b > P_a$)
- All physical shocks in ideal gas are compressive

### Rankine-Hugoniot Jump Conditions

**[Pons Eqs. 4.1-4.2, translated]**

For a shock normal to the $x$-axis, with normal vector:
```math
n^\nu = \gamma_s(V_s, 1, 0, 0)
```

where $V_s$ is the shock velocity **[Pons Eq. 4.4]** and:
```math
\gamma_s = \frac{1}{\sqrt{1 - V_s^2}}
```
is the shock Lorentz factor **[Pons Eq. 4.5, $W_s \to \gamma_s$]**.

The jump conditions **[Pons Eqs. 4.1-4.2]**:
```math
[\rho u^\mu] n_\mu = 0
```

```math
[T^{\mu\nu}] n_\nu = 0
```

where $[F] = F_a - F_b$ denotes the jump across the shock **[Pons Eq. 4.3]**.

**Physical meaning**:
- First equation: Rest-mass flux is continuous
- Second equation: Energy-momentum flux is continuous

### Invariant Mass Flux

**[Pons Eq. 4.7, translated]**:
```math
j \equiv \gamma_s D_a (V_s - v_a^x) = \gamma_s D_b (V_s - v_b^x)
```

**Convention** **[Pons, note after Eq. 4.7]**: $j > 0$ for right-going shocks, $j < 0$ for left-going shocks.

**Physical meaning**: $j$ is the rest-mass flux in the lab frame, an invariant quantity that characterizes the shock strength.

### Jump Conditions in Terms of Conserved Variables

**Normal velocity jump** **[Pons Eq. 4.8, translated]**:
```math
[v^x] = -\frac{j}{\gamma_s} \left[ \frac{1}{D} \right]
```

**Pressure jump** **[Pons Eq. 4.9, translated]**:
```math
[P] = \frac{j}{\gamma_s} \left[ \frac{S^x}{D} \right]
```

**Tangential momentum per unit mass (continuous)** **[Pons Eqs. 4.10-4.11]**:
```math
\left[ \frac{S^y}{D} \right] = 0
```

```math
\left[ \frac{S^z}{D} \right] = 0
```

**Energy-pressure relation** **[Pons Eq. 4.12, translated]**:
```math
[v^x P] = \frac{j}{\gamma_s} \left[ \frac{\tau}{D} \right]
```

**Note** **[Pons, discussion after Eq. 4.12]**: These are derived from Eqs. (4.1-4.2) using the non-zero mass flux. For tangential discontinuities ($j = 0$), only pressure and normal velocity continuity remain.

### Tangential Velocity Invariant

**[Pons Eqs. 4.10-4.11, and discussion, translated]**

From $S^{y,z}/D = h \gamma v^{y,z} = \text{const}$:
```math
h \gamma v^y = \text{constant across shock}
```

```math
h \gamma v^z = \text{constant across shock}
```

**Physical interpretation** **[Pons, discussion after Eq. 4.11]**:
- Same invariant as in rarefaction waves
- Tangential velocity direction is preserved
- Magnitude adjusts due to jumps in $h$ and $\gamma$

### Post-Shock Normal Velocity

**[Pons, derivation after Eq. 4.12, translated]**

From the jump conditions and using $S^x = (\tau + P + D)v^x$:
```math
v_b^x = \frac{h_a \gamma_a v_a^x + \frac{\gamma_s(P_b - P_a)}{j}}{h_a \gamma_a + (P_b - P_a)\left(\frac{\gamma_s v_a^x}{j} + \frac{1}{\rho_a \gamma_a}\right)}
```

**Form** **[Pons, discussion after this equation, translated]**: Appears similar to zero tangential velocity case, but $\gamma_a$ contains the tangential velocity through $\gamma_a = 1/\sqrt{1 - v_a^2}$ where $v_a^2 = (v_a^x)^2 + (v_a^y)^2 + (v_a^z)^2$.

### Post-Shock Tangential Velocities

**[Pons Eq. 4.13, translated]**:
```math
v_b^{y,z} = h_a \gamma_a v_a^{y,z} \left[ \frac{1 - (v_b^x)^2}{h_b^2 + (h_a \gamma_a v_a^{y,z})^2} \right]^{1/2}
```

**Derivation**: From $h_a \gamma_a v_a^{y,z} = h_b \gamma_b v_b^{y,z}$ and $\gamma_b^2 = 1/(1 - (v_b^x)^2 - (v_b^y)^2 - (v_b^z)^2)$, solve for $v_b^{y,z}$.

### Shock Velocity

**[Pons Eq. 4.14, translated]**

From the mass flux definition $j = \gamma_s D_a (V_s - v_a^x)$, solving for $V_s$:
```math
V_s^\pm = \frac{\rho_a^2 \gamma_a^2 v_a^x \pm |j| \sqrt{j^2 + \rho_a^2 \gamma_a^2 (1 - (v_a^x)^2)}}{\rho_a^2 \gamma_a^2 + j^2}
```

where $V_s^+$ is for right-going shocks and $V_s^-$ for left-going shocks.

### Taub Adiabat

**[Pons, discussion before Eq. 4.15]**

The **Taub adiabat** (relativistic Hugoniot adiabat) relates thermodynamic quantities across the shock.

**Derivation**: Multiply Eq. (4.2) by $(h u_\mu)_a$ and $(h u_\mu)_b$ separately, then add.

**Result** **[Pons Eq. 4.15]**:
```math
[h^2] = \left( \frac{h_b}{\rho_b} + \frac{h_a}{\rho_a} \right) [P]
```

For an ideal gas with constant $\gamma_c$, eliminate $\rho_b$ using EOS to get **[Pons Eq. 4.16, $\gamma \to \gamma_c$]**:
```math
h_b^2 \left(1 + \frac{(\gamma_c - 1)(P_a - P_b)}{\gamma_c P_b}\right) - \frac{(\gamma_c - 1)(P_a - P_b)}{\gamma_c P_b} h_b + \frac{h_a(P_a - P_b)}{\rho_a} - h_a^2 = 0
```

This is a quadratic equation for $h_b$ **[Pons, discussion after Eq. 4.16]**. One root is always negative and must be discarded; the positive root is the physical solution.

### Mass Flux from Pressure

**[Pons Eq. 4.17]**

Multiplying Eq. (4.2) by $n_\mu$ and using the mass flux definition:
```math
j^2 = \frac{-[P]}{[h/\rho]}
```

**Algorithm** **[Pons, discussion after Eq. 4.17]**:
1. Given $P_b$, solve Taub adiabat (Eq. 4.16) for $h_b$
2. Use EOS to get $\rho_b$
3. Compute $[h/\rho]$ and $[P]$
4. Obtain $j^2(P_b)$

This, combined with earlier equations, gives all post-shock quantities as functions of $P_b$ alone.

### Shock Function

**[Pons, discussion after Eq. 4.17]**

The relation between post-shock normal velocity and pressure:
```math
v_b^x = \mathcal{S}_{\rightleftharpoons}^a(P_b)
```

**Physical behavior** **[Pons, discussion of Figure 2]**:
- Shocks moving toward state $a$: $P_b > P_a$ (compressive)
- Tangential velocity restricts post-shock normal velocity
- Larger $v^t$ leads to smaller achievable $v_b^x$
- Wave speeds $V_s$ approach $v_a^x$ as $v_a^t \to \sqrt{1 - (v_a^x)^2}$

---

## Complete Riemann Solver Algorithm

### Problem Setup

**[Pons §1]**

Initial left and right states:
- Left: $(P_L, \rho_L, v_L^x, v_L^y, v_L^z)$
- Right: $(P_R, \rho_R, v_R^x, v_R^y, v_R^z)$

Discontinuity at $x = 0$ at $t = 0$.

**Goal**: Find intermediate states $L_*$ and $R_*$ separated by a contact discontinuity.

### Wave Structure

**[Pons §5, after Eq. 5.1]**

The decay produces three waves:
```math
I \rightarrow L \; \mathcal{W}_\leftarrow \; L_* \; \mathcal{C} \; R_* \; \mathcal{W}_\rightarrow \; R
```

where:
- $\mathcal{W}$: shock ($\mathcal{S}$) or rarefaction ($\mathcal{R}$)
- $\mathcal{C}$: contact discontinuity

**Physical picture** **[Pons §5, intro]**:
- Left wave connects $L$ to $L_*$, propagates toward $L$
- Right wave connects $R$ to $R_*$, propagates toward $R$
- Contact moves with local fluid velocity $v_*^x$

### Wave Discrimination

**[Pons Eq. 5.2, $p \to P$]**

For a wave connecting states $a$ (ahead) and $b$ (behind):
```math
\mathcal{W} = \begin{cases}
\mathcal{R} & \text{if } P_b \leq P_a \\
\mathcal{S} & \text{if } P_b > P_a
\end{cases}
```

**Physical basis**:
- Shocks are compressive: $P_b > P_a$, $\rho_b > \rho_a$
- Rarefactions are expansive: $P_b < P_a$, $\rho_b < \rho_a$
- This criterion works for ideal gas and causal EOS

### Jump Conditions at Contact Discontinuity

**[Pons Eqs. 1.1-1.2, $p \to P$]**:

```math
P_{L*} = P_{R*} = P_*
```

```math
v_{L*}^x = v_{R*}^x = v_*^x
```

The tangential velocities $v_{L*}^{y,z} \neq v_{R*}^{y,z}$ and densities $\rho_{L*} \neq \rho_{R*}$ are **discontinuous** across the contact.

**Physical meaning** **[Pons §5]**: The contact is a material interface moving with the fluid. Pressure must be continuous (force balance), and normal velocity must be continuous (no flow across), but tangential velocity and thermodynamic state can jump.

### Solution Procedure

**[Pons §5]**

#### 1. Define composite functions

For the left wave ($L \rightarrow L_*$):
```math
v_{L*}^x(P) = \begin{cases}
\mathcal{R}_\leftarrow^L(P) & \text{if } P \leq P_L \\
\mathcal{S}_\leftarrow^L(P) & \text{if } P > P_L
\end{cases}
```

For the right wave ($R \rightarrow R_*$):
```math
v_{R*}^x(P) = \begin{cases}
\mathcal{R}_\rightarrow^R(P) & \text{if } P \leq P_R \\
\mathcal{S}_\rightarrow^R(P) & \text{if } P > P_R
\end{cases}
```

**Physical interpretation** **[Pons §5]**: Each function maps the intermediate pressure to the corresponding normal velocity in the intermediate state, accounting for whether the wave is a shock or rarefaction.

#### 2. Solve for $P_*$

**[Pons Eq. 1.2]**

Find the root of:
```math
v_{L*}^x(P) - v_{R*}^x(P) = 0
```

**Solution methods**:
- Newton-Raphson (requires derivatives)
- Bisection (robust, slower)
- Brent's method (recommended: robust + fast)

**Bracketing**: For ideal gas, $P_* \in [\min(P_L, P_R), \max(P_L, P_R) + \Delta P]$ where $\Delta P$ depends on velocity jumps.

#### 3. Compute $v_*^x$

Once $P_*$ is known **[Pons, discussion after Eq. 1.2]**:
```math
v_*^x = v_{L*}^x(P_*) = v_{R*}^x(P_*)
```

#### 4. Compute remaining quantities in $L_*$ and $R_*$

**[Pons §5]**

Using the appropriate wave relations (rarefaction or shock) and the EOS:
- **Densities**: $\rho_{L*}$, $\rho_{R*}$ (from Taub adiabat for shocks, isentropic relation for rarefactions)
- **Tangential velocities**: $v_{L*}^{y,z}$, $v_{R*}^{y,z}$ (from Eq. 3.11 or 4.13)
- **Enthalpies**: $h_{L*}$, $h_{R*}$ (from EOS)
- **Internal energies**: $\epsilon_{L*}$, $\epsilon_{R*}$ (from EOS)

#### 5. Determine wave positions

**For rarefactions** **[Pons §3]**:
- Head: $\xi_{\text{head}} = \xi(P_I, v_I, v_I^t)$ at state $I$ (from Eq. 3.6)
- Tail: $\xi_{\text{tail}} = \xi(P_{I*}, v_{I*}, v_{I*}^t)$ at state $I_*$ (from Eq. 3.6)
- Between head and tail: continuous variation

**For shocks** **[Pons §4]**:
- Position: $V_s$ from Eq. (4.14)
- Discontinuous jump at $x = V_s t$

#### 6. Compute Riemann solution at any point $x/t$

**[Pons §5]**

Based on the wave structure, determine which state applies at any self-similarity coordinate $\xi = x/t$:
- If $\xi < \xi_{\text{left wave}}$: state $L$
- If $\xi_{\text{left wave}} < \xi < 0$ (approx): state $L_*$
- If $0 < \xi < \xi_{\text{right wave}}$ (approx): state $R_*$
- If $\xi > \xi_{\text{right wave}}$: state $R$

(More precise: compare with wave speeds/rarefaction fan boundaries)

### Solution for Godunov Flux

**[SRGSPH Eq. 423]**

For SRGSPH, we need the solution at the interface ($x = 0$, so $\xi = 0$):

**If** $v_*^x > 0$: the contact moves right, use $L_*$ state
**If** $v_*^x < 0$: the contact moves left, use $R_*$ state
**If** $v_*^x = 0$: use either $L_*$ or $R_*$ (pressure and normal velocity are equal)

The Riemann solution provides:
```math
P_{ij}^* = P_* \quad \text{(pressure at interface)}
```

```math
\mathbf{v}_{ij}^* = (v_*^x, v_*^y, v_*^z) \quad \text{(velocity at interface)}
```

**Note**: For MUSCL reconstruction **[SRGSPH §2.7.3]**, the input states to the Riemann solver are the reconstructed states, not the particle values directly.

---

## Primitive Variable Recovery

**[SRGSPH §2.8]**

After time integration, we have updated conserved quantities $\langle \mathbf{S}_i \rangle$ and $\langle e_i \rangle$. We need to recover primitive variables $(n, \mathbf{v}, P, \epsilon)$.

### Quartic Equation for Lorentz Factor

**[SRGSPH, derivation in §2.8]**

From the relations:
```math
\mathbf{S} = \gamma H \mathbf{v} \quad \text{[Eq. 104]}
```
```math
e = \gamma H - \frac{P}{Nc^2} \quad \text{[Eq. 107]}
```

and using the EOS, we obtain:
```math
(\langle \gamma_i \rangle^2 - 1)(X \langle e_i \rangle \langle \gamma_i \rangle - 1)^2 - \langle \mathbf{S}_i \rangle^2 (X \langle \gamma_i \rangle^2 - 1)^2 = 0
```

where **[SRGSPH, after Eq. 409]**:
```math
X \equiv \frac{\gamma_c}{\gamma_c - 1}
```

This is a quartic equation in $\langle \gamma_i \rangle$ that must be solved numerically.

**Solution methods**:
- Newton-Raphson iteration (fast, requires good initial guess)
- Ferrari's formula for exact quartic roots
- Brent's method for root finding

### Velocity Recovery

**[SRGSPH Eq. 416]**

Once $\langle \gamma_i \rangle$ is obtained:
```math
\langle \mathbf{v}_i \rangle = \frac{X \langle \gamma_i \rangle^2 - 1}{\langle \gamma_i \rangle (X \langle e_i \rangle \langle \gamma_i \rangle - 1)} \langle \mathbf{S}_i \rangle
```

### Other Primitive Variables

**[SRGSPH §2.8]**

From the number density (known from volume):
```math
\langle n_i \rangle = \frac{\langle N_i \rangle}{\langle \gamma_i \rangle}
```

From the energy relation (rearranging Eq. 107):
```math
\langle H_i \rangle = \frac{X \langle e_i \rangle \langle \gamma_i \rangle - 1}{\langle \gamma_i \rangle^2}
```

From the enthalpy-pressure relation (using EOS):
```math
\langle P_i \rangle = \langle n_i \rangle (\langle H_i \rangle - 1) \frac{\gamma_c - 1}{\gamma_c}
```

From the enthalpy definition (Eq. 112):
```math
\langle u_i \rangle = \left(\langle H_i \rangle - 1 - \frac{\langle P_i \rangle}{\langle n_i \rangle}\right) c^2
```

**Physical interpretation**: This recovery process inverts the conserved-to-primitive map, which is nonlinear due to the Lorentz factor coupling.

---

## Time Integration

### Euler Method

**[SRGSPH Eqs. 425-428]**:

```math
\langle \nu_i \mathbf{S}_i \rangle^{n+1} = \langle \nu_i \mathbf{S}_i \rangle^n + \langle \nu_i \dot{\mathbf{S}}_i \rangle \Delta t
```

```math
\langle \nu_i e_i \rangle^{n+1} = \langle \nu_i e_i \rangle^n + \langle \nu_i \dot{e}_i \rangle \Delta t
```

```math
\mathbf{x}_i^{n+1} = \mathbf{x}_i^n + \langle \mathbf{v}_i \rangle \Delta t
```

**Note**: More accurate integrators (RK2, RK4) can be used, but Euler is sufficient for demonstration.

### Time Step (CFL Condition)

**[SRGSPH Eqs. 432-433]**:

```math
\Delta t = C_{\text{CFL}} \min_i \left[ \frac{h_i}{c_{s,i}} \right]
```

Typical value: $C_{\text{CFL}} = 0.3$

**Physical meaning** **[SRGSPH, discussion after Eq. 433]**: The time step must be small enough that information (traveling at sound speed $c_s$) doesn't propagate more than one smoothing length $h$ per step. This is simpler than the elaborate criteria in some SRSPH papers, but works well in practice.

### Higher-Order Methods (MUSCL)

**[SRGSPH §2.7.3]**

For second-order spatial accuracy, gradients are computed and limited.

**Gradient computation** **[SRGSPH, discussion before Eq. 387]**:

Compute $\left. \frac{\partial f}{\partial s} \right|_{\mathbf{x}_i}$ for $f = n, P, \mathbf{v}$ using SPH gradient formulas (not detailed here, use standard SPH gradients).

**Monotonicity constraint (shock detector)** **[SRGSPH Eq. 390]**:

Set gradient to zero if either condition holds:
```math
C_{\text{shock}} \mathbf{e}_{ij} \cdot (\mathbf{v}_i - \mathbf{v}_j) > \min(c_{s,i}, c_{s,j})
```

```math
\left| \log_{10}\left(\frac{P_i}{P_j}\right) \right| > C_{\text{c.d.}}
```

where $\mathbf{e}_{ij} = (\mathbf{x}_i - \mathbf{x}_j)/|\mathbf{x}_i - \mathbf{x}_j|$ **[SRGSPH Eq. 392]**.

Typical values **[SRGSPH Eq. 397]**: $C_{\text{shock}} = 3$, $C_{\text{c.d.}} = 1$

**Physical meaning**:
- First condition: Detects velocity convergence (potential shock)
- Second condition: Detects large pressure jumps (shock or contact discontinuity)
- Third condition (in Eq. 399): Also switches to first order if pressure or density become negative

**Reconstructed states for Riemann problem**:

At the interface between particles $i$ and $j$:
```math
f_{i,\text{Riemann}} = f_i + \left. \frac{\partial f}{\partial s} \right|_{\mathbf{x}_i} \cdot \frac{\mathbf{x}_j - \mathbf{x}_i}{2}
```

```math
f_{j,\text{Riemann}} = f_j - \left. \frac{\partial f}{\partial s} \right|_{\mathbf{x}_j} \cdot \frac{\mathbf{x}_j - \mathbf{x}_i}{2}
```

These reconstructed states are used as input to the Riemann solver instead of the cell-averaged values.

---

## Summary of Key Physical Insights

### Why Relativistic Hydrodynamics is Fundamentally Different

1. **Lorentz factor coupling** **[Pons §1]**: All velocity components are coupled through $\gamma = (1 - v^2)^{-1/2}$ where $v^2$ includes all components. This is a manifestation of the geometric structure of Minkowski spacetime.

2. **Tangential velocity effects** **[Pons §3, §4]**: Unlike Newtonian hydrodynamics, tangential velocities affect the solution through the Riemann invariant $h \gamma v^t = \text{const}$. This arises from the relativistic momentum definition.

3. **Enthalpy coupling** **[Pons, discussion after Eq. 3.9]**: The specific enthalpy $h$ couples with velocities, important even in slow but thermodynamically relativistic flows ($h > 1$). This reflects the equivalence of mass and energy.

4. **No universal rest frame for Riemann problems** **[Pons §1]**: Unlike isolated shocks or rarefactions, one cannot transform away all tangential velocities in a Riemann problem. The solution intrinsically depends on the full velocity field.

5. **Complex eigenstructure** **[Pons Eq. 3.6]**: Characteristic speeds depend on both normal and tangential velocities, reflecting the hyperbolic structure of relativistic conservation laws.

### Why SPH Needs Riemann Solvers

**[SRGSPH §1, §2]**

1. **Shock capturing**: Traditional SPH with artificial viscosity requires tuning parameters for each problem. Riemann solvers automatically introduce the correct amount of dissipation.

2. **Conservation**: The anti-symmetric form of SRGSPH equations (Eqs. 371, 373) ensures exact conservation of total momentum and energy.

3. **Accuracy**: Convolution integrals (Eqs. 148, 252) provide $\mathcal{O}(h^2)$ accuracy compared to $\mathcal{O}(h)$ for standard SPH.

4. **Multidimensional capability**: The exact Riemann solver with tangential velocities enables genuinely multidimensional simulations without operator splitting.

---

## Implementation Notes

### Riemann Solver Implementation

**[Pons §5, §6]**

1. **Input**: States $L$ and $R$ with $(P, \rho, v^x, v^y, v^z)$
2. **Output**: Interface state $(P^*, v^{x*}, v^{y*}, v^{z*})$
3. **Algorithm**:
   - Compute sound speeds $c_{s,L}$ and $c_{s,R}$ (Eq. 2.9)
   - Set up functions $v_{L*}^x(P)$ and $v_{R*}^x(P)$ (§3 for rarefactions, §4 for shocks)
   - Solve $v_{L*}^x(P) - v_{R*}^x(P) = 0$ for $P^*$ (Eq. 1.2)
   - Evaluate interface velocity $v^{x*} = v_{L*}^x(P^*)$
   - Determine which intermediate state ($L_*$ or $R_*$) is at interface
   - Compute tangential velocities from appropriate invariant
   - Return interface state

### Computational Considerations

**[Pons §6]**

- **Rarefaction ODE**: Integrate using adaptive RK4 or similar (Eq. 3.10)
- **Root finding**: Brent's method recommended for $P_*$ (robust and fast)
- **Quartic solver** **[SRGSPH §2.8]**: Newton iteration or Ferrari's method for Lorentz factor
- **Smoothing length** **[SRGSPH, after Eq. 239]**: Iterate until $|h^{n+1} - h^n|/h^n < 10^{-6}$

**Efficiency** **[Pons §6]**: The computational cost of the exact Riemann solver with tangential velocities is comparable to the purely normal flow version, making it practical for multidimensional simulations.

---

## Symbol Translation Table

For reference when reading the original papers:

| Quantity | SRGSPH (this doc) | Pons et al. |
|----------|-------------------|-------------|
| Lorentz factor | $\gamma$ | $W$ |
| Shock Lorentz factor | $\gamma_s$ | $W_s$ |
| Pressure | $P$ | $p$ |
| Adiabatic index | $\gamma_c$ | $\gamma$ |
| Specific enthalpy | $h$ | $h$ |
| Rest-mass density | $\rho$ | $\rho$ |
| Kernel function | $W(\mathbf{x}, h)$ | (not in Pons) |

---

## References

1. **[SRGSPH]**: Kitajima, K., Inutsuka, S., & Seno, I. (2025). Special Relativistic Smoothed Particle Hydrodynamics Based on Riemann Solver. Journal of Computational Physics (accepted).

2. **[Pons]**: Pons, J.A., Martí, J.M., & Müller, E. (2000). The exact solution of the Riemann problem with non-zero tangential velocities in relativistic hydrodynamics. Journal of Fluid Mechanics, 422, 125-139.

3. Inutsuka, S. (2002). Reformulation of Smoothed Particle Hydrodynamics with Riemann Solver. Journal of Computational Physics, 179, 238-267.

4. Taub, A.H. (1948). Relativistic Rankine-Hugoniot Relations. Physical Review, 74, 328.

5. Taub, A.H. (1978). Relativistic Fluid Mechanics. Annual Review of Fluid Mechanics, 10, 32.

