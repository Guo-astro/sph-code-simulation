# GSPH with Self-Gravity: Complete Lagrangian Derivation

## The Correct Approach

You are absolutely right — self-gravity **must** be included in the Lagrangian from the beginning, and then we perform the variation to obtain the correct equations. This is the only way to guarantee:

1. **Variational consistency**
2. **Exact conservation laws**
3. **Proper coupling between pressure and gravity when h = h(ρ)**

Let me derive this correctly from first principles.

---

## 1. The Complete Lagrangian

### 1.1 Components

The total Lagrangian for a self-gravitating fluid is:

$$L = T - U_{\text{internal}} - U_{\text{grav}}$$

where:
- $T$ = kinetic energy
- $U_{\text{internal}}$ = internal (thermal) energy
- $U_{\text{grav}}$ = gravitational potential energy

### 1.2 GSPH Discretization

Following Inutsuka (2002), we discretize each term:

**Kinetic Energy:**
$$T = \sum_i \frac{1}{2} m_i \dot{\mathbf{x}}_i^2$$

**Internal Energy (GSPH form with convolution):**
$$U_{\text{internal}} = \sum_i m_i \int u(\mathbf{x}) W(\mathbf{x} - \mathbf{x}_i, h[\mathbf{x}]) \, d\mathbf{x}$$

Or in the simpler standard form:
$$U_{\text{internal}} = \sum_i m_i u_i(\rho_i, s_i)$$

**Gravitational Potential Energy:**
$$U_{\text{grav}} = -\frac{1}{2} \sum_i \sum_j G m_i m_j \phi(r_{ij}, h_{ij})$$

where $\phi(r, h)$ is the softened gravitational potential kernel, and $h_{ij}$ can be the average $(h_i + h_j)/2$ or some other symmetric combination.

### 1.3 The Complete Lagrangian

$$\boxed{L = \sum_i \frac{1}{2} m_i \dot{\mathbf{x}}_i^2 - \sum_i m_i u_i(\rho_i) + \frac{1}{2} \sum_i \sum_j G m_i m_j \phi(r_{ij}, h_{ij})}$$

**Critical point:** When $h = h(\rho)$, both $u_i(\rho_i)$ and $\phi(r_{ij}, h_{ij})$ depend implicitly on all particle positions through the density!

---

## 2. The Variation with Variable h

### 2.1 The Density Constraint

The density is defined as:
$$\rho_i = \sum_j m_j W(r_{ij}, h_i)$$

with the constraint:
$$h_i = \eta \left( \frac{m_i}{\rho_i} \right)^{1/d}$$

This means:
$$\frac{\partial h_i}{\partial \rho_i} = -\frac{h_i}{d \rho_i}$$

### 2.2 Variation of Density

When we vary particle position $\mathbf{x}_k$, the density at particle $i$ changes:

$$\frac{\partial \rho_i}{\partial \mathbf{x}_k} = \sum_j m_j \frac{\partial W(r_{ij}, h_i)}{\partial \mathbf{x}_k}$$

This involves two contributions:
1. Direct: change in $r_{ij}$ when $k = i$ or $k = j$
2. Indirect: change in $h_i$ through $\rho_i$

$$\frac{\partial W(r_{ij}, h_i)}{\partial \mathbf{x}_k} = \frac{\partial W}{\partial r_{ij}} \frac{\partial r_{ij}}{\partial \mathbf{x}_k} + \frac{\partial W}{\partial h_i} \frac{\partial h_i}{\partial \rho_i} \frac{\partial \rho_i}{\partial \mathbf{x}_k}$$

### 2.3 Solving for the Density Derivative

Collecting terms with $\partial \rho_i / \partial \mathbf{x}_k$:

$$\frac{\partial \rho_i}{\partial \mathbf{x}_k} = \sum_j m_j \frac{\partial W}{\partial r_{ij}} \frac{\partial r_{ij}}{\partial \mathbf{x}_k} + \left( \sum_j m_j \frac{\partial W}{\partial h_i} \right) \frac{\partial h_i}{\partial \rho_i} \frac{\partial \rho_i}{\partial \mathbf{x}_k}$$

$$\frac{\partial \rho_i}{\partial \mathbf{x}_k} \left[ 1 - \frac{\partial h_i}{\partial \rho_i} \sum_j m_j \frac{\partial W_{ij}}{\partial h_i} \right] = \sum_j m_j \frac{\partial W}{\partial r_{ij}} \frac{\partial r_{ij}}{\partial \mathbf{x}_k}$$

Define:
$$\boxed{\Omega_i \equiv 1 - \frac{\partial h_i}{\partial \rho_i} \sum_j m_j \frac{\partial W(r_{ij}, h_i)}{\partial h_i}}$$

Then:
$$\frac{\partial \rho_i}{\partial \mathbf{x}_k} = \frac{1}{\Omega_i} \sum_j m_j \nabla_k W_{ij}$$

---

## 3. Variation of the Internal Energy

### 3.1 The Chain Rule

$$\frac{\partial U_{\text{internal}}}{\partial \mathbf{x}_k} = \sum_i m_i \frac{\partial u_i}{\partial \rho_i} \frac{\partial \rho_i}{\partial \mathbf{x}_k}$$

Using $P = \rho^2 \partial u / \partial \rho$ (first law of thermodynamics):

$$= \sum_i m_i \frac{P_i}{\rho_i^2} \frac{\partial \rho_i}{\partial \mathbf{x}_k}$$

### 3.2 Substituting the Ω Factor

$$= \sum_i m_i \frac{P_i}{\rho_i^2} \cdot \frac{1}{\Omega_i} \sum_j m_j \nabla_k W_{ij}$$

$$= \sum_i \frac{m_i P_i}{\rho_i^2 \Omega_i} \sum_j m_j \nabla_k W_{ij}$$

---

## 4. Variation of the Gravitational Energy

### 4.1 The Softened Potential

$$U_{\text{grav}} = -\frac{1}{2} \sum_i \sum_j G m_i m_j \phi(r_{ij}, h_{ij})$$

where $h_{ij} = (h_i + h_j)/2$ (or similar symmetric form).

### 4.2 The Variation

$$\frac{\partial U_{\text{grav}}}{\partial \mathbf{x}_k} = -\frac{1}{2} \sum_i \sum_j G m_i m_j \frac{\partial \phi(r_{ij}, h_{ij})}{\partial \mathbf{x}_k}$$

This has **two contributions**:

**Direct (position dependence):**
$$\frac{\partial \phi}{\partial r_{ij}} \frac{\partial r_{ij}}{\partial \mathbf{x}_k}$$

**Indirect (through h = h(ρ)):**
$$\frac{\partial \phi}{\partial h_{ij}} \frac{\partial h_{ij}}{\partial \mathbf{x}_k}$$

### 4.3 The Indirect Term

$$\frac{\partial h_{ij}}{\partial \mathbf{x}_k} = \frac{1}{2} \left( \frac{\partial h_i}{\partial \mathbf{x}_k} + \frac{\partial h_j}{\partial \mathbf{x}_k} \right)$$

$$= \frac{1}{2} \left( \frac{\partial h_i}{\partial \rho_i} \frac{\partial \rho_i}{\partial \mathbf{x}_k} + \frac{\partial h_j}{\partial \rho_j} \frac{\partial \rho_j}{\partial \mathbf{x}_k} \right)$$

This couples gravity to the density field through h!

### 4.4 The Complete Gravitational Variation

$$\frac{\partial U_{\text{grav}}}{\partial \mathbf{x}_k} = -\sum_j G m_k m_j \frac{\partial \phi}{\partial r} \hat{\mathbf{r}}_{kj} - \frac{1}{2} \sum_i \sum_j G m_i m_j \frac{\partial \phi}{\partial h_{ij}} \frac{\partial h_{ij}}{\partial \mathbf{x}_k}$$

The first term is the **standard gravitational force**.

The second term is an **additional correction** when h varies — this is the **gravitational analog of the Ω correction**!

---

## 5. The Euler-Lagrange Equation

### 5.1 The Full Equation of Motion

$$m_k \ddot{\mathbf{x}}_k = -\frac{\partial U_{\text{internal}}}{\partial \mathbf{x}_k} - \frac{\partial U_{\text{grav}}}{\partial \mathbf{x}_k}$$

### 5.2 Pressure Force (with Ω)

From Section 3:

$$\mathbf{F}_k^{\text{pressure}} = -\sum_j m_j \left( \frac{P_k}{\rho_k^2 \Omega_k} \nabla_k W_{kj}^k + \frac{P_j}{\rho_j^2 \Omega_j} \nabla_k W_{kj}^j \right)$$

### 5.3 Gravitational Force (with h-correction)

$$\mathbf{F}_k^{\text{grav}} = \sum_j G m_k m_j \frac{\partial \phi}{\partial r} \hat{\mathbf{r}}_{kj} + \text{(h-derivative terms)}$$

The h-derivative terms can be written as:

$$\mathbf{F}_k^{\text{grav},h} = \frac{1}{2} \sum_i \sum_j G m_i m_j \frac{\partial \phi}{\partial h_{ij}} \frac{\partial h_{ij}}{\partial \mathbf{x}_k}$$

### 5.4 The Complete Equation

$$\boxed{m_k \ddot{\mathbf{x}}_k = -\sum_j m_j \left( \frac{P_k}{\rho_k^2 \Omega_k} \nabla W_{kj}^k + \frac{P_j}{\rho_j^2 \Omega_j} \nabla W_{kj}^j \right) + \sum_j G m_k m_j g(r_{kj}, h_{kj}) \hat{\mathbf{r}}_{kj} + \mathbf{F}_k^{h\text{-grav}}}$$

---

## 6. The Gravitational h-Correction Term

### 6.1 Explicit Form

The gravitational h-correction is:

$$\mathbf{F}_k^{h\text{-grav}} = \frac{1}{2} \sum_i \sum_j G m_i m_j \frac{\partial \phi(r_{ij}, h_{ij})}{\partial h_{ij}} \cdot \frac{\partial h_{ij}}{\partial \mathbf{x}_k}$$

Using $\partial h_i / \partial \rho_i = -h_i / (d \rho_i)$ and $\partial \rho_i / \partial \mathbf{x}_k = (1/\Omega_i) \sum_l m_l \nabla_k W_{il}$:

$$\mathbf{F}_k^{h\text{-grav}} = -\frac{1}{4d} \sum_i \sum_j G m_i m_j \frac{\partial \phi}{\partial h_{ij}} \left( \frac{h_i}{\rho_i \Omega_i} \sum_l m_l \nabla_k W_{il} + \frac{h_j}{\rho_j \Omega_j} \sum_l m_l \nabla_k W_{jl} \right)$$

### 6.2 Simplified Form

This can be rewritten as an effective "gravitational pressure":

$$\mathbf{F}_k^{h\text{-grav}} = -\sum_j m_j \left( \frac{\zeta_k}{\Omega_k} \nabla_k W_{kj}^k + \frac{\zeta_j}{\Omega_j} \nabla_k W_{kj}^j \right)$$

where:
$$\zeta_i = -\frac{h_i}{2d\rho_i} \sum_j G m_j \frac{\partial \phi(r_{ij}, h_{ij})}{\partial h_{ij}}$$

This is the **Price & Monaghan (2007)** gravitational grad-h term!

---

## 7. Conservation Properties

### 7.1 Momentum Conservation

The complete equation of motion has the form:

$$m_k \ddot{\mathbf{x}}_k = \sum_j \mathbf{F}_{kj}$$

where $\mathbf{F}_{kj} = -\mathbf{F}_{jk}$ (antisymmetric).

Therefore:
$$\frac{d}{dt} \sum_k m_k \dot{\mathbf{x}}_k = \sum_k \sum_j \mathbf{F}_{kj} = 0$$

**Momentum is exactly conserved.**

### 7.2 Energy Conservation

Since the equations derive from a Lagrangian:

$$\frac{dE}{dt} = \frac{d}{dt} \left( T + U_{\text{internal}} + U_{\text{grav}} \right) = 0$$

**Energy is exactly conserved** (in the absence of dissipation).

### 7.3 Angular Momentum Conservation

Since $\mathbf{F}_{kj} \parallel \mathbf{r}_{kj}$, angular momentum is conserved:

$$\frac{d}{dt} \sum_k m_k \mathbf{x}_k \times \dot{\mathbf{x}}_k = 0$$

---

## 8. Why This Matters for Hydrostatic Equilibrium

### 8.1 The Equilibrium Condition

In hydrostatic equilibrium:
$$\nabla P = \rho \mathbf{g}$$

The SPH version must satisfy:
$$\sum_j m_j \frac{P}{\rho^2 \Omega} \nabla W = \rho \cdot \left( \sum_j G m_j g(r, h) \hat{\mathbf{r}} + \text{h-correction} \right)$$

### 8.2 Without Proper Lagrangian Derivation

If we naively add gravity to SPH without going through the Lagrangian:

1. **Pressure force** might have Ω correction
2. **Gravity** has no corresponding correction
3. **Energy is not conserved**
4. **Hydrostatic balance is inconsistent**

### 8.3 With Proper Lagrangian Derivation

Both pressure and gravity get their h-corrections from the same variational principle:

1. **Pressure:** $P/(\rho^2 \Omega)$ instead of $P/\rho^2$
2. **Gravity:** Additional $\zeta/\Omega$ term from $\partial \phi / \partial h$
3. **Energy is exactly conserved**
4. **Hydrostatic balance is variationally consistent**

---

## 9. Practical Implementation

### 9.1 The Gravitational Potential Derivative

For the Wendland C4 kernel potential $\phi(r, h)$:

$$\frac{\partial \phi}{\partial h} = \frac{\partial}{\partial h} \left[ \frac{\text{poly}(r/h)}{h} \right]$$

This must be computed alongside the force kernel $g(r, h)$.

### 9.2 The ζ Term

Define the gravitational h-correction coefficient:

$$\zeta_i = -\frac{h_i}{2d\rho_i} \sum_j G m_j \frac{\partial \phi(r_{ij}, h_{ij})}{\partial h_{ij}}$$

This is computed during the gravity tree walk.

### 9.3 The Complete Force

$$\mathbf{F}_k = \underbrace{-\sum_j m_j \frac{P_k + P_j}{\rho_k \rho_j} \frac{1}{\bar{\Omega}} \nabla W}_{\text{pressure with Ω}} + \underbrace{\sum_j G m_k m_j g(r, h) \hat{\mathbf{r}}}_{\text{gravity}} + \underbrace{-\sum_j m_j \frac{\zeta_k + \zeta_j}{\rho_k \rho_j} \frac{1}{\bar{\Omega}} \nabla W}_{\text{gravity h-correction}}$$

---

## 10. Connection to Your Code

### 10.1 Current Implementation

In `g_fluid_force.cpp`, you have:

```cpp
const real omega_i = p_i.gradh;
const real omega_j = p_j.gradh;

const vec_t f = dw_i * (p_j.mass * pstar * rho2_inv_i * omega_i)
              + dw_j * (p_j.mass * pstar * rho2_inv_j * omega_j);
```

This correctly includes the **pressure Ω correction**.

### 10.2 What's Missing?

The gravitational h-correction term $\zeta/\Omega$ should be added to `gravity_force.cpp` for full variational consistency:

```cpp
// Compute ∂φ/∂h for each pair
real dphi_dh = compute_dphi_dh(r, h_ij);

// Accumulate ζ_i
zeta_i += G * m_j * dphi_dh;

// Later, add to force:
// F_grav += zeta_i / (rho_i * Omega_i) * sum_j m_j * grad_W_ij
```

### 10.3 Does It Matter in Practice?

For **hydrostatic equilibrium**, the gravitational h-correction is typically smaller than the pressure Ω correction because:

1. Gravity softening is often fixed ($\partial \phi / \partial h = 0$)
2. Or the kernel is designed so $\partial \phi / \partial h$ is small

But for **full variational consistency**, it should be included.

---

## 11. Summary

### The Key Insight

You are correct: **self-gravity must be in the Lagrangian from the start**. The proper derivation gives:

1. **Pressure force** with Ω correction: $P/(\rho^2 \Omega)$
2. **Gravitational force** with h-correction: $\zeta/\Omega$ term
3. **Exact conservation** of momentum, energy, angular momentum
4. **Consistent hydrostatic equilibrium**

### The Ω Factor Appears Everywhere

Because everything derives from the same Lagrangian with $h = h(\rho)$:

$$\Omega = 1 - \frac{\partial h}{\partial \rho} \sum_j m_j \frac{\partial W}{\partial h}$$

This factor multiplies **both** the pressure term **and** the gravitational h-correction, ensuring variational consistency.

### Why It Works

The Lagrangian approach guarantees that:
- All forces derive from a single variational principle
- Conservation laws are automatic (Noether's theorem)
- Hydrostatic equilibrium is self-consistent
- There are no "mismatched" corrections between pressure and gravity

---

## References

1. Inutsuka, S. (2002). "Reformulation of Smoothed Particle Hydrodynamics with Riemann Solver." JCP, 179, 238-267.

2. Springel, V. & Hernquist, L. (2002). "Cosmological smoothed particle hydrodynamics simulations: the entropy equation." MNRAS, 333, 649.

3. Price, D.J. & Monaghan, J.J. (2007). "An energy-conserving formalism for adaptive gravitational force softening in SPH and N-body codes." MNRAS, 374, 1347.

4. Hopkins, P.F. (2015). "A new class of accurate, mesh-free hydrodynamic simulation methods." MNRAS, 450, 53.
