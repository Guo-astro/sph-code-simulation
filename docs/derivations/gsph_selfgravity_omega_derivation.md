# GSPH with Self-Gravity and Ω Correction: A First-Principles Derivation

## Abstract

This document derives from first principles why the grad-h (Ω) correction is **essential** for self-gravitating hydrostatic equilibrium but **negligible** for shock tube simulations. We start from the Lagrangian formulation of GSPH (Inutsuka 2002), extend it to include self-gravity, and demonstrate mathematically why the Ω modification emerges naturally from variational consistency when the smoothing length varies with density.

---

## 1. The GSPH Lagrangian Framework

### 1.1 Fundamental Definitions

Following Inutsuka (2002), we define the density field from particle positions:

$$\rho(\mathbf{x}) = \sum_j m_j W(\mathbf{x} - \mathbf{x}_j, h[\mathbf{x}])$$

This definition satisfies the **exact identity**:

$$1 = \sum_j \frac{m_j}{\rho(\mathbf{x})} W(\mathbf{x} - \mathbf{x}_j, h)$$

which is the key to GSPH's consistency (unlike standard SPH which only approximately satisfies this).

### 1.2 The GSPH Lagrangian (Constant h)

For a barotropic fluid with $P = P(\rho)$, the Lagrangian is:

$$L = \sum_i m_i \left[ \frac{1}{2} \dot{\mathbf{x}}_i^2 - \int u(\mathbf{x}) W(\mathbf{x} - \mathbf{x}_i, h) \, d\mathbf{x} \right]$$

where $u$ is the specific internal energy. The convolution $\int u \, W \, d\mathbf{x}$ is second-order accurate:

$$\langle u \rangle(\mathbf{x}_i) = u(\mathbf{x}_i) + \frac{h_{\text{eff}}^2}{4} \nabla^2 u + O(h^4)$$

### 1.3 Euler-Lagrange Equation

The Euler-Lagrange equation:

$$\frac{d}{dt} \frac{\partial L}{\partial \dot{\mathbf{x}}_i} - \frac{\partial L}{\partial \mathbf{x}_i} = 0$$

yields the GSPH equation of motion:

$$\ddot{\mathbf{x}}_i = -\sum_j m_j \int \frac{P(\mathbf{x})}{\rho^2(\mathbf{x})} \left[ \frac{\partial}{\partial \mathbf{x}_i} - \frac{\partial}{\partial \mathbf{x}_j} \right] W(\mathbf{x} - \mathbf{x}_i, h) W(\mathbf{x} - \mathbf{x}_j, h) \, d\mathbf{x}$$

The antisymmetric appearance of $i$ and $j$ guarantees **exact conservation** of linear and angular momentum.

---

## 2. Extension to Self-Gravity

### 2.1 Adding Gravitational Potential Energy

The gravitational potential energy is:

$$U_{\text{grav}} = -\frac{1}{2} \sum_i \sum_j G m_i m_j \phi(|\mathbf{x}_i - \mathbf{x}_j|, h_{ij})$$

where $\phi(r, h)$ is the softened potential kernel (e.g., Hernquist-Katz or Wendland C4).

The complete Lagrangian becomes:

$$\boxed{L = \sum_i m_i \left[ \frac{1}{2} \dot{\mathbf{x}}_i^2 - \langle u \rangle_i \right] + \frac{1}{2} \sum_i \sum_j G m_i m_j \phi(r_{ij}, h_{ij})}$$

### 2.2 Gravitational Force

The gravitational acceleration is:

$$\mathbf{g}_i = -\frac{\partial U_{\text{grav}}}{\partial \mathbf{x}_i} = -\sum_j G m_j g(r_{ij}, h_{ij}) \mathbf{r}_{ij}$$

where $g(r, h) = -\frac{1}{r} \frac{\partial \phi}{\partial r}$ is the force kernel.

### 2.3 Complete Equation of Motion

The full equation of motion with self-gravity is:

$$m_i \ddot{\mathbf{x}}_i = \mathbf{F}_i^{\text{pressure}} + \mathbf{F}_i^{\text{gravity}}$$

where:
- **Pressure force** (from GSPH Riemann solver):
  $$\mathbf{F}_i^{\text{pressure}} = -\sum_j m_j P^*_{ij} \left[ V_{ij}^2(h_i) \nabla_i W_{ij}^i + V_{ij}^2(h_j) \nabla_i W_{ij}^j \right]$$

- **Gravitational force**:
  $$\mathbf{F}_i^{\text{gravity}} = -\sum_j G m_i m_j g(r_{ij}, h_{ij}) \hat{\mathbf{r}}_{ij}$$

---

## 3. The Variable Smoothing Length Problem

### 3.1 Why h Must Vary

In self-gravitating systems, the density contrast can span orders of magnitude:
- Core of a star: $\rho \sim 10^3 \rho_{\text{edge}}$
- Accretion disk: $\rho$ varies by factors of $10^4$

To maintain resolution (constant neighbor count $N_{\text{ngb}}$), we require:

$$h_i = \eta \left( \frac{m_i}{\rho_i} \right)^{1/d}$$

where $d$ is the spatial dimension and $\eta \sim 1$.

### 3.2 The Implicit Density Equation

With variable $h = h(\rho)$, the density becomes an **implicit** function:

$$\rho_i = \sum_j m_j W(r_{ij}, h(\rho_i))$$

Both sides depend on $\rho_i$! This must be solved iteratively or treated carefully.

### 3.3 The Hidden Derivative

When we take variations of the Lagrangian with respect to particle positions, we must account for:

$$\frac{\partial \rho}{\partial \mathbf{x}_i} = \sum_k m_k \left[ \nabla_i W_{ik} + \frac{\partial W_{ik}}{\partial h} \frac{\partial h}{\partial \rho} \frac{\partial \rho}{\partial \mathbf{x}_i} \right]$$

This is a **self-referential equation**! Solving for $\partial \rho / \partial \mathbf{x}_i$:

$$\frac{\partial \rho}{\partial \mathbf{x}_i} = \frac{1}{\Omega} \sum_k m_k \nabla_i W_{ik}$$

where:

$$\boxed{\Omega \equiv 1 - \frac{\partial h}{\partial \rho} \sum_j m_j \frac{\partial W(r_{ij}, h)}{\partial h}}$$

---

## 4. The Ω-Corrected GSPH Equations

### 4.1 Variational Derivation

Starting from the Lagrangian with $h = h(\rho)$:

$$L = \sum_i m_i \left[ \frac{1}{2} \dot{\mathbf{x}}_i^2 - u(\rho_i, s_i) \right] + U_{\text{grav}}$$

The Euler-Lagrange equation gives:

$$m_i \ddot{\mathbf{x}}_i = -\frac{\partial}{\partial \mathbf{x}_i} \sum_j m_j u_j(\rho_j)$$

Using the chain rule with $u = u(\rho)$ and $P = \rho^2 \partial u / \partial \rho$:

$$= -\sum_j m_j \frac{\partial u_j}{\partial \rho_j} \frac{\partial \rho_j}{\partial \mathbf{x}_i}$$

$$= -\sum_j m_j \frac{P_j}{\rho_j^2} \cdot \frac{1}{\Omega_j} \sum_k m_k \nabla_i W_{jk}$$

### 4.2 The Corrected Momentum Equation

The Ω-corrected GSPH momentum equation becomes:

$$\boxed{\ddot{\mathbf{x}}_i = -\sum_j m_j \left( \frac{P_i^*}{\rho_i^2 \Omega_i} \nabla_i W_{ij}^i + \frac{P_j^*}{\rho_j^2 \Omega_j} \nabla_i W_{ij}^j \right) + \mathbf{g}_i}$$

where $P^*$ comes from the Riemann solver.

### 4.3 The Corrected Energy Equation

Similarly, the energy equation:

$$\dot{u}_i = \frac{P_i}{\rho_i^2 \Omega_i} \sum_j m_j (\mathbf{v}_i - \mathbf{v}_j^*) \cdot \nabla_i W_{ij}^i$$

---

## 5. Why Ω is Essential for Hydrostatic Equilibrium

### 5.1 The Hydrostatic Balance Condition

In equilibrium, the pressure gradient must exactly balance gravity:

$$\nabla P = \rho \mathbf{g}$$

In SPH, this becomes:

$$\sum_j m_j \frac{P}{\rho^2} \nabla W_{ij} = \rho_i \mathbf{g}_i$$

### 5.2 The Error Without Ω Correction

**Without** the Ω correction, the SPH pressure gradient is:

$$\left( \nabla P \right)_{\text{SPH}} = \sum_j m_j \frac{P_j}{\rho_j^2} \nabla W(r_{ij}, h_j)$$

But when $h = h(\rho)$, the kernel gradient picks up a spurious term:

$$\nabla W(r, h(\rho)) = \frac{\partial W}{\partial r} \hat{\mathbf{r}} + \frac{\partial W}{\partial h} \nabla h$$

The second term is **not physical** — it's an artifact of the resolution change, not real pressure!

### 5.3 Quantifying the Error

For $h \propto \rho^{-1/3}$ (3D):

$$\nabla h = -\frac{h}{3\rho} \nabla \rho$$

The spurious force is:

$$\mathbf{F}_{\text{spurious}} \sim \sum_j m_j \frac{P}{\rho^2} \frac{\partial W}{\partial h} \left( -\frac{h}{3\rho} \nabla \rho \right)$$

For a polytropic Lane-Emden sphere with $\rho \propto r^{-2}$ in the envelope:

$$|\mathbf{F}_{\text{spurious}}| \sim \frac{P h}{\rho R} \sim \frac{h}{R} |\mathbf{F}_{\text{pressure}}|$$

This is an **O(h/R) systematic error** in the force balance!

### 5.4 The Consequence: Core Collapse

In the core of a self-gravitating sphere:
- $\rho$ is maximum → $h$ is minimum
- Moving outward: $\partial h / \partial r > 0$
- The spurious term **reduces** the effective pressure gradient
- Gravity is unbalanced → **core collapses**

### 5.5 How Ω Fixes This

The Ω correction modifies the pressure term:

$$\frac{P}{\rho^2} \to \frac{P}{\rho^2 \Omega}$$

This **exactly compensates** for the spurious $\partial W / \partial h$ contribution.

**Proof:** The Ω factor is defined precisely such that:

$$\frac{1}{\Omega} \nabla W = \nabla W \Big|_{h=\text{const}} + \text{correction for } \nabla h$$

gives the correct physical gradient, as if $h$ were constant.

### 5.6 Explicit Calculation

For the Wendland C4 kernel in 3D with $h = \eta (m/\rho)^{1/3}$:

$$\Omega = 1 + \frac{h}{3\rho} \sum_j m_j \frac{\partial W}{\partial h}$$

Typically, $\Omega \approx 0.7$ to $1.3$ depending on local density gradient.

In the core: $\Omega < 1$ → pressure force **enhanced**
At the edge: $\Omega \approx 1$ → standard SPH

This asymmetric correction restores the exact hydrostatic balance:

$$\frac{P}{\rho^2 \Omega} \nabla W = \frac{P_{\text{true}}}{\rho^2} \nabla W \Big|_{\text{physical}}$$

---

## 6. Why Ω is Negligible for Shock Tubes

### 6.1 Uniform Density Regions

In a Sod shock tube, the initial conditions are piecewise uniform:

$$\rho = \begin{cases} \rho_L & x < 0 \\ \rho_R & x > 0 \end{cases}$$

Within each uniform region:
- $\nabla \rho = 0$ → $\nabla h = 0$
- $\frac{\partial h}{\partial \rho} \sum_j m_j \frac{\partial W}{\partial h} \approx 0$
- Therefore: $\Omega \approx 1$

### 6.2 No Equilibrium to Maintain

The Sod shock tube is a **transient** problem:
- Initial discontinuity → shock + rarefaction + contact
- No static force balance required
- Dynamics determined by wave propagation, not equilibrium

A 1% error in force magnitude slightly shifts wave speeds, but doesn't cause catastrophic failure.

### 6.3 Discontinuities Handled by Riemann Solver

At the shock and contact discontinuity:
- $h$ changes abruptly
- But the **Riemann solver** captures the physics correctly
- $P^*$ and $v^*$ are computed from the jump conditions
- The Ω correction is a smooth O(h) effect — in the noise

### 6.4 Error Tolerance

| Property | Self-Gravity Hydrostatic | Sod Shock Tube |
|----------|-------------------------|----------------|
| Force balance required | Yes, to 0.1% | No |
| Density profile | Smooth, continuous | Piecewise uniform |
| $\nabla h$ | Nonzero everywhere | ~0 in uniform regions |
| $\Omega - 1$ | O(0.1) to O(0.3) | O(0.01) |
| Error accumulation | Cumulative over time | Local, transient |
| Failure mode without Ω | Core collapse | Slight wave speed error |

---

## 7. Mathematical Summary

### 7.1 The Critical Difference

Define the dimensionless grad-h parameter:

$$\epsilon_{\Omega} \equiv |1 - \Omega| \sim \left| \frac{h}{\rho} \frac{\partial \rho}{\partial r} \right|$$

**Self-gravitating hydrostatic:**
$$\epsilon_{\Omega} \sim \frac{h}{R_{\text{system}}} \sim 0.1 \text{ to } 0.3$$

**Sod shock tube:**
$$\epsilon_{\Omega} \sim 0 \text{ (in uniform regions)}$$

### 7.2 Force Error Scaling

Without Ω correction, the relative force error is:

$$\frac{\delta F}{F} \sim \epsilon_{\Omega}$$

**Hydrostatic case:** 10-30% force error → equilibrium impossible → collapse

**Shock tube:** 0-1% force error → slightly shifted wave speeds → acceptable

### 7.3 Time Integration Effect

The position error after time $t$:

$$\delta x \sim \frac{\delta F}{m} t^2 \sim \epsilon_{\Omega} \frac{F}{m} t^2$$

**Hydrostatic (must hold for $t \to \infty$):**
Any nonzero $\epsilon_{\Omega}$ causes unbounded drift → collapse

**Shock tube (simulation time $t \sim R/c_s$):**
$$\delta x \sim \epsilon_{\Omega} \cdot h \sim 0 \text{ for uniform regions}$$

---

## 8. Implementation in Your Code

From `g_fluid_force.cpp`:

```cpp
// Grad-h correction (Springel & Hernquist 2002)
const real omega_i = p_i.gradh;  // Note: this is 1/Ω
const real omega_j = p_j.gradh;

// Standard GSPH force with pstar from well-balanced Riemann solver
const vec_t f = dw_i * (p_j.mass * pstar * rho2_inv_i * omega_i)
              + dw_j * (p_j.mass * pstar * rho2_inv_j * omega_j);
```

Here `omega_i = 1/Ω_i`, so multiplying by `omega_i` effectively divides the pressure by Ω, giving the correct grad-h compensated force.

---

## 9. Conclusion

The Ω (grad-h) correction is **mathematically required** when:

1. **$h$ varies with $\rho$** (adaptive resolution)
2. **Continuous density gradients exist** (smooth profiles)
3. **Force balance must be maintained** (equilibrium)

It is **negligible** when:

1. **Density is piecewise uniform** (shock tubes)
2. **h is approximately constant** locally
3. **Dynamics dominate** (no equilibrium required)

For self-gravitating systems like Lane-Emden polytropes, the Ω correction is **the difference between a stable hydrostatic sphere and catastrophic core collapse**. For Sod shock tubes, it's a small perturbation that can be safely ignored.

---

## References

1. Inutsuka, S. (2002). "Reformulation of Smoothed Particle Hydrodynamics with Riemann Solver." JCP, 179, 238-267.

2. Springel, V. & Hernquist, L. (2002). "Cosmological smoothed particle hydrodynamics simulations: the entropy equation." MNRAS, 333, 649.

3. Hopkins, P.F. (2013). "A general class of Lagrangian smoothed particle hydrodynamics methods and implications for fluid mixing problems." MNRAS, 428, 2840.

4. Price, D.J. & Monaghan, J.J. (2007). "An energy-conserving formalism for adaptive gravitational force softening in SPH and N-body codes." MNRAS, 374, 1347.

5. Hernquist, L. & Katz, N. (1989). "TREESPH: A unification of SPH with the hierarchical tree method." ApJS, 70, 419.

---

## Appendix: Computing Ω

For $h = \eta (m/\rho)^{1/d}$ where $d$ is dimension:

$$\frac{\partial h}{\partial \rho} = -\frac{h}{d \rho}$$

Therefore:

$$\Omega_i = 1 + \frac{h_i}{d \rho_i} \sum_j m_j \frac{\partial W(r_{ij}, h_i)}{\partial h_i}$$

For the Wendland C4 kernel $W(r, h) = \frac{\sigma}{h^d} w(q)$ where $q = r/h$:

$$\frac{\partial W}{\partial h} = -\frac{d}{h} W - \frac{q}{h} W'(q)$$

Substituting:

$$\Omega_i = 1 - \frac{1}{d} \sum_j \frac{m_j}{\rho_i} \left[ d \cdot W_{ij} + q_{ij} W'(q_{ij}) \right] \frac{1}{h_i}$$

This is computed during the density summation loop and stored as `p_i.gradh = 1/Ω_i` for efficiency.
