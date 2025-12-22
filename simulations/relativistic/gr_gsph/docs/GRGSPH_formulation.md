# GRGSPH: General Relativistic Godunov Smoothed Particle Hydrodynamics

## A First-Principles Derivation

**Based on:**
- Kitajima, Inutsuka & Seno (2025) - SRGSPH formulation
- Liptai & Price (2019) - GRSPH framework

---

## 1. Introduction

This document presents the derivation of **GRGSPH (General Relativistic Godunov SPH)**, which combines:
1. Kitajima's Riemann solver approach for accurate shock capturing
2. Liptai & Price's general relativistic hydrodynamics framework

The key innovation is replacing artificial viscosity with a Riemann solver while maintaining full general relativistic covariance.

---

## 2. Notation and Conventions

### 2.1 Index Conventions

| Index Type | Range | Description |
|------------|-------|-------------|
| Greek (μ, ν, ...) | 0,1,2,3 | Spacetime indices |
| Latin (i, j, k, ...) | 1,2,3 | Spatial indices |
| Subscripts a, b | - | Particle labels |

### 2.2 Fundamental Variables

| Symbol | Definition | Description |
|--------|------------|-------------|
| $n$ | - | Baryon number density (rest frame) |
| $N$ | $\Gamma n$ | Baryon number density (coordinate frame) |
| $N^*$ | $\sqrt{-g}\, n\, U^0$ | Conserved baryon number density |
| $\nu$ | - | Baryon number per SPH particle |
| $\rho$ | $m_b n$ | Rest mass density |
| $u$ | - | Specific internal energy |
| $P$ | - | Pressure |
| $H$ | $1 + u + P/(nc^2)$ | Enthalpy per baryon |
| $w$ | $1 + u + P/\rho$ | Specific enthalpy |
| $\mathbf{S}$ | $H\Gamma\mathbf{V}$ | Momentum per baryon |
| $e$ | see §3.3 | Energy per baryon |

### 2.3 Units and Signature

- Natural units: $c = G = 1$
- Metric signature: $(-,+,+,+)$
- Adiabatic index: $\gamma_c$ (typically 5/3)

---

## 3. Spacetime Structure

### 3.1 The 3+1 Decomposition (ADM Formalism)

The spacetime metric is decomposed as:

$$ds^2 = -\alpha^2 dt^2 + \gamma_{ij}(dx^i + \beta^i dt)(dx^j + \beta^j dt)$$

where:
- $\alpha$ = **lapse function** (relates proper time to coordinate time)
- $\beta^i$ = **shift vector** (relates spatial coordinates between slices)
- $\gamma_{ij}$ = **spatial 3-metric** (geometry of spatial hypersurfaces)

The metric determinant:
$$\sqrt{-g} = \alpha\sqrt{\gamma}$$

where $\gamma = \det(\gamma_{ij})$.

### 3.2 Four-Velocity

The four-velocity is:
$$U^\mu = \frac{dx^\mu}{d\tau}$$

The coordinate velocity:
$$v^\mu \equiv \frac{dx^\mu}{dt} = \frac{U^\mu}{U^0}$$

where $v^0 = 1$ by definition.

### 3.3 Eulerian Observer Frame

The velocity measured by a local Eulerian observer (zero angular momentum observer, ZAMO):
$$V^i = \frac{v^i + \beta^i}{\alpha}$$

The **generalized Lorentz factor**:
$$\Gamma = \frac{1}{\sqrt{1 - V^i V_i}} = \frac{1}{\sqrt{1 - \gamma_{ij}V^i V^j}}$$

Relation to $U^0$:
$$U^0 = \frac{\Gamma}{\alpha}$$

---

## 4. Conserved Variables in GR

### 4.1 Conserved Baryon Number Density

$$N^* \equiv \sqrt{-g}\, n\, U^0 = \sqrt{\gamma}\, \Gamma\, n$$

This is the conserved quantity that replaces $\rho^*$ in mass-based formulations.

### 4.2 Momentum Per Baryon (Covariant)

$$S_i = H \Gamma V_i = w \Gamma V_i$$

where the spatial index is lowered with $\gamma_{ij}$:
$$S_i = \gamma_{ij} S^j$$

### 4.3 Energy Per Baryon

$$e = S_i v^i + \frac{\alpha(1+u)}{\Gamma}$$

Alternative form:
$$e = w\Gamma(\alpha - V^i\beta_i) - \frac{\alpha P}{\Gamma n}$$

---

## 5. Fundamental Equations of GR Hydrodynamics

### 5.1 Stress-Energy Tensor

For a perfect fluid:
$$T^{\mu\nu} = (\rho + \rho u + P)U^\mu U^\nu + P g^{\mu\nu}$$
$$= \rho w U^\mu U^\nu + P g^{\mu\nu}$$
$$= n H (U^0)^2 v^\mu v^\nu + P g^{\mu\nu}$$

### 5.2 Conservation Laws

From $\nabla_\mu T^{\mu\nu} = 0$ and baryon conservation $\nabla_\mu(n U^\mu) = 0$:

**Continuity:**
$$\frac{dN^*}{dt} = -N^* \frac{\partial v^i}{\partial x^i}$$

**Momentum:**
$$\frac{dS_i}{dt} = -\frac{1}{N^*}\frac{\partial(\sqrt{-g}P)}{\partial x^i} + f_i$$

**Energy:**
$$\frac{de}{dt} = -\frac{1}{N^*}\frac{\partial(\sqrt{-g}P v^i)}{\partial x^i} + \Lambda$$

### 5.3 Metric Source Terms

The **gravitational force** (momentum source):
$$f_i = \frac{\sqrt{-g}}{2N^*} T^{\mu\nu}\frac{\partial g_{\mu\nu}}{\partial x^i}$$

The **gravitational work** (energy source):
$$\Lambda = -\frac{\sqrt{-g}}{2N^*} T^{\mu\nu}\frac{\partial g_{\mu\nu}}{\partial t}$$

Expanded form:
$$\frac{\sqrt{-g}T^{\mu\nu}}{2N^*} = \frac{1}{2}\left[H U^0 v^\mu v^\nu + \frac{P g^{\mu\nu}}{n U^0}\right]$$

---

## 6. Volume-Based SPH Formulation

### 6.1 Particle Volume (Kitajima Approach)

Instead of focusing on density, we define the **particle volume field**:
$$V_p(\mathbf{x}) = \left[\sum_j W(\mathbf{x}-\mathbf{x}_j, h(\mathbf{x}))\right]^{-1}$$

### 6.2 Variable Smoothing Length

The smoothing length is defined self-consistently:
$$h(\mathbf{x}) = \eta \left[V_p^*(\mathbf{x})\right]^{1/d}$$

where $d$ is the number of spatial dimensions and:
$$V_p^*(\mathbf{x}) = \left[\sum_j W(\mathbf{x}-\mathbf{x}_j, C_{\rm smooth} h(\mathbf{x}))\right]^{-1}$$

**Recommended parameters:**
- $\eta = 1.0$
- $C_{\rm smooth} = 2.0$

The $C_{\rm smooth} > 1$ factor ensures the smoothing length varies slowly, allowing us to neglect $\nabla h$ terms.

### 6.3 Number Density from Volume

$$N^*(\mathbf{x}) = \nu(\mathbf{x}) V_p^{-1}(\mathbf{x})$$

where $\nu(\mathbf{x})$ is the baryon number field (interpolated from particles).

### 6.4 Gaussian Kernel

$$W(\mathbf{x}, h) = \left[\frac{1}{h\sqrt{\pi}}\right]^d \exp\left[-\frac{\mathbf{x}^2}{h^2}\right]$$

The kernel gradient:
$$\nabla W(\mathbf{x}, h) = -\frac{2\mathbf{x}}{h^2} W(\mathbf{x}, h)$$

### 6.5 Convolution Definition

Physical quantities at particle $i$ are defined by convolution:
$$\langle N_i \rangle \equiv N(\mathbf{x}_i)$$
$$\langle \nu_i f_i \rangle \equiv \int \nu(\mathbf{x}) f(\mathbf{x}) W(\mathbf{x}-\mathbf{x}_i, h(\mathbf{x})) d^3x$$

---

## 7. Derivation of GRGSPH Equations

### 7.1 Momentum Equation Derivation

Starting from the convolution of the momentum equation:
$$\langle \nu_i \dot{S}_k^i \rangle = \int \nu(\mathbf{x}) \frac{dS_k(\mathbf{x})}{dt} W(\mathbf{x}-\mathbf{x}_i, h(\mathbf{x})) d^3x$$

Substituting the GR momentum equation:
$$= -\int V_p(\mathbf{x}) \sqrt{-g}(\mathbf{x}) \frac{\partial P}{\partial x^k} W_i d^3x + \langle \nu_i f_k^i \rangle$$

where $W_i = W(\mathbf{x}-\mathbf{x}_i, h(\mathbf{x}))$.

Integration by parts:
$$= +\int \sqrt{-g} P \frac{\partial}{\partial x^k}[V_p W_i] d^3x + \langle \nu_i f_k^i \rangle$$

Using $\partial_k V_p = -V_p^2 \sum_j \partial_k W_j$:
$$= +\sum_j \int \sqrt{-g} P V_p^2 \left[W_j \nabla_k W_i - W_i \nabla_k W_j\right] d^3x + \langle \nu_i f_k^i \rangle$$

### 7.2 Evaluation of Convolution Integrals

The key integral (following Kitajima):
$$\int f(\mathbf{x}) V_p^2(\mathbf{x}) [\nabla_i - \nabla_j] W_i W_j \, d^3x$$

Using the Gaussian kernel product identity:
$$W(\mathbf{x}-\mathbf{x}_i, h) W(\mathbf{x}-\mathbf{x}_j, h) \propto \exp\left[-\frac{2(\mathbf{x} - \frac{\mathbf{x}_i+\mathbf{x}_j}{2})^2 + \frac{(\mathbf{x}_i-\mathbf{x}_j)^2}{2}}{h^2}\right]$$

This evaluates to:
$$\approx f_{ij} V^2_{ij} \left[\nabla_i W(\mathbf{x}_i-\mathbf{x}_j, \sqrt{2}h_i) - \nabla_j W(\mathbf{x}_i-\mathbf{x}_j, \sqrt{2}h_j)\right]$$

where:
$$V^2_{ij} = \frac{1}{2}\left(V^2_{ij}(h_i) + V^2_{ij}(h_j)\right)$$

---

## 8. Final GRGSPH Equations

### 8.1 Equation of Motion

$$\boxed{\langle \nu_i \dot{S}_k^i \rangle = -\sum_j P^*_{ij} \overline{\sqrt{-g}}_{ij} V^2_{ij,\text{interp}} \left[\nabla_k W_{ij}(\sqrt{2}h_i) - \nabla_k W_{ij}(\sqrt{2}h_j)\right] + \nu_i f_k^i}$$

### 8.2 Energy Equation

$$\boxed{\langle \nu_i \dot{e}_i \rangle = -\sum_j P^*_{ij} v^{*k}_{ij} \overline{\sqrt{-g}}_{ij} V^2_{ij,\text{interp}} \left[\nabla_k W_{ij}(\sqrt{2}h_i) - \nabla_k W_{ij}(\sqrt{2}h_j)\right] + \nu_i \Lambda_i}$$

### 8.3 Definitions

| Term | Definition |
|------|------------|
| $P^*_{ij}$ | Pressure from Riemann solver |
| $v^*_{ij}$ | Velocity from Riemann solver |
| $V^2_{ij}$ | $\frac{1}{2}(V^2_{ij}(h_i) + V^2_{ij}(h_j))$ |
| $\overline{\sqrt{-g}}_{ij}$ | $\frac{1}{2}(\sqrt{-g_i} + \sqrt{-g_j})$ |
| $W_{ij}(h)$ | $W(\mathbf{x}_i - \mathbf{x}_j, h)$ |
| $\nabla_k W_{ij}(h)$ | $\partial W_{ij}/\partial x^k_i$ |

### 8.4 Conservation Properties

The equations satisfy antisymmetry under exchange of $i \leftrightarrow j$:
- $P^*_{ij} = P^*_{ji}$
- $v^*_{ij} = v^*_{ji}$
- $V^2_{ij} = V^2_{ji}$
- $\nabla_i W_{ij}(\sqrt{2}h_i) - \nabla_j W_{ij}(\sqrt{2}h_j)$ is antisymmetric

This guarantees:
- **Total momentum conservation** (in absence of external forces)
- **Total energy conservation** (for stationary metrics)

---

## 9. Riemann Problem in GR

### 9.1 Local Frame Projection

For each particle pair $(i,j)$, project to the line-of-sight direction:

**Unit vector along line of sight:**
$$\mathbf{e}_{ij} = \frac{\mathbf{x}_i - \mathbf{x}_j}{|\mathbf{x}_i - \mathbf{x}_j|}$$

**Line-of-sight velocities** (projection onto $\mathbf{e}_{ij}$):
$$V^*_a = \gamma_{kl}\, e^k_{ab}\, V^l_a \quad \text{(for particle } a \text{)}$$

Or equivalently using lowered index on $e$:
$$V^*_a = e_{ab,k}\, V^k_a$$

where $e_{ab,k} = \gamma_{kl}\, e^l_{ab}$ and $V^k_a$ is the Eulerian velocity of particle $a$.

For the pair $(i, j)$:
$$V^*_i = \gamma_{kl}\, e^k_{ij}\, V^l_i, \quad V^*_j = \gamma_{kl}\, e^k_{ij}\, V^l_j$$

**Line-of-sight Lorentz factors:**
$$\Gamma^*_a = \frac{1}{\sqrt{1 - (V^*_a)^2}}$$

### 9.2 Justification: Why Use SR Riemann Solver in GR?

**The Question:** Kitajima's Riemann solver is derived for special relativity (flat spacetime). How can we use it in curved spacetime?

**The Answer:** This is justified by the **equivalence principle** and **local flatness**, and is the same approach used by all major GR hydrodynamics codes (HARM, WhiskyTHC, IllinoisGRMHD, etc.).

#### 9.2.1 Theoretical Basis

**1. Local Flatness Theorem:**
At any point $p$ in spacetime, we can choose local inertial coordinates (freely falling frame) where:
- The metric becomes Minkowski: $g_{\mu\nu}(p) = \eta_{\mu\nu}$
- Christoffel symbols vanish: $\Gamma^\alpha_{\mu\nu}(p) = 0$
- First derivatives of metric vanish: $\partial_\lambda g_{\mu\nu}(p) = 0$

**2. GR Equations in Conservative Form:**
The GR hydrodynamic equations can be written as:
$$\partial_t \mathbf{U} + \partial_i \mathbf{F}^i = \mathbf{S}$$

where:
- $\mathbf{U} = \sqrt{\gamma}(D, S_j, \tau)^T$ — conserved variables
- $\mathbf{F}^i$ — fluxes (same structure as SR when using Eulerian variables)
- $\mathbf{S}$ — source terms (contain ALL metric derivatives)

**3. Separation of Physics:**
- The **flux terms** $\mathbf{F}^i$ describe fluid dynamics (pressure waves, shocks)
- The **source terms** $\mathbf{S}$ describe gravitational effects (metric gradients)

The Riemann problem addresses ONLY the flux terms. Gravitational effects enter separately through $f_i$ and $\Lambda$.

#### 9.2.2 The Local Frame Approach

For each particle pair $(i, j)$:

**Step 1: Identify the local region**
The interaction occurs over distance $|\mathbf{x}_i - \mathbf{x}_j| \sim h$

**Step 2: Use Eulerian observer frame**
The velocities $V^i = (v^i + \beta^i)/\alpha$ are already measured by the local ZAMO (zero angular momentum observer). This is the natural "local inertial frame."

**Step 3: Check validity condition**
The approximation is valid when:
$$h \ll L_{\text{curvature}} \sim |\mathcal{R}_{\mu\nu\rho\sigma}|^{-1/2}$$

For Schwarzschild: $L_{\text{curvature}} \sim r^{3/2}/\sqrt{M}$
- At $r = 10M$: $L \sim 30M$ → need $h \ll 30M$ ✓
- At $r = 3M$: $L \sim 5M$ → need $h \ll 5M$ ✓
- At horizon $r = 2M$: $L \sim 3M$ → need $h \ll 3M$ ✓ (if resolved)

**Step 4: Solve SR Riemann problem**
In this local frame, use the special relativistic Riemann solver.

**Step 5: GR effects added separately**
The metric gradient terms $f_i$, $\Lambda$ are computed and added to the RHS.

#### 9.2.3 Why This Works

| Aspect | How it's Handled |
|--------|-----------------|
| Local fluid dynamics | SR Riemann solver (flux terms) |
| Gravitational redshift | Through lapse $\alpha$ in $V^i$, time integration |
| Frame dragging | Through shift $\beta^i$ in $V^i$ |
| Spatial curvature | Through 3-metric $\gamma_{ij}$ in projections |
| Gravitational forces | Source terms $f_i$, $\Lambda$ |

#### 9.2.4 Precedent in Literature

This approach is standard in computational GR:

1. **Valencia Formulation** (Banyuls et al. 1997): Write GR hydro in flux-conservative form, use SR-like Riemann solvers

2. **HARM Code** (Gammie et al. 2003): Uses HLL solver derived from SR characteristics

3. **WhiskyTHC** (Radice et al. 2012): Uses HLLC solver with SR wave speeds

4. **Liptai & Price (2019)**: Their artificial viscosity is also based on SR signal speeds

**Quote from Martí & Müller (2015) Living Review:**
> "The local Riemann problem at cell interfaces is solved using approximate solvers derived for special relativistic hydrodynamics, which is justified by the equivalence principle."

---

### 9.3 1D Relativistic Riemann Problem

Solve the 1D **special relativistic** Riemann problem with:

**Left state (particle i):**
$$(n_L, P_L, V^*_L) = (n_i, P_i, V^*_i)$$

**Right state (particle j):**
$$(n_R, P_R, V^*_R) = (n_j, P_j, V^*_j)$$

**Output:** Interface pressure $P^*_{ij}$ and velocity $v^*_{ij}$

### 9.4 Riemann Solver Options

1. **Exact solver** (Pons et al. 2000) - Most accurate, more expensive
2. **HLLC solver** - Good balance of accuracy and speed
3. **HLL solver** - Simplest, more diffusive

### 9.5 MUSCL Reconstruction (2nd Order)

For higher accuracy, extrapolate primitives to the interface:

$$f_i^L = f_i + \frac{1}{2}\left.\frac{\partial f}{\partial s}\right|_i \cdot |\mathbf{x}_i - \mathbf{x}_j|$$
$$f_j^R = f_j - \frac{1}{2}\left.\frac{\partial f}{\partial s}\right|_j \cdot |\mathbf{x}_i - \mathbf{x}_j|$$

where $f \in \{n, P, V^*\}$ and $s$ is the coordinate along $\mathbf{e}_{ij}$.

### 9.6 Monotonicity Constraints (Limiter)

Switch to 1st order (set gradients to zero) when:

**Shock detection:**
$$C_{\rm shock}\, \mathbf{e}_{ij} \cdot (\mathbf{v}_i - \mathbf{v}_j) > \min(c_{s,i}, c_{s,j})$$

**Contact discontinuity detection:**
$$\left|\log_{10}\left(\frac{P_i}{P_j}\right)\right| > C_{\rm c.d.}$$

**Recommended parameters:**
- $C_{\rm shock} = 3$
- $C_{\rm c.d.} = 1$

---

## 10. Sound Speed in GR

The relativistic sound speed:
$$c_s = \sqrt{\frac{\gamma_c P}{n H}} = \sqrt{\frac{(\gamma_c - 1)(H-1)}{H}}$$

For ideal gas with $P = (\gamma_c - 1)nu$:
$$c_s = \sqrt{\frac{\gamma_c(\gamma_c - 1)u}{1 + \gamma_c u}}$$

---

## 11. Primitive Variable Recovery

### 11.1 The Inversion Problem

Given conserved variables $(N^*, S_i, e)$ and the metric, recover primitives $(n, u, v^i, P)$.

### 11.2 Express Everything in Terms of Enthalpy

**Lorentz factor:**
$$\Gamma(w) = \sqrt{1 + \frac{S^i S_i}{w^2}}$$

where $S^i = \gamma^{ij} S_j$.

**Number density:**
$$n(w) = \frac{N^*}{\sqrt{\gamma}\,\Gamma(w)}$$

**Pressure:**
$$P(w) = \frac{N^*}{\alpha\sqrt{\gamma}}\left[w\Gamma(w)\alpha - e - S_i\beta^i\right]$$

**Velocity:**
$$v_i(w) = \frac{\alpha S_i}{w\Gamma(w)} - \beta_i$$
$$v^i(w) = \gamma^{ij} v_j(w)$$

### 11.3 Newton-Raphson Iteration

For ideal gas EOS:
$$w = 1 + \frac{\gamma_c}{\gamma_c - 1}\frac{P}{\rho} = 1 + \frac{\gamma_c}{\gamma_c - 1}\frac{P}{m_b n}$$

Define residual:
$$f(w) = w(n(w), P(w)) - w$$

Iterate:
$$w_{n+1} = w_n - \frac{f(w_n)}{f'(w_n)}$$

until $|w_{n+1} - w_n|/w_{n+1} < \epsilon_w$ (typically $\epsilon_w = 10^{-12}$).

### 11.4 Derivative for Newton-Raphson

$$f'(w) = \frac{\gamma_c}{\gamma_c - 1}\left(1 - \frac{S^i S_i P}{w^3 n \Gamma^2}\right) - 1$$

---

## 12. Time Integration

### 12.1 CFL Condition

$$\Delta t = C_{\rm CFL} \min_i \left[\frac{h_i}{\alpha_i(c_{s,i} + |V_i|)}\right]$$

**Recommended:** $C_{\rm CFL} = 0.3$

Note: The factor of $\alpha_i$ accounts for time dilation near the black hole.

### 12.2 Simple Euler Method

$$\langle \nu_i S_k^i \rangle^{n+1} = \langle \nu_i S_k^i \rangle^n + \langle \nu_i \dot{S}_k^i \rangle \Delta t$$
$$\langle \nu_i e_i \rangle^{n+1} = \langle \nu_i e_i \rangle^n + \langle \nu_i \dot{e}_i \rangle \Delta t$$
$$x^{i,n+1}_a = x^{i,n}_a + v^i_a \Delta t$$

### 12.3 Recommended: Leapfrog/Verlet

For better energy conservation, use symplectic integrators:

**Kick-Drift-Kick:**
```
S^{n+1/2} = S^n + (dS/dt)^n · Δt/2
x^{n+1} = x^n + v(S^{n+1/2}) · Δt
S^{n+1} = S^{n+1/2} + (dS/dt)^{n+1} · Δt/2
```

### 12.4 Substeps for External Forces

Following Liptai & Price, use substeps for the metric source terms:
- SPH forces: timestep $\Delta t_{\rm sph}$
- Metric forces: timestep $\Delta t_{\rm ext} = \Delta t_{\rm sph}/m$

where $m = \lceil \Delta t_{\rm sph}/\Delta t_{\rm ext} \rceil$.

---

## 13. Specific Metrics

### 13.1 Minkowski (Special Relativity)

$$\alpha = 1, \quad \beta^i = 0, \quad \gamma_{ij} = \delta_{ij}$$
$$\sqrt{-g} = 1$$

Source terms vanish: $f_i = 0$, $\Lambda = 0$.

### 13.2 Schwarzschild (Non-rotating Black Hole)

In Cartesian coordinates with $r = \sqrt{x^2 + y^2 + z^2}$:

$$\alpha = \sqrt{1 - \frac{2M}{r}}$$
$$\beta^i = 0$$
$$\gamma_{ij} = \delta_{ij} + \frac{2M}{r(r-2M)}\frac{x_i x_j}{r^2}$$

### 13.3 Kerr (Rotating Black Hole)

In Cartesian Boyer-Lindquist coordinates (see Liptai & Price Appendix for full expressions):

$$r = \sqrt{\frac{R^2 - a^2 + \sqrt{(R^2-a^2)^2 + 4a^2z^2}}{2}}$$

where $R^2 = x^2 + y^2 + z^2$ and $a$ is the spin parameter.

---

## 14. Algorithm Summary

```
GRGSPH Algorithm
================

INITIALIZATION:
1. Set particle positions, baryon numbers ν_i
2. Initialize primitives (n, P, v)
3. Compute initial conserved variables (N*, S, e)

MAIN LOOP:
1. SMOOTHING LENGTH:
   - Iterate h(x) = η[V_p*(x)]^(1/d) until converged

2. DENSITY/VOLUME:
   - Compute V_p,i for each particle
   - Compute N*_i = ν_i / V_p,i

3. PRIMITIVE RECOVERY:
   - For each particle: (N*, S, e) → (n, u, P, v)
   - Use Newton-Raphson on enthalpy

4. RIEMANN SOLVER (for each pair i,j):
   a. Compute line-of-sight velocities V*_i, V*_j
   b. Apply MUSCL reconstruction (if 2nd order)
   c. Check limiter conditions
   d. Solve Riemann problem → P*_ij, v*_ij

5. FORCE CALCULATION:
   - Compute kernel gradients ∇W_ij(√2 h)
   - Compute V²_ij,interp
   - Accumulate SPH forces
   - Add GR source terms f_i, Λ

6. TIME INTEGRATION:
   - Update S, e, x
   - Apply boundary conditions

7. REPEAT from step 1
```

---

## 15. Comparison: GRSPH vs GRGSPH

| Aspect | GRSPH (Liptai & Price) | GRGSPH (This Work) |
|--------|------------------------|-------------------|
| **Shock Handling** | Artificial viscosity ($\alpha_{\rm AV}$) | Riemann solver |
| **Dissipation Control** | Tunable parameters | Automatic from physics |
| **Kernel Argument** | $h$ | $\sqrt{2}h$ (from convolution) |
| **Density Approach** | Mass-based ($\rho^*$) | Volume-based ($V_p$) |
| **Contact Discontinuities** | Requires conductivity term | Natural handling |
| **Computational Cost** | Lower | Higher (Riemann solve) |
| **Accuracy at Shocks** | Depends on $\alpha_{\rm AV}$ tuning | Consistent high accuracy |
| **Conservation** | Guaranteed | Guaranteed by antisymmetry |

---

## 16. Implementation Notes

### 16.1 Numerical Stability

1. **Floor values:** Set minimum $n_{\rm min}$, $P_{\rm min}$ to prevent division by zero
2. **Velocity limiter:** Ensure $|V| < 1$ (subluminal)
3. **Enthalpy bounds:** $w \geq 1$ always

### 16.2 Parallelization

- Particle loop: embarrassingly parallel
- Pair loop: requires careful handling of race conditions
- Riemann solver: independent per pair

### 16.3 Boundary Conditions

- **Inflow/Outflow:** Inject/remove particles with specified primitives
- **Reflecting:** Mirror particle positions and velocities
- **Periodic:** Standard ghost particle approach
- **Black hole:** Remove particles crossing horizon ($r < r_H$)

---

## 17. References

1. Kitajima, K., Inutsuka, S., & Seno, I. (2025). "Special Relativistic Smoothed Particle Hydrodynamics Based on Riemann Solver." *Journal of Computational Physics*.

2. Liptai, D., & Price, D. J. (2019). "General relativistic smoothed particle hydrodynamics." *MNRAS*, 485(1), 819-842.

3. Monaghan, J. J., & Price, D. J. (2001). "Variational principles for relativistic smoothed particle hydrodynamics." *MNRAS*, 328(2), 381-392.

4. Pons, J. A., Martí, J. M., & Müller, E. (2000). "The exact solution of the Riemann problem with non-zero tangential velocities in relativistic hydrodynamics." *Journal of Fluid Mechanics*, 422, 125-139.

5. Inutsuka, S. (2002). "Reformulation of Smoothed Particle Hydrodynamics with Riemann Solver." *Journal of Computational Physics*, 179(1), 238-267.

---

## Appendix A: Kernel Gradient in Cartesian Coordinates

For the Gaussian kernel:
$$\nabla_k W(\mathbf{x}_i - \mathbf{x}_j, h) = -\frac{2(x^k_i - x^k_j)}{h^2} W(\mathbf{x}_i - \mathbf{x}_j, h)$$

With smoothing length $\sqrt{2}h$:
$$\nabla_k W_{ij}(\sqrt{2}h) = -\frac{(x^k_i - x^k_j)}{h^2} W_{ij}(\sqrt{2}h)$$

---

## Appendix B: V² Interpolation

The volume integral term:
$$V^2_{ij}(h) = \int \frac{1}{N^{*2}(\mathbf{x})} \left(\frac{\sqrt{2}}{h\sqrt{\pi}}\right)^d \exp\left[-\frac{2(\mathbf{x} - \frac{\mathbf{x}_i+\mathbf{x}_j}{2})^2}{h^2}\right] d^3x$$

This is approximated by interpolating $N^{*-2}$ (linear or cubic spline) and integrating analytically. See Inutsuka (2002) for details.

Symmetric average:
$$V^2_{ij} = \frac{1}{2}\left(V^2_{ij}(h_i) + V^2_{ij}(h_j)\right)$$

---

## Appendix C: Stress-Energy Tensor Components

For computing source terms, the relevant combination is:

$$\frac{\sqrt{-g} T^{\mu\nu}}{2N^*} = \frac{1}{2}\left[H U^0 v^\mu v^\nu + \frac{P g^{\mu\nu}}{n U^0}\right]$$

Explicitly:
$$\frac{\sqrt{-g} T^{00}}{2N^*} = \frac{1}{2}\left[H U^0 + \frac{P g^{00}}{n U^0}\right]$$
$$\frac{\sqrt{-g} T^{0i}}{2N^*} = \frac{1}{2}\left[H U^0 v^i + \frac{P g^{0i}}{n U^0}\right]$$
$$\frac{\sqrt{-g} T^{ij}}{2N^*} = \frac{1}{2}\left[H U^0 v^i v^j + \frac{P g^{ij}}{n U^0}\right]$$
