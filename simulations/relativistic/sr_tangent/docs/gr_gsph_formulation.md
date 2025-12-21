# GR-GSPH: General Relativistic Godunov SPH Formulation

## Overview

This document presents the mathematical formulation for GR-GSPH (General Relativistic Godunov SPH), which combines:
1. **GRSPH framework** (Liptai & Price 2019): 3+1 decomposition, conserved variables, metric terms
2. **Godunov approach** (Kitajima et al. 2025): Exact Riemann solver for interface states instead of artificial viscosity

---

## 1. Metric and 3+1 Decomposition

### 1.1 Line Element

The spacetime metric in 3+1 form:

$$ds^2 = -\alpha^2 dt^2 + \gamma_{ij}(dx^i + \beta^i dt)(dx^j + \beta^j dt)$$

where:
- $\alpha$ = lapse function
- $\beta^i$ = shift vector
- $\gamma_{ij}$ = spatial 3-metric

### 1.2 Key Metric Relations

**4-metric determinant:**
$$\sqrt{-g} = \alpha \sqrt{\gamma}$$

**Contravariant spatial metric:**
$$\gamma^{ij}\gamma_{jk} = \delta^i_k$$

**4-velocity normalization:**
$$g_{\mu\nu} U^\mu U^\nu = -1$$

### 1.3 Velocities

**Coordinate velocity:**
$$v^i = \frac{dx^i}{dt}$$

**Eulerian observer velocity (physical velocity):**
$$V^i = \frac{1}{\alpha}(v^i + \beta^i)$$

**Lorentz factor:**
$$\Gamma = \frac{1}{\sqrt{1 - V_i V^i}} = \frac{1}{\sqrt{1 - \gamma_{ij} V^i V^j}}$$

**4-velocity components:**
$$U^0 = \frac{\Gamma}{\alpha}, \quad U^i = \Gamma(V^i - \frac{\beta^i}{\alpha})$$

---

## 2. Conserved Variables

### 2.1 GRSPH Conserved Variables

Following Liptai & Price (2019):

**Conserved baryon density:**
$$\rho^* = \sqrt{-g} \, n \, U^0 = \sqrt{\gamma} \, n \, \Gamma$$

**Covariant momentum:**
$$p_i = U^0 w g_{i\mu} v^\mu = \Gamma w (\alpha V_i + g_{it}/\alpha)$$

For spatial indices with no time-space metric coupling ($g_{it} = 0$):
$$p_i = \Gamma w \alpha V_i = \Gamma w \gamma_{ij} V^j \alpha$$

**Total specific energy:**
$$e = U^0 w - P\sqrt{-g}/\rho^*$$

where $w = nh = n(1 + u + P/n)$ is the enthalpy density.

### 2.2 Entropy Variable (Alternative to Energy)

For shock problems, evolve entropy:
$$K = \frac{P}{n^\Gamma}$$

Evolution:
$$\frac{dK}{dt} = \frac{U^0 K}{u}\left[\text{dissipation terms}\right]$$

---

## 3. SPH Equations (GRSPH Form)

### 3.1 Continuity Equation

$$\frac{d\rho^*_a}{dt} = -\sum_b m_b (v^i_a - v^i_b) D^a_i$$

where the kernel gradient operator:
$$D^a_i = \frac{\sqrt{-g_a}}{\Omega_a} \frac{\partial W_{ab}(h_a)}{\partial x^i}$$

### 3.2 Momentum Equation

**Current GRSPH (with artificial viscosity):**
$$\frac{dp_{i,a}}{dt} = -\sum_b m_b \left[\frac{P_a + q_a}{(\rho^*_a)^2} D^a_i + \frac{P_b + q_b}{(\rho^*_b)^2} D^b_i\right] + f_i$$

where $q_a$ is artificial viscosity and $f_i$ is the gravitational source term:
$$f_i = \frac{\sqrt{-g}}{2\rho^*} T^{\mu\nu} \frac{\partial g_{\mu\nu}}{\partial x^i}$$

### 3.3 Energy Equation

$$\frac{de_a}{dt} = -\sum_b m_b \left[\frac{(P_a + q_a) v^i_a}{(\rho^*_a)^2} D^a_i + \frac{(P_b + q_b) v^i_b}{(\rho^*_b)^2} D^b_i\right] + \Lambda + \Pi_\text{cond}$$

---

## 4. GR-GSPH Modification: Godunov Approach

### 4.1 Core Idea

**Replace artificial viscosity $q_a$ with Riemann solver-based fluxes.**

At each particle pair (a,b), solve a local Riemann problem to find interface state $(P^*, v^{*}, n^*, v^{t*})$, then use these to compute fluxes.

### 4.2 Local Riemann Problem Setup

**Step 1: Project to line-of-sight frame**

Unit vector from a to b:
$$\hat{r}^i_{ab} = \frac{x^i_a - x^i_b}{|r_{ab}|}$$

Line-of-sight velocities:
$$V^*_a = \hat{r}^i_{ab} V_{i,a}, \quad V^*_b = \hat{r}^i_{ab} V_{i,b}$$

Tangent velocity (perpendicular component):
$$V^{t2}_a = V_i V^i - (V^*_a)^2$$

**Step 2: Define left/right states for Riemann problem**

Left state (particle a if approaching):
- $P_L = P_a$
- $n_L = n_a$
- $v^x_L = V^*_a$ (line-of-sight velocity)
- $v^t_L = \sqrt{V^{t2}_a}$ (tangent velocity magnitude)

Right state (particle b):
- $P_R = P_b$
- $n_R = n_b$
- $v^x_R = V^*_b$
- $v^t_R = \sqrt{V^{t2}_b}$

### 4.3 Solving the Riemann Problem

#### 4.3.1 Metric Considerations

For GR, the local Riemann problem is solved in the **local Eulerian observer frame** where:
- The metric appears locally Minkowski
- All velocity additions use relativistic formulas
- Sound speed: $c_s = \sqrt{\Gamma_\text{ad} P / (n w)}$

#### 4.3.2 Exact Solver (as in SR-GSPH)

**Root finding for $P^*$:**
$$f(P^*) = v^x_L(P^*) - v^x_R(P^*) = 0$$

**Shock relations** (when $P^* > P_\text{ahead}$):

From Taub adiabat:
$$H_b = \frac{-b + \sqrt{b^2 - 4ac}}{2a}$$

with coefficients:
- $a = 1 + \frac{(\Gamma_\text{ad} - 1)(P_a - P_b)}{\Gamma_\text{ad} P_b}$
- $b = -\frac{(\Gamma_\text{ad} - 1)(P_a - P_b)}{\Gamma_\text{ad} P_b}$
- $c = \frac{H_a (P_a - P_b)}{n_a} - H_a^2$

Mass flux:
$$j^2 = -\frac{P_b - P_a}{\frac{H_b}{n_b} - \frac{H_a}{n_a}}$$

**Rarefaction relations** (when $P^* < P_\text{ahead}$):

K-invariant:
$$A = h W v^t = \text{const}$$

Isentropic:
$$n_b = \left(\frac{P_b}{K_s}\right)^{1/\Gamma_\text{ad}}$$

Velocity via Gauss-Legendre integration:
$$v^x_b = \tanh\left[\text{arctanh}(v^x_a) + \text{sign} \cdot \int_{P_a}^{P_b} I(P) \, dP\right]$$

where:
$$I(P) = \frac{\sqrt{h^2 + A^2(1 - c_s^2)}}{(h^2 + A^2) \cdot n \cdot c_s}$$

### 4.4 Interface Flux Computation

Once $(P^*, v^{x*}, n^*, v^{t*})$ are found:

**Momentum flux contribution:**
$$F^p_{ab} = P^* \cdot A_{ab}$$

where $A_{ab}$ is the effective interface area (related to kernel gradient).

**Energy flux contribution:**
$$F^e_{ab} = P^* v^{x*} \cdot A_{ab}$$

### 4.5 Modified SPH Equations (GR-GSPH)

**Momentum equation:**
$$\frac{dp_{i,a}}{dt} = -\sum_b m_b \left[\frac{P^*_{ab}}{(\rho^*_a)^2} D^a_i + \frac{P^*_{ab}}{(\rho^*_b)^2} D^b_i\right] + f_i$$

Note: The interface pressure $P^*_{ab}$ replaces both $(P_a + q_a)$ and $(P_b + q_b)$.

**Energy equation:**
$$\frac{de_a}{dt} = -\sum_b m_b \left[\frac{P^*_{ab} v^{*}_{ab}}{(\rho^*_a)^2} D^a_i + \frac{P^*_{ab} v^{*}_{ab}}{(\rho^*_b)^2} D^b_i\right] + \Lambda$$

---

## 5. Practical Implementation Considerations

### 5.1 When to Solve Riemann Problem

Only for **approaching** particles:
$$\hat{r}^i_{ab}(v_a^i - v_b^i) < 0$$

For separating particles, use standard SPH pressure terms.

### 5.2 Fallback to HLLC

For robustness, use HLLC approximate solver when:
- Exact solver fails to converge
- Extreme velocity/pressure ratios
- Very weak waves (HLLC is faster)

### 5.3 Signal Speed for Timestep

Use exact Riemann solver wave speeds for timestep constraint:
$$v_\text{sig} = \max(|V_s^L|, |V_s^R|)$$

where $V_s^{L,R}$ are the shock/rarefaction wave speeds.

### 5.4 Gravitational Source Terms

The metric gradient term $f_i$ remains unchanged:
$$f_i = \frac{\sqrt{-g}}{2\rho^*} T^{\mu\nu} \frac{\partial g_{\mu\nu}}{\partial x^i}$$

This captures gravitational effects from curved spacetime.

---

## 6. Algorithm Summary

### 6.1 GR-GSPH Timestep

```
1. For each particle pair (a,b):
   a. Check if approaching: r_ab · (v_a - v_b) < 0
   b. If approaching:
      i.   Project velocities to line-of-sight frame
      ii.  Construct left/right states for Riemann problem
      iii. Solve Riemann problem for (P*, v*, n*, vt*)
      iv.  Store interface state for this pair
   c. If separating:
      - Use average pressure: P*_ab = (P_a + P_b)/2

2. Compute accelerations:
   a. Pressure gradient from interface states
   b. Gravitational source terms from metric gradients

3. Update conserved variables:
   a. Momentum: dp_i/dt
   b. Energy: de/dt (or entropy dK/dt)
   c. Position: dx^i/dt

4. Recover primitives:
   a. Newton-Raphson solve for enthalpy w
   b. Get n, P, v^i from conserved variables
```

### 6.2 Riemann Solver Call

```cpp
struct RiemannState {
    double P, n, vx, vt;  // Primitive variables
    double h, W;          // Derived quantities
};

struct RiemannResult {
    double P_star;        // Interface pressure
    double vx_star;       // Interface normal velocity
    double n_star;        // Interface density
    double vt_star;       // Interface tangent velocity
    double S_L, S_R;      // Wave speeds (for timestep)
};

RiemannResult solve_riemann_gr(
    const RiemannState& left,
    const RiemannState& right,
    double gamma_ad,
    const Metric& metric  // Local metric information
);
```

---

## 7. Comparison: GRSPH vs GR-GSPH

| Aspect | GRSPH (Liptai & Price) | GR-GSPH (Proposed) |
|--------|------------------------|---------------------|
| Shock capturing | Artificial viscosity $q$ | Exact Riemann solver |
| Contact resolution | Artificial conductivity | Natural from Riemann solution |
| Tangent velocity | May degrade at shocks | Preserved via K-invariant |
| Computational cost | Lower (no iteration) | Higher (iterative Riemann solve) |
| Accuracy | Good for smooth flows | Better for discontinuities |
| Robustness | Very stable | May need fallback solver |

---

## 8. Key Advantages of GR-GSPH

### 8.1 Better Shock Capturing

The exact Riemann solver provides:
- Correct jump conditions across shocks
- No tunable parameters (unlike $\alpha_\text{AV}$, $\beta_\text{AV}$)
- Proper entropy generation

### 8.2 Tangent Velocity Preservation

The K-invariant formulation ensures:
$$A = h W v^t = \text{const through rarefactions}$$

This prevents spurious acceleration/deceleration in tangent direction.

### 8.3 Natural Contact Discontinuity Treatment

The contact wave is explicitly captured by the Riemann solver:
- Pressure continuous: $P_L^* = P_R^*$
- Normal velocity continuous: $v^{x*}_L = v^{x*}_R$
- Tangent velocity and density can jump

### 8.4 No Artificial Viscosity Heating

In GRSPH, artificial viscosity causes excessive heating in:
- Geodesic Bondi flow
- Sonic point Bondi flow

GR-GSPH avoids this by only introducing entropy at true shocks.

---

## 9. Implementation Roadmap

### Phase 1: SR-GSPH Verification
1. Implement exact Riemann solver for SR (done in Kitajima tests)
2. Verify against known shock tube solutions
3. Test tangent velocity problems

### Phase 2: GR-GSPH Core
1. Extend Riemann solver for local metric frame
2. Integrate with GRSPH conserved variable framework
3. Implement cons2prim recovery with new flux terms

### Phase 3: GR-GSPH Tests
1. Schwarzschild shock tubes
2. Bondi accretion (compare heating to GRSPH)
3. Kerr metric tests

### Phase 4: Optimization
1. HLLC fallback for robustness
2. Caching of Riemann solutions for efficiency
3. Parallel Riemann solves

---

## 10. Detailed Implementation Plan

Based on the existing SR-GSPH codebase structure.

### 10.1 File Structure

The GR-GSPH module would extend the existing `srgsph/` directory:

```
src/
├── srgsph/                    # Existing SR-GSPH
│   ├── sr_exact_riemann.cpp   # Exact Riemann solver ✓
│   ├── sr_fluid_force.cpp     # HLLC + force computation ✓
│   ├── sr_primitive_recovery.cpp
│   └── sr_timestep.cpp
├── grgsph/                    # New GR-GSPH module
│   ├── gr_metric.hpp          # Metric classes (Schwarzschild, Kerr)
│   ├── gr_metric.cpp
│   ├── gr_riemann.hpp         # GR-aware Riemann solver wrapper
│   ├── gr_riemann.cpp
│   ├── gr_fluid_force.hpp     # GR fluid force with metric source terms
│   ├── gr_fluid_force.cpp
│   ├── gr_primitive_recovery.hpp
│   ├── gr_primitive_recovery.cpp
│   ├── gr_timestep.hpp
│   └── gr_timestep.cpp
```

### 10.2 Phase 1: Metric Infrastructure

**New file: `gr_metric.hpp/cpp`**

```cpp
namespace sph {
namespace grgsph {

// 3+1 metric decomposition
struct Metric31 {
    real alpha;           // Lapse function
    vec_t beta;           // Shift vector
    real gamma_ij[3][3];  // Spatial 3-metric
    real sqrt_gamma;      // sqrt(det(γ_ij))

    // Derived quantities
    real gamma_inv[3][3]; // Inverse spatial metric

    // Compute from position
    virtual void compute(const vec_t& pos) = 0;

    // Metric derivatives for source terms
    virtual void compute_derivatives(const vec_t& pos,
        real dg_tt[3], real dg_ij[3][3][3]) = 0;
};

class SchwarzschildMetric : public Metric31 {
    real M;  // Black hole mass
public:
    SchwarzschildMetric(real mass) : M(mass) {}
    void compute(const vec_t& pos) override;
    void compute_derivatives(const vec_t& pos, ...) override;
};

class KerrMetric : public Metric31 {
    real M;  // Black hole mass
    real a;  // Spin parameter
public:
    KerrMetric(real mass, real spin) : M(mass), a(spin) {}
    void compute(const vec_t& pos) override;
    void compute_derivatives(const vec_t& pos, ...) override;
};

} // namespace grgsph
} // namespace sph
```

### 10.3 Phase 2: GR-Aware Riemann Solver

**New file: `gr_riemann.hpp/cpp`**

Key insight: The Riemann problem is solved in the **local Eulerian observer frame** where the metric is locally Minkowski. This means:
1. Transform velocities from coordinate frame to Eulerian frame
2. Solve SR Riemann problem (reuse existing `sr_exact_riemann.cpp`)
3. Transform solution back to coordinate frame

```cpp
namespace sph {
namespace grgsph {

// Result of GR Riemann solve
struct GRRiemannResult {
    real P_star;      // Interface pressure
    real v_x_star;    // Interface velocity (Eulerian frame)
    real v_t_star;    // Interface tangent velocity
    real S_L, S_R;    // Wave speeds for timestep
};

/**
 * Solve Riemann problem in curved spacetime
 *
 * 1. Extract local metric at interface
 * 2. Transform velocities to Eulerian observer frame: V^i = (v^i + β^i)/α
 * 3. Compute Lorentz factor: Γ = 1/√(1 - γ_ij V^i V^j)
 * 4. Solve SR Riemann problem in Eulerian frame
 * 5. Return interface state in Eulerian frame (for flux computation)
 */
GRRiemannResult solve_riemann_gr(
    const GRParticle& left,
    const GRParticle& right,
    const Metric31& metric_L,
    const Metric31& metric_R,
    real gamma_ad
);

} // namespace grgsph
} // namespace sph
```

### 10.4 Phase 3: Conserved Variables and Primitive Recovery

**New file: `gr_primitive_recovery.hpp/cpp`**

Following GRSPH (Liptai & Price):

```cpp
namespace sph {
namespace grgsph {

// Conserved variables
struct GRConserved {
    real rho_star;     // = √γ · n · Γ
    vec_t p_i;         // Covariant momentum
    real e;            // Energy (or K for entropy)
};

// Primitive variables
struct GRPrimitive {
    real n;            // Rest-frame density
    real P;            // Pressure
    vec_t V;           // Eulerian velocity V^i
    real Gamma;        // Lorentz factor
};

/**
 * Recover primitives from conserved variables
 * Uses Newton-Raphson on enthalpy w
 *
 * Algorithm (from GRSPH Appendix B):
 * 1. Compute |p|² = γ^ij p_i p_j
 * 2. Initial guess: w₀ from previous timestep
 * 3. Newton iteration: f(w) = w(ρ,P) - w = 0
 * 4. Recover: Γ = √(1 + |p|²/w²), ρ = ρ*/Γ√γ, etc.
 */
GRPrimitive cons2prim(
    const GRConserved& U,
    const Metric31& metric,
    real gamma_ad,
    real w_prev  // Previous enthalpy (initial guess)
);

} // namespace grgsph
} // namespace sph
```

### 10.5 Phase 4: GR Fluid Force

**New file: `gr_fluid_force.hpp/cpp`**

```cpp
namespace sph {
namespace grgsph {

class GRFluidForce : public Integrator {
    // Metric for computing spacetime curvature
    std::unique_ptr<Metric31> m_metric;

public:
    void calculation(std::shared_ptr<Simulation> sim) override;

private:
    /**
     * Compute gravitational source term f_i
     *
     * f_i = (√-g / 2ρ*) · T^μν · ∂g_μν/∂x^i
     *
     * This captures:
     * - Gravitational attraction from metric curvature
     * - Frame-dragging effects (Kerr)
     */
    vec_t compute_metric_source(
        const GRParticle& p,
        const Metric31& metric
    );

    /**
     * Main force loop:
     * 1. For each particle pair:
     *    a. Compute local metrics at both positions
     *    b. Solve GR Riemann problem for interface state
     *    c. Compute pressure gradient force
     * 2. Add gravitational source term f_i
     */
    void compute_forces(std::vector<GRParticle>& particles);
};

} // namespace grgsph
} // namespace sph
```

### 10.6 Phase 5: Time Integration

Following GRSPH hybrid leapfrog:

```cpp
// Split acceleration into short-range (metric) and long-range (SPH)
//
// f_i^ext = metric source term (cheap, can substep)
// f_i^sph = pressure gradient (expensive, no substep)
//
// Hybrid algorithm:
// 1. Kick: p_i += (Δt/2) * f_i^sph
// 2. Substeps for external force:
//    a. Implicit momentum: p_i += (Δt_ext/2) * f_i^ext(p_i, x)
//    b. Implicit drift: x += (Δt_ext/2) * [v(p,x) + v(p,x_new)]
//    c. Kick: p_i += (Δt_ext/2) * f_i^ext
// 3. Final kick: p_i += (Δt/2) * f_i^sph
```

### 10.7 Test Cases

| Test | Description | Key Physics |
|------|-------------|-------------|
| Schwarzschild shock tube | Shock tube in Schwarzschild metric | SR + gravitational redshift |
| Geodesic Bondi | Pressureless accretion | Metric source terms, no viscous heating |
| Sonic point Bondi | Transonic accretion | Pressure + gravity balance |
| Kerr orbit | Circular orbit in Kerr | Frame dragging |
| GR blast wave | Spherical blast in Schwarzschild | Strong shocks in curved space |

### 10.8 Key Differences from GRSPH

| Aspect | GRSPH (original) | GR-GSPH (proposed) |
|--------|------------------|---------------------|
| Shock capturing | Artificial viscosity $q$ | Exact Riemann solver |
| Dissipation | $\alpha_{AV}$, $\alpha_u$ tunable | None (from jump conditions) |
| Heating in Bondi | Excessive (needs $\alpha_{AV}=0.1$) | Minimal (only at true shocks) |
| Contact discontinuity | Artificial conductivity | Natural from Riemann solver |
| Tangent velocity | May degrade | Preserved via K-invariant |

---

## References

1. Liptai, D., & Price, D. J. (2019). "General relativistic smoothed particle hydrodynamics." *MNRAS*, 485, 819-842.

2. Kitajima, Y., et al. (2025). "Special Relativistic Godunov SPH." *arXiv:2510.18251v1*.

3. Pons, J. A., Martí, J. M., & Müller, E. (2000). "The exact solution of the Riemann problem with non-zero tangential velocities." *J. Fluid Mech.*, 422, 125-139.

4. Mignone, A., & Bodo, G. (2005). "An HLLC Riemann solver for relativistic flows." *MNRAS*, 364, 126-136.
