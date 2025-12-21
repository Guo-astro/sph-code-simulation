# GR-GSPH Formulation: General Relativistic Godunov SPH

## Overview

GR-GSPH (General Relativistic Godunov SPH) combines:
1. **GRSPH** (General Relativistic SPH) framework from Liptai & Price (2019)
2. **Godunov** SPH method with exact Riemann solver from Kitajima et al. (2025)

The key advantage is that the exact Riemann solver replaces artificial viscosity, providing physically correct entropy generation at shocks only.

---

## Are We Solving Einstein's Equations?

**NO.** We are NOT solving Einstein's field equations:

$$G_{\mu\nu} = 8\pi T_{\mu\nu}$$

Instead, we use a **prescribed (fixed) spacetime metric**. The metric $g_{\mu\nu}$ is given analytically (e.g., Minkowski, Schwarzschild, Kerr), and we evolve only the **hydrodynamic variables** (density, velocity, pressure) in this fixed curved background.

This is called the **test-fluid approximation**:
- Fluid mass is negligible compared to the central object (e.g., black hole)
- The fluid's stress-energy doesn't backreact on the spacetime
- Appropriate for: accretion disks, jets, neutron star surfaces, etc.

If we were solving Einstein's equations dynamically, we would need:
- Evolution of the metric variables (lapse, shift, 3-metric)
- Constraint equations (Hamiltonian, momentum)
- Gauge choices (1+log slicing, Gamma-driver shift)
- This is "numerical relativity" - a much more complex problem

---

## Why 3+1 Decomposition?

The 3+1 (ADM) decomposition is how we describe curved spacetime for numerical evolution. It splits 4D spacetime into 3D spatial slices evolving in time.

### The 3+1 Metric

The spacetime line element is written as:

$$ds^2 = -\alpha^2 dt^2 + \gamma_{ij}(dx^i + \beta^i dt)(dx^j + \beta^j dt)$$

where:
- **$\alpha$** = lapse function: time dilation between coordinate time and proper time
- **$\beta^i$** = shift vector: how spatial coordinates drift between time slices
- **$\gamma_{ij}$** = spatial 3-metric: geometry of each time slice

### Why We Need It

Even though we don't solve Einstein's equations, we need 3+1 variables because:

1. **Conserved variables** are defined with respect to the normal observer:
   - $N = n\Gamma\sqrt{\gamma}$ (lab-frame baryon density, where $\Gamma$ is Lorentz factor)
   - $S_i = n h \Gamma^2 V_i$ (covariant momentum)
   - $e = n h \Gamma^2 - P - n\Gamma$ (energy density)

2. **Metric appears in the equations**:
   - Lorentz factor: $\Gamma = \alpha u^0 = (1 - \gamma_{ij}V^iV^j)^{-1/2}$
   - Volume element: $\sqrt{-g} = \alpha\sqrt{\gamma}$
   - Gravitational source: $f_i = \frac{\sqrt{-g}}{2\rho^*} T^{\mu\nu} \partial_i g_{\mu\nu}$

3. **Source terms require metric derivatives**:
   - $\partial_i \alpha$: gravitational "acceleration"
   - $\partial_i \beta^j$: frame-dragging effects
   - $\partial_i \gamma_{jk}$: spatial curvature effects

---

## Minkowski Spacetime: The SR Limit

In Minkowski (flat) spacetime:

$$\alpha = 1, \quad \beta^i = 0, \quad \gamma_{ij} = \delta_{ij}$$

All metric derivatives vanish:

$$\partial_i \alpha = 0, \quad \partial_i \beta^j = 0, \quad \partial_i \gamma_{jk} = 0$$

Therefore:
- **No gravitational source terms**: $f_i = 0$
- **No frame-dragging**: $\beta^i = 0$
- **Euclidean geometry**: $\gamma_{ij} = \delta_{ij}$

**GR-GSPH reduces EXACTLY to SR-GSPH in Minkowski spacetime.**

---

## Mathematical Equivalence: GR-GSPH (Minkowski) = SR-GSPH

### Conserved Variables

| Variable | SR-GSPH | GR-GSPH (Minkowski) |
|----------|---------|---------------------|
| Lab density | $N = n\Gamma$ | $N = n\Gamma\sqrt{\gamma} = n\Gamma$ |
| Momentum | $S_i = nhW^2v_i$ | $S_i = nhW^2V_i$ |
| Energy | $e = nhW^2 - P - nW$ | $e = nhW^2 - P - nW$ |

Where $W = \Gamma$ (Lorentz factor), and in Minkowski: $V^i = v^i$ (coordinate velocity = 3-velocity).

### Evolution Equations

**SR-GSPH** (Kitajima et al. 2025, Eq. 58-59):
$$\frac{dS_i}{dt} = -\sum_b \nu_b \left[ P^*_{ab} V_a^2 \Omega_a \nabla_a W_{ab} + P^*_{ab} V_b^2 \Omega_b \nabla_b W_{ab} \right]$$

$$\frac{de}{dt} = v^*_i \frac{dS_i}{dt}$$

**GR-GSPH** (with Minkowski metric):
$$\frac{dS_i}{dt} = -\sum_b \nu_b \left[ P^*_{ab} V_a^2 \Omega_a \nabla_a W_{ab} + P^*_{ab} V_b^2 \Omega_b \nabla_b W_{ab} \right] + \underbrace{f_i}_{=0}$$

The gravitational source term vanishes in Minkowski:
$$f_i = \frac{\sqrt{-g}}{2\rho^*} T^{\mu\nu} \partial_i g_{\mu\nu} = 0$$

### Riemann Solver

Both use the same exact SR Riemann solver (Kitajima formulation):

1. **Input**: Left/Right states $(v^x, n, P, c_s, v^t)$
2. **Algorithm**: Newton-Raphson iteration for $P^*$, using shock/rarefaction wave relations
3. **Output**: Interface state $(P^*, v^x_*, v^t_*)$

The tangent velocity follows the K-invariant (Pons et al. 2000):
$$K = h W v^t = \text{const across rarefactions}$$

### Code Comparison

**SR-GSPH** (`sr_fluid_force.cpp:791-793`):
```cpp
const vec_t force = grad_W_i * (-P_star * V_i * V_i * omega_i)
                  + grad_W_j * (-P_star * V_j * V_j * omega_j);
const real power = inner_product(v_star_vec, force);
```

**GR-GSPH** (`gr_fluid_force.cpp:352-356`):
```cpp
const vec_t force = grad_W_i * (-P_star * V_i * V_i * omega_i)
                  + grad_W_j * (-P_star * V_j * V_j * omega_j);
const real power = inner_product(v_star_vec, force);
```

**Identical force calculation!** The only difference is GR-GSPH adds gravitational source terms which are zero in Minkowski.

---

## Riemann Problem Structure

The special relativistic Riemann problem has the wave structure:

```
     Left         Rarefaction    Contact     Shock       Right
     State           Fan       Discontinuity  Wave       State

    (n_L,v_L,P_L)  ←----→    |    *    |  ----→    (n_R,v_R,P_R)
                                (P*,v*)
                   λ_-         λ_cd=v*        λ_+
```

Wave types depend on pressure jump:
- **Shock**: $P^* > P_{ahead}$ (compression)
- **Rarefaction**: $P^* < P_{ahead}$ (expansion)

For Rosswog Test 1: $(n_L, P_L) = (10, 40/3)$, $(n_R, P_R) = (1, 10^{-6})$
- Left: Rarefaction wave
- Right: Shock wave
- Solution type: **RS** (Rarefaction-Shock)

---

## Gravitational Source Terms (for non-Minkowski)

In curved spacetime, the gravitational source term is:

$$f_i = \frac{\sqrt{-g}}{2\rho^*} T^{\mu\nu} \partial_i g_{\mu\nu}$$

Where the stress-energy tensor in 3+1 form:
$$T^{00} = \frac{\rho h \Gamma^2}{\alpha^2} - \frac{P}{\alpha^2}$$
$$T^{0i} = \frac{\rho h \Gamma^2 V^i}{\alpha}$$
$$T^{ij} = \rho h \Gamma^2 V^i V^j + P \gamma^{ij}$$

This captures:
- Gravitational attraction (from $\partial_i \alpha$)
- Frame-dragging (from $\partial_i \beta^j$ for Kerr)
- Spatial curvature effects (from $\partial_i \gamma_{jk}$)

---

## Schwarzschild Radial Shock Tube Test

### Test Configuration

A 1D shock tube aligned with the radial direction from a Schwarzschild black hole:

- **Domain**: $r \in [3M, 15M]$ (outside the horizon $r_s = 2M$)
- **Discontinuity**: $r = 6M$
- **Initial conditions**: Rosswog Test 1
  - Left ($r < 6M$): $(n, v, P) = (10, 0, 40/3)$
  - Right ($r > 6M$): $(n, v, P) = (1, 10^{-6})$

### Schwarzschild Metric

In Cartesian-like coordinates (GRSPH Appendix A):

$$\alpha = \sqrt{1 - \frac{2M}{r}}, \quad \beta^i = 0, \quad \sqrt{-g} = \alpha$$

Key properties:
- Lapse $\alpha < 1$ causes gravitational time dilation
- Near $r = 2M$: $\alpha \to 0$ (horizon)
- Metric derivatives give gravitational source terms

### Expected Physics

1. **Gravitational infall**: Particles drift toward the BH (negative velocity)
2. **Asymmetric wave propagation**:
   - Outgoing waves are slowed (gravitational redshift)
   - Ingoing waves are accelerated (blueshift)
3. **Comparison with Minkowski**: Without source terms, waves would propagate symmetrically

### Why There's No Analytical Solution

**There is NO closed-form analytical solution** for the Riemann problem in Schwarzschild spacetime because:

1. **Position-dependent wave speeds**: Local sound speed $c_s^{\text{local}} = \alpha(r) \cdot c_s^{\text{proper}}$
2. **Not self-similar**: Cannot use $\xi = (r - r_0)/t$ as in flat spacetime
3. **Continuous source terms**: Gravitational acceleration modifies the solution at every timestep

### Verification Approach

Since no analytical solution exists, we verify correctness by:

1. **Minkowski limit test**: GR-GSPH with Minkowski metric should match exact SR solution
   - Result: L2 density error ~5%, pressure error ~3% ✓
2. **GR effects test**: Schwarzschild should differ from Minkowski
   - Result: Mean velocity shift ~-0.3 (infall toward BH) ✓

Run verification:
```bash
python3 simulations/relativistic/gr_gsph/scripts/grgsph_verification.py
```

### Running the Test

```bash
./build/sph simulations/relativistic/gr_gsph/config/presets/gr_schwarzschild_shock.json
python3 simulations/relativistic/gr_gsph/scripts/animate_schwarzschild_shock.py
```

Results: `simulations/relativistic/gr_gsph/results/schwarzschild_shock/`

---

## Dynamic Spacetime: BSSN vs Post-Minkowski

### Fixed Spacetime (Current Implementation)

We use a **prescribed metric** (Minkowski, Schwarzschild, Kerr) and only evolve the fluid. This is the test-fluid approximation.

### Dynamic Spacetime Options

**Option 1: BSSN (Grid-based)**
- Baumgarte-Shapiro-Shibata-Nakamura formulation
- Evolves $(\alpha, \beta^i, \tilde{\gamma}_{ij}, K_{ij}, ...)$
- Requires Eulerian grid + constraint equations
- NOT natural for SPH (would need hybrid grid+particle)

**Option 2: Post-Minkowski Expansion (Recommended for SPH)**
- Expand metric: $g_{\mu\nu} = \eta_{\mu\nu} + h_{\mu\nu}^{(1)} + h_{\mu\nu}^{(2)} + ...$
- $h^{(1)}$ = Newtonian potential
- $h^{(2)}$ = 1PN (post-Newtonian) corrections
- Compute metric from particle positions via Poisson-like equations
- Self-consistent: fluid → metric → fluid

This is essentially **post-Newtonian SPH** - suitable for compact binary inspirals, neutron star mergers, etc.

---

## Implementation Summary

| Feature | SR-GSPH | GR-GSPH |
|---------|---------|---------|
| Metric | Minkowski (implicit) | Any (Minkowski, Schwarzschild, Kerr) |
| Einstein Eqs | No | No |
| Conserved vars | $N, S_i, e$ | $N, S_i, e$ (with $\sqrt{\gamma}$) |
| Force | $-P^* V^2 \Omega \nabla W$ | $-P^* V^2 \Omega \nabla W + f_i$ |
| Riemann solver | Exact (Kitajima) | Exact (same) |
| Source terms | None | $f_i$ from metric derivatives |
| Minkowski limit | N/A | Reduces to SR-GSPH |

---

## References

1. **Kitajima, Inutsuka, Seno (2025)**: "Godunov SPH for special relativistic hydrodynamics" [arXiv:2510.18251v1]
2. **Liptai & Price (2019)**: "General relativistic smoothed particle hydrodynamics" [MNRAS 485, 819]
3. **Pons, Martí, Müller (2000)**: "Exact solution of the Riemann problem with non-zero tangential velocities" [J. Fluid Mech. 422, 125]
4. **Rosswog (2010)**: "Conservative, special-relativistic SPH" [J. Comp. Phys. 229, 8591]
5. **Rezzolla & Zanotti (2013)**: "Relativistic Hydrodynamics" (textbook)
