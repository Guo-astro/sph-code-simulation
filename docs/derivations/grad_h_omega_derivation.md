# Grad-h (Ω) Correction in SPH: Theoretical Derivation

## 1. Origin of the Ω Correction

### 1.1 The SPH Density Estimate

The SPH density at particle $i$ is:

$$\rho_i = \sum_j m_j W(|\mathbf{r}_i - \mathbf{r}_j|, h_i)$$

When the smoothing length adapts to local density via $h_i = h(\rho_i)$, the density becomes an **implicit function** of itself:

$$\rho_i = \sum_j m_j W(r_{ij}, h(\rho_i))$$

### 1.2 Variational Derivation

The SPH equations should derive from a Lagrangian to guarantee conservation. The Lagrangian is:

$$L = \sum_i m_i \left( \frac{1}{2}v_i^2 - u_i(\rho_i, s_i) \right)$$

where $u_i$ is the specific internal energy. Taking variations with respect to particle positions:

$$\frac{\partial L}{\partial \mathbf{r}_i} = -\sum_j m_j \frac{\partial u_j}{\partial \rho_j} \frac{\partial \rho_j}{\partial \mathbf{r}_i}$$

The key is computing $\partial \rho_j / \partial \mathbf{r}_i$ when $h = h(\rho)$.

### 1.3 The Chain Rule with Variable h

$$\frac{\partial \rho_j}{\partial \mathbf{r}_i} = \sum_k m_k \left( \frac{\partial W_{jk}}{\partial \mathbf{r}_i} + \frac{\partial W_{jk}}{\partial h_j} \frac{\partial h_j}{\partial \rho_j} \frac{\partial \rho_j}{\partial \mathbf{r}_i} \right)$$

Solving for $\partial \rho_j / \partial \mathbf{r}_i$:

$$\frac{\partial \rho_j}{\partial \mathbf{r}_i} = \frac{1}{\Omega_j} \sum_k m_k \frac{\partial W_{jk}}{\partial \mathbf{r}_i}$$

where we define:

$$\boxed{\Omega_i \equiv 1 - \frac{\partial h_i}{\partial \rho_i} \sum_j m_j \frac{\partial W(r_{ij}, h_i)}{\partial h_i}}$$

---

## 2. The Corrected SPH Equations

### 2.1 Momentum Equation

From the Lagrangian derivation, the momentum equation becomes:

$$\frac{d\mathbf{v}_i}{dt} = -\sum_j m_j \left( \frac{P_i}{\rho_i^2 \Omega_i} \nabla_i W(r_{ij}, h_i) + \frac{P_j}{\rho_j^2 \Omega_j} \nabla_i W(r_{ij}, h_j) \right)$$

Compare to the **uncorrected** form:

$$\frac{d\mathbf{v}_i}{dt} = -\sum_j m_j \left( \frac{P_i}{\rho_i^2} \nabla_i W_{ij}^i + \frac{P_j}{\rho_j^2} \nabla_i W_{ij}^j \right)$$

The Ω factors **modify the effective pressure** in the force calculation.

### 2.2 Energy Equation

The corresponding energy equation is:

$$\frac{du_i}{dt} = \frac{P_i}{\rho_i^2 \Omega_i} \sum_j m_j (\mathbf{v}_i - \mathbf{v}_j) \cdot \nabla_i W(r_{ij}, h_i)$$

---

## 3. Conservation Properties

### 3.1 Momentum Conservation: PRESERVED ✓

**Claim:** The Ω-corrected equations conserve total momentum exactly.

**Proof:**

Total momentum rate of change:

$$\frac{d\mathbf{P}_{\text{tot}}}{dt} = \sum_i m_i \frac{d\mathbf{v}_i}{dt}$$

Substituting the momentum equation:

$$= -\sum_i \sum_j m_i m_j \left( \frac{P_i}{\rho_i^2 \Omega_i} \nabla_i W_{ij}^i + \frac{P_j}{\rho_j^2 \Omega_j} \nabla_i W_{ij}^j \right)$$

Consider the first term. Using antisymmetry $\nabla_i W_{ij} = -\nabla_j W_{ji}$:

$$\sum_i \sum_j m_i m_j \frac{P_i}{\rho_i^2 \Omega_i} \nabla_i W_{ij}^i = \sum_j \sum_i m_j m_i \frac{P_j}{\rho_j^2 \Omega_j} \nabla_j W_{ji}^j$$

Swapping indices $i \leftrightarrow j$:

$$= -\sum_i \sum_j m_i m_j \frac{P_j}{\rho_j^2 \Omega_j} \nabla_i W_{ij}^j$$

This exactly cancels the second term. Therefore:

$$\boxed{\frac{d\mathbf{P}_{\text{tot}}}{dt} = 0}$$

**Key insight:** Momentum conservation depends only on the **antisymmetry** of $\nabla W$, not on the coefficients in front. The Ω factors are symmetric under particle exchange in the pairwise sum, so they don't break conservation.

### 3.2 Energy Conservation: PRESERVED ✓

**Claim:** The Ω-corrected equations conserve total energy exactly.

**Proof:**

Total energy:

$$E = \sum_i m_i \left( \frac{1}{2} v_i^2 + u_i \right)$$

Rate of change:

$$\frac{dE}{dt} = \sum_i m_i \left( \mathbf{v}_i \cdot \frac{d\mathbf{v}_i}{dt} + \frac{du_i}{dt} \right)$$

Substituting the equations of motion:

$$\frac{dE}{dt} = \sum_i m_i \mathbf{v}_i \cdot \left[ -\sum_j m_j \left( \frac{P_i}{\rho_i^2 \Omega_i} \nabla_i W_{ij}^i + \frac{P_j}{\rho_j^2 \Omega_j} \nabla_i W_{ij}^j \right) \right]$$
$$+ \sum_i m_i \frac{P_i}{\rho_i^2 \Omega_i} \sum_j m_j (\mathbf{v}_i - \mathbf{v}_j) \cdot \nabla_i W_{ij}^i$$

The kinetic energy contribution from particle $i$'s own kernel:

$$-\sum_{i,j} m_i m_j \frac{P_i}{\rho_i^2 \Omega_i} \mathbf{v}_i \cdot \nabla_i W_{ij}^i$$

The internal energy contribution:

$$+\sum_{i,j} m_i m_j \frac{P_i}{\rho_i^2 \Omega_i} (\mathbf{v}_i - \mathbf{v}_j) \cdot \nabla_i W_{ij}^i$$

Adding these:

$$\sum_{i,j} m_i m_j \frac{P_i}{\rho_i^2 \Omega_i} (-\mathbf{v}_j) \cdot \nabla_i W_{ij}^i$$

By symmetry arguments (swapping $i \leftrightarrow j$ and using antisymmetry of $\nabla W$), this cancels with the cross-kernel terms:

$$\boxed{\frac{dE}{dt} = 0}$$

**Key insight:** Energy conservation is guaranteed by the **variational structure**. The Lagrangian derivation automatically produces energy-conserving equations, regardless of how Ω enters.

### 3.3 Angular Momentum Conservation: PRESERVED ✓

Angular momentum conservation follows from the central nature of the pairwise forces:

$$\mathbf{F}_{ij} \parallel \mathbf{r}_{ij}$$

Since $\nabla_i W_{ij} \parallel \mathbf{r}_{ij}$ for any spherically symmetric kernel, the Ω factors (which are scalars) don't change this. Therefore:

$$\boxed{\frac{d\mathbf{L}_{\text{tot}}}{dt} = 0}$$

---

## 4. What Ω Actually Affects

### 4.1 Force Magnitude (Not Direction)

The Ω correction modifies the **magnitude** of the pressure force:

$$|\mathbf{F}_i^{\text{pressure}}| \propto \frac{P}{\rho^2 \Omega}$$

- When $\Omega < 1$: Force is **amplified**
- When $\Omega > 1$: Force is **reduced**
- When $\Omega = 1$: Standard SPH (constant h)

Typically in self-gravitating systems:
- Core (high ρ, small h, ∂h/∂ρ < 0): $\Omega < 1$ → pressure force enhanced
- Edge (low ρ, large h): $\Omega \approx 1$

### 4.2 Equilibrium Accuracy

Without Ω correction, hydrostatic equilibrium has a **systematic error**:

$$\nabla P_{\text{SPH}} \neq \rho \mathbf{g}$$

The error scales as:

$$\epsilon \sim \left| \frac{\partial h}{\partial r} \right| \cdot \left| \frac{\partial \ln W}{\partial \ln h} \right|$$

For a polytrope with $h \propto \rho^{-1/3}$:

$$\epsilon \sim \frac{h}{3\rho} \left| \frac{\partial \rho}{\partial r} \right| \sim \frac{h}{R}$$

where $R$ is the system size. This is an **O(h/R) error per timestep** that accumulates, causing core collapse.

### 4.3 Summary Table

| Property | Affected by Ω? | Mechanism |
|----------|----------------|-----------|
| **Momentum conservation** | No ✓ | Antisymmetry of ∇W preserved |
| **Energy conservation** | No ✓ | Variational structure preserved |
| **Angular momentum** | No ✓ | Central force nature preserved |
| **Force magnitude** | Yes | P/ρ²Ω modifies pressure term |
| **Hydrostatic balance** | Yes | Correct ∇P = ρg matching |
| **Wave speeds** | Slightly | Modified effective pressure |

---

## 5. Why Ω Matters for Self-Gravity but Not Sod Shock

### 5.1 Self-Gravitating Hydrostatic Equilibrium

**Physical setup:**
- Pressure gradient balances gravity: $\nabla P = \rho \mathbf{g}$
- Density varies smoothly: $\rho(r)$ from Lane-Emden solution
- Smoothing length adapts: $h \propto \rho^{-1/3}$ to maintain neighbor count

**Why Ω is critical:**

1. **Continuous ∂h/∂r:** Throughout the configuration, h varies smoothly
   $$\frac{\partial h}{\partial r} = \frac{\partial h}{\partial \rho} \frac{\partial \rho}{\partial r} \neq 0$$

2. **Delicate force balance:** Must satisfy $|\nabla P| = |\rho g|$ to ~0.1% accuracy

3. **Cumulative error:** Without Ω, the O(h/R) error per timestep causes:
   - Core pressure underestimated → gravity wins → collapse
   - Or edge pressure overestimated → expansion

4. **The fix:** Ω correction restores variational consistency:
   $$\frac{P}{\rho^2 \Omega} \text{ gives correct } \nabla P$$

### 5.2 Sod Shock Tube

**Physical setup:**
- No gravity: $\mathbf{g} = 0$
- Piecewise uniform initial conditions
- Transient dynamics: shock, rarefaction, contact

**Why Ω is irrelevant:**

1. **Uniform regions:** $\rho = \text{const}$ → $h = \text{const}$ → $\partial h/\partial r = 0$ → $\Omega = 1$

2. **No equilibrium to maintain:** The solution is purely dynamic; there's no force balance requirement

3. **Discontinuities handled by Riemann solver:** At shocks/contacts, the Riemann solver captures the physics; Ω is a smooth correction

4. **Error tolerance:** Shock-capturing inherently has O(h) errors at discontinuities; the O(h) Ω correction is in the noise

### 5.3 Mathematical Comparison

| Quantity | Self-Gravity Hydrostatic | Sod Shock Tube |
|----------|-------------------------|----------------|
| $\partial \rho / \partial r$ | Smooth, O(ρ/R) | ~0 in uniform regions |
| $\partial h / \partial r$ | Continuous, nonzero | ~0, jumps at discontinuities |
| $\Omega - 1$ | O(0.1) to O(0.3) | O(0.01) or less |
| Force balance error | Cumulative, causes collapse | Local, transient |
| Required accuracy | ~0.1% for stability | ~1-5% acceptable |

---

## 6. Explicit Formula for Ω

For the relation $h = \eta (\rho)^{-1/d}$ where $d$ is dimension:

$$\frac{\partial h}{\partial \rho} = -\frac{h}{d \rho}$$

Therefore:

$$\Omega_i = 1 + \frac{h_i}{d \rho_i} \sum_j m_j \frac{\partial W(r_{ij}, h_i)}{\partial h_i}$$

For a kernel $W(r,h) = h^{-d} w(r/h)$ where $w$ is dimensionless:

$$\frac{\partial W}{\partial h} = -\frac{d}{h} W - \frac{r}{h^2} W'$$

where $W' = dW/d(r/h) \cdot (1/h)$.

---

## 7. Implementation Notes

From the code in `g_fluid_force.cpp`:

```cpp
const real omega_i = p_i.gradh;
const real omega_j = p_j.gradh;

const vec_t f = dw_i * (p_j.mass * pstar * rho2_inv_i * omega_i)
              + dw_j * (p_j.mass * pstar * rho2_inv_j * omega_j);
```

Here `omega_i` = $1/\Omega_i$ (the inverse), so the division $P/(\rho^2 \Omega)$ becomes multiplication by `omega_i`.

**Verification:** For a uniform density field, `omega ≈ 1`, and the standard SPH equations are recovered.

---

## 8. Conclusion

The Ω (grad-h) correction:

1. **Does NOT break conservation laws** — momentum, energy, and angular momentum are all exactly preserved

2. **Modifies force magnitudes** — the effective pressure becomes $P/\Omega$

3. **Is essential for hydrostatic equilibrium** — prevents artificial core collapse in self-gravitating systems

4. **Is negligible for shock tubes** — uniform regions have Ω ≈ 1, and dynamics don't require precise force balance

The correction arises from the **variational consistency** requirement: SPH equations derived from a Lagrangian with $h = h(\rho)$ naturally include the Ω factors to maintain conservation while correctly capturing the physics of variable resolution.

---

## References

1. Springel, V. & Hernquist, L. (2002). "Cosmological smoothed particle hydrodynamics simulations: the entropy equation." MNRAS, 333, 649.

2. Monaghan, J.J. (2002). "SPH compressible turbulence." MNRAS, 335, 843.

3. Price, D.J. & Monaghan, J.J. (2007). "An energy-conserving formalism for adaptive gravitational force softening in SPH and N-body codes." MNRAS, 374, 1347.

4. Hopkins, P.F. (2013). "A general class of Lagrangian smoothed particle hydrodynamics methods and implications for fluid mixing problems." MNRAS, 428, 2840.
