# Why GSPH Needs Ω Correction but Standard SPH Does Not (Revised)

## The Confusion

You correctly identified a flaw in my reasoning:

> "If standard SPH cancels the spurious force, won't the acceleration be weaker, so core collapse happens faster?"

This is an excellent point. My previous explanation was wrong. Let me reconsider from scratch.

**The real question is:** What is fundamentally different between standard SPH and GSPH that makes one stable and the other unstable without Ω?

---

## 1. Re-examining the Problem

### What Actually Happens in Simulations:

| Method | Without Ω | With Ω |
|--------|-----------|--------|
| **Standard SPH** | Stable ✓ | Stable ✓ |
| **GSPH** | Core collapse ✗ | Stable ✓ |

### The Real Question:

Why does GSPH **specifically** suffer from core collapse without Ω, while standard SPH remains stable?

---

## 2. The Correct Analysis: It's About the Density Sum, Not Just Forces

### 2.1 The Density Equation

Both methods compute density as:
$$\rho_i = \sum_j m_j W(r_{ij}, h_i)$$

When $h = h(\rho)$, this is an **implicit equation** — $\rho_i$ appears on both sides.

### 2.2 How Forces Use Density

**Standard SPH:**
$$\mathbf{a}_i = -\sum_j m_j \left( \frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2} \right) \nabla W_{ij}$$

Uses **particle densities** $\rho_i$, $\rho_j$ directly.

**GSPH (Inutsuka formulation):**
$$\mathbf{a}_i = -\sum_j m_j \int \frac{P(\mathbf{x})}{\rho^2(\mathbf{x})} \left[ \nabla W_i \cdot W_j - W_i \cdot \nabla W_j \right] d\mathbf{x}$$

Uses the **continuous density field** $\rho(\mathbf{x})$ inside the integral.

---

## 3. The Key Difference: Where the Ω Factor Acts

### 3.1 Standard SPH: Ω in Density Only

In standard SPH, the density sum is:
$$\rho_i = \sum_j m_j W(r_{ij}, h_i)$$

When we vary particle positions with $h = h(\rho)$, the Ω factor appears in:
$$\frac{\partial \rho_i}{\partial \mathbf{x}_k} = \frac{1}{\Omega_i} \sum_j m_j \nabla_k W_{ij}$$

But in the **force equation**, we use $\rho_i$ as a **fixed value** (computed from the density sum). The gradient $\nabla W_{ij}$ is evaluated at fixed particle positions.

**The Ω factor affects how density responds to motion, but the force formula itself uses the standard kernel gradient.**

### 3.2 GSPH: Ω in the Force Integral

In GSPH, the force involves an **integral over continuous fields**:
$$\int \frac{P(\mathbf{x})}{\rho^2(\mathbf{x})} \nabla W(\mathbf{x} - \mathbf{x}_i, h[\mathbf{x}]) \, d\mathbf{x}$$

Here, $h[\mathbf{x}]$ varies **inside the integral**. The gradient $\nabla W$ must account for:
$$\nabla W(\mathbf{x} - \mathbf{x}_i, h[\mathbf{x}]) = \frac{\partial W}{\partial r} \nabla r + \frac{\partial W}{\partial h} \nabla h$$

**The spurious $\partial W / \partial h$ term appears directly in the force integral.**

---

## 4. The Physical Mechanism

### 4.1 Standard SPH: Discrete Particle Interactions

Standard SPH computes forces between **discrete particles**:
- Particle $i$ feels pressure from particle $j$
- The force depends on $P_i$, $P_j$, $\rho_i$, $\rho_j$, and $\nabla W_{ij}$
- The kernel gradient $\nabla W_{ij}$ is a function of separation $r_{ij}$ only

When $h$ varies:
- The kernel shape changes
- But **for a given pair**, the force still points along $\hat{\mathbf{r}}_{ij}$
- The h-variation just modifies **how strongly** nearby particles interact

**There's no spurious force direction** — only a modified interaction strength.

### 4.2 GSPH: Continuous Field Gradients

GSPH computes forces from **continuous pressure fields**:
- The force is an integral over the kernel overlap region
- The integrand involves $\nabla W$ with $h = h(\mathbf{x})$
- The gradient of a spatially-varying $h$ creates a **spurious force component**

When $h$ varies:
$$\nabla W = \frac{\partial W}{\partial r} \hat{\mathbf{r}} + \frac{\partial W}{\partial h} \nabla h$$

The second term is **not parallel to $\hat{\mathbf{r}}$** in general!

**This creates a spurious force component that doesn't exist in standard SPH.**

---

## 5. Detailed Derivation

### 5.1 Standard SPH Force

$$\mathbf{a}_i^{\text{SPH}} = -\sum_j m_j \left( \frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2} \right) \nabla_i W(r_{ij}, \bar{h})$$

where $\bar{h} = (h_i + h_j)/2$ or similar.

The gradient is:
$$\nabla_i W(r_{ij}, \bar{h}) = W'(r_{ij}, \bar{h}) \hat{\mathbf{r}}_{ij}$$

This **always points along $\hat{\mathbf{r}}_{ij}$**, regardless of how $h$ varies!

The h-variation only affects:
1. The magnitude $W'(r, h)$
2. Which neighbors are included (compact support)

**There is no spurious tangential or radial-bias force.**

### 5.2 GSPH Force (Simplified Discrete Form)

In practice, GSPH uses:
$$\mathbf{a}_i^{\text{GSPH}} = -\sum_j m_j P^*_{ij} \left( \frac{\nabla W_{ij}^i}{\rho_i^2} + \frac{\nabla W_{ij}^j}{\rho_j^2} \right)$$

where $W_{ij}^i = W(r_{ij}, h_i)$ and $W_{ij}^j = W(r_{ij}, h_j)$.

Now consider what happens when $h_i \neq h_j$:
$$\nabla W_{ij}^i = W'(r_{ij}, h_i) \hat{\mathbf{r}}_{ij}$$
$$\nabla W_{ij}^j = W'(r_{ij}, h_j) \hat{\mathbf{r}}_{ij}$$

Both still point along $\hat{\mathbf{r}}_{ij}$! So where's the problem?

### 5.3 The Real Issue: Asymmetric Weighting

The problem is the **different weighting** of the two terms.

Define:
$$\mathbf{a}_i = -\sum_j m_j P^* \left( \frac{W'_i}{\rho_i^2} + \frac{W'_j}{\rho_j^2} \right) \hat{\mathbf{r}}_{ij}$$

In equilibrium, we need:
$$\sum_j m_j P^* \left( \frac{W'_i}{\rho_i^2} + \frac{W'_j}{\rho_j^2} \right) \hat{\mathbf{r}}_{ij} = -\mathbf{g}_i$$

But $W'_i = W'(r, h_i)$ and $W'_j = W'(r, h_j)$ have **different magnitudes** when $h_i \neq h_j$!

### 5.4 The Asymmetry Creates Net Inward Force

In a self-gravitating sphere:
- Core: high $\rho$ → small $h$
- Edge: low $\rho$ → large $h$

For a particle pair with $i$ closer to core, $j$ closer to edge:
- $h_i < h_j$
- $W'(r, h_i) > W'(r, h_j)$ (smaller h → steeper kernel → larger gradient)

The force term:
$$\frac{W'_i}{\rho_i^2} + \frac{W'_j}{\rho_j^2}$$

is **dominated by the high-density side** (which has larger $W'$ despite smaller $1/\rho^2$).

This creates a **systematic bias toward the core**.

---

## 6. Standard SPH: Why It's Immune

### 6.1 The Shared Kernel

In standard SPH:
$$\mathbf{a}_i = -\sum_j m_j \left( \frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2} \right) W'(\bar{h}) \hat{\mathbf{r}}_{ij}$$

The kernel gradient $W'(\bar{h})$ is evaluated at the **average** smoothing length.

**Both particles contribute through the same $W'$!**

The asymmetry in $P_i/\rho_i^2$ vs $P_j/\rho_j^2$ is a **physical** asymmetry (pressure gradient), not a numerical artifact.

### 6.2 The Cancellation

Consider the sum over all neighbors:
$$\mathbf{a}_i = -\sum_j m_j \left( \frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2} \right) W'_{ij} \hat{\mathbf{r}}_{ij}$$

The $P_i/\rho_i^2$ term:
$$-\frac{P_i}{\rho_i^2} \sum_j m_j W'_{ij} \hat{\mathbf{r}}_{ij}$$

For a particle in a smooth density field, this sum is approximately:
$$\sum_j m_j W'_{ij} \hat{\mathbf{r}}_{ij} \approx -\nabla \rho / \text{(some factor)}$$

The $P_j/\rho_j^2$ term similarly gives the pressure gradient.

**There's no systematic bias because both terms use the same kernel!**

---

## 7. The Correct Interpretation

### 7.1 Standard SPH

- Uses **one kernel per pair** (averaged or symmetric)
- h-variation changes **interaction strength**, not **force direction**
- No systematic radial bias
- Ω correction improves accuracy but isn't essential for stability

### 7.2 GSPH

- Uses **two kernels per pair** (one for each particle's contribution)
- h-variation creates **asymmetric weighting**
- Smaller h (high density) → steeper kernel → dominates the sum
- **Systematic bias toward high-density regions**
- Ω correction is **essential** to restore symmetry

---

## 8. Mathematical Proof of the Bias

### 8.1 GSPH Without Ω

Force on particle $i$:
$$\mathbf{a}_i = -\sum_j m_j P^* \left( \frac{W'(r, h_i)}{\rho_i^2} + \frac{W'(r, h_j)}{\rho_j^2} \right) \hat{\mathbf{r}}_{ij}$$

For a kernel $W(r, h) = h^{-d} w(r/h)$:
$$W'(r, h) = \frac{\partial W}{\partial r} = h^{-d-1} w'(r/h)$$

So:
$$\frac{W'(r, h_i)}{W'(r, h_j)} = \left( \frac{h_j}{h_i} \right)^{d+1} \frac{w'(r/h_i)}{w'(r/h_j)}$$

When $h_i < h_j$ (particle $i$ in denser region):
$$\frac{W'(r, h_i)}{W'(r, h_j)} > 1$$

**The denser side contributes more to the force!**

### 8.2 The Net Effect

In equilibrium, the pressure force should balance gravity:
$$\mathbf{a}^P + \mathbf{g} = 0$$

But the GSPH pressure force is **biased toward high-density regions**.

For a particle at intermediate radius:
- Inner neighbors (high $\rho$, small $h$): contribute **stronger outward force**
- Outer neighbors (low $\rho$, large $h$): contribute **weaker inward force**

Wait... this would mean the pressure is **stronger**, which would cause **expansion**, not collapse!

### 8.3 Let Me Reconsider...

Actually, the force direction is always $\hat{\mathbf{r}}_{ij}$, pointing from $j$ to $i$.

For particle $i$ in the core:
- Most neighbors are at larger $r$ (outward)
- The kernel gradient points **inward** (toward $i$)
- Pressure force pushes $i$ **outward**

The question is: is this outward push too strong or too weak?

If $W'(r, h_i) > W'(r, h_j)$:
- The $i$-term (using $h_i$) dominates
- The $i$-term is $\propto 1/\rho_i^2$ (small, because $\rho_i$ is large)
- The $j$-term is $\propto 1/\rho_j^2$ (large, because $\rho_j$ is small)

So actually: **the outer particle's contribution ($j$-term) dominates** because $1/\rho_j^2 \gg 1/\rho_i^2$ even though $W'_i > W'_j$.

Hmm, this is getting complicated. Let me think differently.

---

## 9. A Cleaner Approach: The Ω Correction Origin

### 9.1 What Ω Actually Corrects

The Ω correction arises from the **variational principle**:
$$\frac{\partial \rho_i}{\partial \mathbf{x}_k} = \frac{1}{\Omega_i} \sum_j m_j \nabla_k W_{ij}$$

In standard SPH, this affects the **density response to motion**.

In GSPH, the force formula explicitly involves **integrals over density fields**, so the Ω factor must appear in the force itself.

### 9.2 Standard SPH: Ω is Implicit

Standard SPH uses:
$$\mathbf{a}_i = -\sum_j m_j \left( \frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2} \right) \nabla W_{ij}$$

The densities $\rho_i$, $\rho_j$ are computed from the density sum. When we evolve the system:
1. Particles move
2. Density is recomputed: $\rho_i = \sum_j m_j W_{ij}$
3. Forces are computed using new $\rho_i$

The Ω factor is **implicitly handled** by recomputing density at each step. The density "self-adjusts" to satisfy $h = h(\rho)$.

### 9.3 GSPH: Ω Must Be Explicit

GSPH computes forces from **continuous fields**. The Lagrangian derivation shows:
$$\mathbf{a}_i \propto \int \frac{P}{\rho^2} \nabla W \, d\mathbf{x}$$

When $h = h(\rho(\mathbf{x}))$, the variation of this integral gives extra terms involving $\partial W / \partial h$.

These terms **don't cancel** in GSPH because the two kernels (at $\mathbf{x}_i$ and $\mathbf{x}_j$) are treated separately.

**The Ω correction explicitly compensates for these terms.**

---

## 10. The Final Answer

### Why Standard SPH Works Without Ω:

1. **Single kernel per pair**: Both particles contribute through the same $\nabla W_{ij}$
2. **Implicit density update**: $\rho$ is recomputed each step, automatically adjusting to $h = h(\rho)$
3. **No continuous field integration**: Forces are discrete sums, not field integrals
4. **Symmetric structure**: The $(P_i/\rho_i^2 + P_j/\rho_j^2)$ form naturally balances

### Why GSPH Fails Without Ω:

1. **Two kernels per pair**: Each particle uses its own $h$ in $\nabla W$
2. **Continuous field structure**: The Lagrangian involves field integrals
3. **Asymmetric h-dependence**: Different kernels have different $\partial W / \partial h$
4. **Variational inconsistency**: Without Ω, the force doesn't derive from a consistent Lagrangian

### The Core Insight:

**Standard SPH's simpler structure is accidentally robust to h-variation.**

**GSPH's more sophisticated structure exposes the inconsistency, requiring explicit Ω correction.**

---

## 11. Analogy

### Standard SPH: A Simple Scale
- Two weights on a balance beam
- If the beam isn't perfectly level, both weights are equally affected
- The scale still balances

### GSPH: A Compound Lever
- Two weights on separate lever arms
- Each arm has different length (different h)
- Small errors in arm lengths create **net torque**
- Requires careful calibration (Ω correction)

---

---

## 12. The Real Answer: Two Kernels vs One Kernel

Looking at your actual GSPH implementation (`g_fluid_force.cpp`):

```cpp
const vec_t dw_i = kernel->dw(r_ij, r, h_i);       // Kernel gradient with h_i
const vec_t dw_j = kernel->dw(r_ij, r, p_j.sml);   // Kernel gradient with h_j

const vec_t f = dw_i * (p_j.mass * pstar * rho2_inv_i * omega_i)
              + dw_j * (p_j.mass * pstar * rho2_inv_j * omega_j);
```

**GSPH uses TWO kernel gradients per particle pair!**

### Standard SPH Force:
$$\mathbf{F}_{ij} = -m_i m_j \left(\frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2}\right) \nabla W(r_{ij}, \bar{h})$$

**One kernel gradient** $\nabla W(r_{ij}, \bar{h})$ shared by both terms.

### GSPH Force:
$$\mathbf{F}_{ij} = -m_j P^* \left(\frac{\nabla W(r_{ij}, h_i)}{\rho_i^2} + \frac{\nabla W(r_{ij}, h_j)}{\rho_j^2}\right)$$

**Two kernel gradients**: $\nabla W(r_{ij}, h_i)$ and $\nabla W(r_{ij}, h_j)$.

---

## 13. Why Two Kernels Cause Problems

### 13.1 The Magnitude Difference

For a kernel $W(r, h) = h^{-d} w(r/h)$:
$$|\nabla W(r, h)| = h^{-(d+1)} |w'(r/h)|$$

When $h_i < h_j$ (particle $i$ in denser region):
$$\frac{|\nabla W(r, h_i)|}{|\nabla W(r, h_j)|} = \left(\frac{h_j}{h_i}\right)^{d+1} \cdot \frac{|w'(r/h_i)|}{|w'(r/h_j)|}$$

The smaller kernel ($h_i$) has a **much steeper gradient**!

### 13.2 The Asymmetric Force

In GSPH, the force on particle $i$ is:
$$\mathbf{F}_i = -\sum_j m_j P^* \left(\frac{|\nabla W_i|}{\rho_i^2} + \frac{|\nabla W_j|}{\rho_j^2}\right) \hat{\mathbf{r}}_{ij}$$

For a particle in the core surrounded by less dense neighbors:
- $h_i < h_j$ → $|\nabla W_i| \gg |\nabla W_j|$
- $\rho_i > \rho_j$ → $1/\rho_i^2 < 1/\rho_j^2$

The two effects compete, but **the kernel gradient ratio wins**:
$$\frac{|\nabla W_i| / \rho_i^2}{|\nabla W_j| / \rho_j^2} = \left(\frac{h_j}{h_i}\right)^{d+1} \left(\frac{\rho_j}{\rho_i}\right)^2 \cdot (\text{kernel shape factor})$$

Using $h \propto \rho^{-1/d}$:
$$\left(\frac{h_j}{h_i}\right)^{d+1} = \left(\frac{\rho_i}{\rho_j}\right)^{(d+1)/d}$$

The ratio becomes:
$$\left(\frac{\rho_i}{\rho_j}\right)^{(d+1)/d - 2} = \left(\frac{\rho_i}{\rho_j}\right)^{(1-d)/d}$$

In 3D ($d=3$): exponent is $(1-3)/3 = -2/3$

So:
$$\frac{\text{i-term}}{\text{j-term}} = \left(\frac{\rho_i}{\rho_j}\right)^{-2/3} = \left(\frac{\rho_j}{\rho_i}\right)^{2/3}$$

**The low-density (outer) term dominates!**

### 13.3 What This Means for Force Balance

For a core particle:
- The outer neighbors' terms ($j$-terms with large $h_j$) dominate
- These terms push **inward** (from outer to inner)
- The net pressure force is **weaker than it should be**
- Gravity wins → **core collapse**

---

## 14. Why Standard SPH Doesn't Have This Problem

### 14.1 Shared Kernel Gradient

In standard SPH:
$$\mathbf{F}_i = -\sum_j m_j \left(\frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2}\right) \nabla W(r_{ij}, \bar{h})$$

**Both pressure terms use the same $\nabla W$!**

The kernel gradient $\nabla W(r_{ij}, \bar{h})$ is evaluated at the **average** $\bar{h} = (h_i + h_j)/2$.

### 14.2 Balanced Weighting

The ratio of terms is simply:
$$\frac{P_i / \rho_i^2}{P_j / \rho_j^2}$$

This is a **physical ratio** determined by the equation of state, not corrupted by kernel gradient differences.

### 14.3 The Ω Effect is Uniform

If we add Ω to standard SPH:
$$\mathbf{F}_i = -\sum_j m_j \left(\frac{P_i}{\rho_i^2 \Omega_i} + \frac{P_j}{\rho_j^2 \Omega_j}\right) \nabla W$$

The Ω factors modify the **pressure weighting**, but **not the kernel gradient**.

Without Ω, the error is:
- In the pressure coefficients only
- Symmetric in structure
- Affects magnitude, not the relative weighting between particles

---

## 15. The Core Insight

### Standard SPH:
```
F = (P_i/ρ_i² + P_j/ρ_j²) × [ONE shared ∇W]
     └─────────────────┘    └───────────────┘
      Physical weighting    Same for both terms
```

The h-variation in $\nabla W$ affects **both terms equally** → no relative bias.

### GSPH:
```
F = P* × (∇W_i/ρ_i² + ∇W_j/ρ_j²)
          └────────┘  └────────┘
          Uses h_i    Uses h_j
          (steep)     (shallow)
```

When $h_i \neq h_j$, the kernels have **different steepnesses** → systematic bias.

---

## 16. Final Summary

| Feature | Standard SPH | GSPH |
|---------|-------------|------|
| Kernel gradients per pair | **1** (shared $\bar{h}$) | **2** (separate $h_i$, $h_j$) |
| h-variation affects | Both terms equally | Each term differently |
| Relative term weighting | Physical ($P/\rho^2$ ratio) | Corrupted by kernel ratio |
| Without Ω | Small error, stable | Large bias, **collapse** |
| With Ω | More accurate | Required for stability |

### The Fundamental Reason:

**GSPH's use of separate kernels for each particle creates an asymmetric weighting that systematically underestimates the pressure force in high-density regions. Standard SPH's shared kernel maintains symmetric weighting regardless of h-variation.**

The Ω correction in GSPH compensates for this asymmetry by adjusting each term's contribution: $\omega_i = 1/\Omega_i$ and $\omega_j = 1/\Omega_j$ rebalance the kernels to restore correct force weighting.

---

## References

1. Inutsuka, S. (2002). JCP, 179, 238.
2. Cha, S.-H. & Whitworth, A.P. (2003). MNRAS, 340, 73.
3. Springel, V. & Hernquist, L. (2002). MNRAS, 333, 649.
4. Hopkins, P.F. (2013). MNRAS, 428, 2840.
