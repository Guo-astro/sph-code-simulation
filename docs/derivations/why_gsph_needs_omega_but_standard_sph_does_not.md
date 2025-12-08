# Why GSPH Needs Ω Correction but Standard SPH Does Not

## The Observation

Your simulation results show:
- **Standard SPH** (without grad-h): Maintains hydrostatic equilibrium ✓
- **GSPH** (without grad-h): Core collapse ✗
- **GSPH** (with grad-h Ω): Maintains hydrostatic equilibrium ✓

This seems paradoxical! Why would the "better" method (GSPH) fail where the "simpler" method (standard SPH) succeeds?

---

## 1. The Key Difference: Force Structure

### 1.1 Standard SPH Force

$$\mathbf{a}_i^{\text{SPH}} = -\sum_j m_j \left( \frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2} \right) \nabla_i W_{ij}$$

This is **symmetric in particles i and j** — both contribute their own $P/\rho^2$ term.

### 1.2 GSPH Force (Simplified)

$$\mathbf{a}_i^{\text{GSPH}} = -\sum_j m_j \frac{P^*_{ij}}{\bar{\rho}^2} \left( \nabla_i W_{ij}^i + \nabla_i W_{ij}^j \right)$$

where $P^*_{ij}$ is the Riemann solver pressure at the interface.

Or in the Inutsuka formulation:
$$\mathbf{a}_i^{\text{GSPH}} = -\sum_j m_j \int \frac{P(\mathbf{x})}{\rho^2(\mathbf{x})} \left[ \nabla W(\mathbf{x}-\mathbf{x}_i) W(\mathbf{x}-\mathbf{x}_j) - W(\mathbf{x}-\mathbf{x}_i) \nabla W(\mathbf{x}-\mathbf{x}_j) \right] d\mathbf{x}$$

This involves **convolution** — integrating over the kernel overlap region.

---

## 2. The Fundamental Difference: How $\nabla W$ Enters

### 2.1 Standard SPH: Single Kernel Gradient

In standard SPH, only **one** kernel gradient appears:
$$\nabla_i W_{ij} = \nabla_i W(|\mathbf{x}_i - \mathbf{x}_j|, h)$$

This gradient is evaluated at the **particle separation** $r_{ij}$, not integrated over space.

When $h$ varies:
$$\nabla_i W_{ij} = \frac{\partial W}{\partial r} \hat{\mathbf{r}}_{ij} + \frac{\partial W}{\partial h} \nabla_i h$$

But here's the key: **the spurious term is the same for both the $P_i/\rho_i^2$ and $P_j/\rho_j^2$ contributions!**

### 2.2 GSPH: Kernel Product and Convolution

In GSPH, the force involves products like:
$$\nabla W(\mathbf{x}-\mathbf{x}_i) \cdot W(\mathbf{x}-\mathbf{x}_j)$$

integrated over $\mathbf{x}$. This is fundamentally different because:
1. The convolution "smears" the pressure field
2. The gradient acts on a **field** $P(\mathbf{x})/\rho^2(\mathbf{x})$, not just particle values
3. The $h$-dependence enters the **integrand**, affecting the entire integral

---

## 3. Error Cancellation in Standard SPH

### 3.1 The Symmetric Structure

Standard SPH force:
$$\mathbf{F}_i = -\sum_j m_i m_j \left( \frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2} \right) \nabla_i W_{ij}$$

Consider the spurious contribution from $\partial W / \partial h \cdot \nabla h$:
$$\delta \mathbf{F}_i = -\sum_j m_i m_j \left( \frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2} \right) \frac{\partial W}{\partial h} \nabla h$$

### 3.2 The Cancellation Mechanism

In hydrostatic equilibrium, $P/\rho^2$ varies smoothly. For nearby particles:
$$\frac{P_i}{\rho_i^2} \approx \frac{P_j}{\rho_j^2} + O(h/R)$$

The spurious force on particle $i$ from particle $j$:
$$\delta \mathbf{F}_{ij} \propto \left( \frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2} \right) \frac{\partial W}{\partial h} \nabla h$$

And on particle $j$ from particle $i$:
$$\delta \mathbf{F}_{ji} \propto \left( \frac{P_j}{\rho_j^2} + \frac{P_i}{\rho_i^2} \right) \frac{\partial W}{\partial h} \nabla h$$

**These are equal and opposite!** The spurious forces cancel pairwise.

### 3.3 Why This Works

The standard SPH formulation has **built-in antisymmetry**:
$$\nabla_i W_{ij} = -\nabla_j W_{ji}$$

This antisymmetry is preserved even when $h$ varies, because both particles see the **same** kernel (evaluated at the same $r_{ij}$).

The spurious $\partial W / \partial h$ term:
- Has the same magnitude for both particles
- Acts in opposite directions (due to $\nabla_i = -\nabla_j$)
- **Cancels in the momentum sum**

---

## 4. Why GSPH Breaks This Cancellation

### 4.1 The Convolution Destroys Symmetry

In GSPH, the force involves:
$$\int \frac{P(\mathbf{x})}{\rho^2(\mathbf{x})} \nabla W(\mathbf{x}-\mathbf{x}_i, h(\mathbf{x})) W(\mathbf{x}-\mathbf{x}_j, h(\mathbf{x})) \, d\mathbf{x}$$

Here, $h(\mathbf{x})$ varies **inside the integral**. The spurious contribution is:
$$\int \frac{P}{\rho^2} \frac{\partial W_i}{\partial h} \nabla h \cdot W_j \, d\mathbf{x}$$

### 4.2 Asymmetric Weighting

The key difference: in GSPH, the kernels $W_i$ and $W_j$ are **different functions** centered at different particles.

The spurious term for particle $i$:
$$\delta \mathbf{F}_i \propto \int \frac{P}{\rho^2} \frac{\partial W_i}{\partial h} \nabla h \cdot W_j \, d\mathbf{x}$$

The spurious term for particle $j$:
$$\delta \mathbf{F}_j \propto \int \frac{P}{\rho^2} \frac{\partial W_j}{\partial h} \nabla h \cdot W_i \, d\mathbf{x}$$

These are **NOT equal and opposite** because:
- $\partial W_i / \partial h \neq \partial W_j / \partial h$ (different centers)
- The weighting by $W_j$ vs $W_i$ is asymmetric
- The integral samples different regions of the $P/\rho^2$ field

### 4.3 The Physical Picture

**Standard SPH:** "Each particle pair shares the same kernel, so they see the same spurious force (which cancels)."

**GSPH:** "Each particle has its own kernel, so the spurious forces don't match up and don't cancel."

---

## 5. Mathematical Proof

### 5.1 Standard SPH: Momentum Conservation Despite h-Variation

Total momentum change:
$$\frac{d\mathbf{P}}{dt} = \sum_i m_i \mathbf{a}_i = -\sum_i \sum_j m_i m_j \left( \frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2} \right) \nabla_i W_{ij}$$

Using $\nabla_i W_{ij} = -\nabla_j W_{ji}$ and swapping indices:
$$= -\frac{1}{2} \sum_i \sum_j m_i m_j \left( \frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2} \right) \left( \nabla_i W_{ij} - \nabla_j W_{ji} \right) = 0$$

**This holds even with the spurious $\partial W / \partial h$ terms**, because they also satisfy $\nabla_i = -\nabla_j$.

### 5.2 Standard SPH: Force Error Analysis

The net acceleration on particle $i$:
$$\mathbf{a}_i = -\sum_j m_j \left( \frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2} \right) \left( \frac{\partial W}{\partial r} \hat{\mathbf{r}}_{ij} + \frac{\partial W}{\partial h} \nabla h \right)$$

The spurious part:
$$\delta \mathbf{a}_i = -\sum_j m_j \left( \frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2} \right) \frac{\partial W}{\partial h} \nabla h$$

Now, $\nabla h$ at the midpoint is approximately:
$$\nabla h \approx \frac{h_i - h_j}{r_{ij}} \hat{\mathbf{r}}_{ij}$$

So the spurious term is **parallel to $\hat{\mathbf{r}}_{ij}$**, just like the physical force!

This means the spurious force:
1. Points in the same direction as the physical force
2. Is symmetric under $i \leftrightarrow j$
3. **Modifies the magnitude but not the structure of the force**

### 5.3 GSPH: Broken Symmetry

In GSPH with the Riemann solver:
$$\mathbf{a}_i = -\sum_j m_j P^*_{ij} \left( \frac{1}{\rho_i^2} \nabla W_{ij}^i + \frac{1}{\rho_j^2} \nabla W_{ij}^j \right)$$

The two gradient terms use **different kernels** ($h_i$ vs $h_j$).

When $h$ varies:
$$\nabla W_{ij}^i = \frac{\partial W^i}{\partial r} \hat{\mathbf{r}} + \frac{\partial W^i}{\partial h_i} \nabla h_i$$
$$\nabla W_{ij}^j = \frac{\partial W^j}{\partial r} \hat{\mathbf{r}} + \frac{\partial W^j}{\partial h_j} \nabla h_j$$

The spurious terms involve **different** $\partial W / \partial h$ (because $h_i \neq h_j$) and **different** $\nabla h$ (because they're evaluated at different locations).

**These don't cancel!**

---

## 6. The Error Magnitude

### 6.1 Standard SPH Without Ω

The spurious force modifies the **magnitude** of the pressure force, but:
- It points in the correct direction
- It's symmetric between particle pairs
- The net effect is a small rescaling of the effective pressure

Error in force balance:
$$\epsilon_{\text{SPH}} \sim \frac{h}{R} \times (\text{small geometric factor})$$

Because of the symmetric structure, much of the error cancels in the sum.

### 6.2 GSPH Without Ω

The spurious forces:
- Point in **different directions** for different terms
- Are **asymmetric** between particle pairs
- Create a **net spurious acceleration** toward high-density regions

Error in force balance:
$$\epsilon_{\text{GSPH}} \sim \frac{h}{R}$$

This error is **systematic** and doesn't cancel — it always points toward the core.

---

## 7. Physical Interpretation

### 7.1 Standard SPH: "Democratic" Pressure

In standard SPH:
$$\frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2}$$

Both particles contribute equally. The h-variation affects both contributions the same way.

**Analogy:** Two people pushing on opposite ends of a spring. If both pushes are slightly wrong by the same amount, the spring still balances.

### 7.2 GSPH: "Riemann" Pressure

In GSPH, the Riemann solver finds $P^*$ at the interface. But the **gradient weighting** uses separate kernels:
$$\frac{P^*}{\rho_i^2} \nabla W^i + \frac{P^*}{\rho_j^2} \nabla W^j$$

When $h_i \neq h_j$, the two gradient terms don't balance properly.

**Analogy:** Two people pushing on a spring, but using different-length levers. Even if the forces are equal, the torques don't balance.

---

## 8. Why Ω Fixes GSPH

### 8.1 The Ω Correction

$$\mathbf{a}_i = -\sum_j m_j P^* \left( \frac{1}{\rho_i^2 \Omega_i} \nabla W^i + \frac{1}{\rho_j^2 \Omega_j} \nabla W^j \right)$$

The Ω factors **reweight** each kernel's contribution to account for its h-dependence.

### 8.2 Restoring Balance

Each term gets its own correction:
- $1/\Omega_i$ compensates for $\partial W^i / \partial h_i$
- $1/\Omega_j$ compensates for $\partial W^j / \partial h_j$

After correction:
$$\frac{1}{\Omega_i} \nabla W^i + \frac{1}{\Omega_j} \nabla W^j \approx \text{(what it would be if } h = \text{const)}$$

---

## 9. Does Standard SPH Need Ω?

### 9.1 Strictly Speaking, Yes

The Lagrangian derivation shows that Ω should appear:
$$\mathbf{a}_i = -\sum_j m_j \left( \frac{P_i}{\rho_i^2 \Omega_i} + \frac{P_j}{\rho_j^2 \Omega_j} \right) \nabla W_{ij}$$

### 9.2 But the Error is Smaller

Because standard SPH uses a **single shared kernel gradient**, the h-variation error:
1. Is symmetric between particles
2. Partially cancels in the sum over neighbors
3. Mainly affects force magnitude, not direction

The remaining error is often small enough to be masked by other SPH errors (kernel bias, tensile instability, etc.).

### 9.3 When Standard SPH Fails

Standard SPH without Ω can still fail in:
- Very steep density gradients
- Long-time integrations where errors accumulate
- High-precision requirements

But for typical astrophysical simulations, the built-in error cancellation makes it "good enough."

---

## 10. Summary

### The Core Insight

| Property | Standard SPH | GSPH |
|----------|-------------|------|
| Kernel gradient | Single shared $\nabla W_{ij}$ | Separate $\nabla W^i$, $\nabla W^j$ |
| h-variation error | Symmetric, partially cancels | Asymmetric, accumulates |
| Force direction | Error ∥ physical force | Error has spurious component |
| Needs Ω for stability? | Usually no | **Yes** |

### Why GSPH is More Sensitive

GSPH is a **higher-order method** that:
1. Captures more physics (Riemann solver)
2. Uses more information (kernel convolution)
3. Is **more sensitive to inconsistencies**

The h-variation inconsistency that standard SPH can tolerate becomes a fatal flaw in GSPH.

### The Price of Accuracy

GSPH's greater accuracy comes with a requirement for **greater consistency**. The Ω correction provides this consistency by properly accounting for variable resolution in the Lagrangian formulation.

---

## 11. Analogy: Numerical Differentiation

### First-Order Difference (Standard SPH)
$$f'(x) \approx \frac{f(x+h) - f(x-h)}{2h}$$

Errors are symmetric: $+O(h^2)$ on both sides, they cancel.

### Higher-Order Method (GSPH)
$$f'(x) \approx \frac{-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)}{12h}$$

More accurate, but if $h$ varies, the coefficients (8, -1) no longer match up correctly. **The method becomes unstable unless you correct for the varying h.**

This is exactly what happens with GSPH: higher-order accuracy requires consistent treatment of variable resolution.

---

## References

1. Springel, V. & Hernquist, L. (2002). MNRAS, 333, 649.
2. Inutsuka, S. (2002). JCP, 179, 238.
3. Hopkins, P.F. (2013). MNRAS, 428, 2840.
