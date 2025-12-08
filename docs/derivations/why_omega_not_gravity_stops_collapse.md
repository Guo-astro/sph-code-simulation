# Why the Ω Correction in Hydrodynamics (Not Gravity) Stops Core Collapse

## Your Observation

In your simulations, you found:
- Gravity kernel convolution (Wendland C4 softening) → **still collapses**
- Ω correction in hydrodynamics → **stops collapse**

This is a deep result. Let me explain why from first principles.

---

## 1. The Fundamental Asymmetry

### 1.1 What Determines Hydrostatic Equilibrium?

In equilibrium:
$$\nabla P = \rho \mathbf{g}$$

Rewriting:
$$\frac{\nabla P}{\rho} = \mathbf{g}$$

The **specific force** (force per unit mass) must balance on each fluid element.

### 1.2 The Key Insight

In SPH, we compute:
- **Pressure acceleration:** $\mathbf{a}^P_i = -\sum_j m_j \frac{P}{\rho^2} \nabla W$
- **Gravitational acceleration:** $\mathbf{a}^g_i = -\sum_j G m_j g(r,h) \hat{\mathbf{r}}$

For equilibrium: $\mathbf{a}^P_i + \mathbf{a}^g_i = 0$

**The error in equilibrium depends on which acceleration is wrong!**

---

## 2. Where Does the h-Variation Error Enter?

### 2.1 Pressure Force Error

The pressure acceleration involves:
$$\mathbf{a}^P_i = -\sum_j m_j \frac{P_j}{\rho_j^2} \nabla W(r_{ij}, h_j)$$

When $h = h(\rho)$, the kernel gradient is:
$$\nabla W(r, h(\rho)) = \frac{\partial W}{\partial r} \hat{\mathbf{r}} + \frac{\partial W}{\partial h} \nabla h$$

The second term is **spurious** — it comes from resolution variation, not physics!

**Error in pressure acceleration:**
$$\delta \mathbf{a}^P \sim \sum_j m_j \frac{P}{\rho^2} \frac{\partial W}{\partial h} \nabla h \sim \frac{P}{\rho} \frac{\nabla h}{h}$$

### 2.2 Gravitational Force Error

The gravitational acceleration involves:
$$\mathbf{a}^g_i = -\sum_j G m_j g(r_{ij}, h_{ij}) \hat{\mathbf{r}}_{ij}$$

When $h$ varies, we get:
$$\delta \mathbf{a}^g \sim \sum_j G m_j \frac{\partial g}{\partial h} \delta h$$

**But this is fundamentally different!**

---

## 3. The Critical Difference: Density Weighting

### 3.1 Pressure Force Structure

The pressure term has the structure:
$$\mathbf{a}^P \sim \frac{P}{\rho^2} \times (\text{kernel gradient})$$

The $1/\rho^2$ factor means this is a **density-weighted integral**. The kernel gradient $\nabla W$ is being used to **estimate a continuum derivative** $\nabla P / \rho$.

When $h$ varies, the kernel gradient **misrepresents the physical gradient** because:
$$\nabla W(r, h(x)) \neq \nabla W(r, h)|_{h=\text{const}}$$

The Ω correction fixes this by accounting for the implicit h-dependence in the density sum.

### 3.2 Gravitational Force Structure

The gravity term has the structure:
$$\mathbf{a}^g \sim G m_j \times g(r, h) \times \hat{\mathbf{r}}$$

This is a **direct pairwise sum** — it's not trying to estimate a continuum gradient!

The softening kernel $g(r, h)$ just determines **how much force** each particle pair contributes. Changing $h$ changes the softening, but:
- The force is still in the correct direction ($\hat{\mathbf{r}}$)
- The total enclosed mass is still correct
- There's no "gradient estimation" that can be corrupted

---

## 4. Mathematical Proof

### 4.1 The Density Gradient Problem

In hydrostatic equilibrium with a polytrope:
$$P = K \rho^\gamma, \quad \nabla P = K \gamma \rho^{\gamma-1} \nabla \rho$$

The SPH pressure gradient is:
$$(\nabla P)_{\text{SPH}} = \rho_i \sum_j m_j \frac{P_j}{\rho_j^2} \nabla W_{ij}$$

Using the SPH identity:
$$\nabla \rho = \sum_j m_j \nabla W_{ij}$$

We need:
$$\sum_j m_j \frac{P_j}{\rho_j^2} \nabla W_{ij} \stackrel{?}{=} \frac{\nabla P}{\rho} = \frac{K \gamma \rho^{\gamma-1} \nabla \rho}{\rho}$$

### 4.2 The Error Without Ω

When $h = h(\rho)$:
$$\nabla W_{ij} = \nabla W|_h + \frac{\partial W}{\partial h} \nabla h$$

The spurious term is:
$$\sum_j m_j \frac{P_j}{\rho_j^2} \frac{\partial W}{\partial h} \nabla h$$

This **does not vanish** even in equilibrium because:
- $\nabla h \neq 0$ (h varies with density)
- $\partial W / \partial h \neq 0$ (kernel depends on h)
- The sum is weighted by $P/\rho^2$, which varies spatially

### 4.3 Why Gravity Doesn't Have This Problem

The gravitational force is:
$$\mathbf{F}^g_i = -\sum_j G m_i m_j g(r_{ij}, h_{ij}) \hat{\mathbf{r}}_{ij}$$

There is **no gradient operator** acting on a field here! The force kernel $g(r, h)$ is just a function of the pair separation.

Even if $h$ varies:
- The force direction is still $\hat{\mathbf{r}}_{ij}$ (correct)
- The magnitude changes smoothly with $h$
- There's no "spurious gradient" because we're not estimating $\nabla \Phi$

**The gravitational potential $\Phi$ and its gradient are computed directly from the mass distribution, not from a kernel-convolved field.**

---

## 5. Physical Intuition

### 5.1 What SPH Actually Computes

**Pressure force:** "What is the pressure gradient here?"
- Uses kernel to estimate a local derivative
- h-variation corrupts this estimate
- Needs Ω correction

**Gravitational force:** "How hard is that mass pulling on me?"
- Direct pairwise calculation
- Softening just prevents singularity
- No derivative estimation to corrupt

### 5.2 Analogy: Measuring vs Computing

Imagine measuring temperature gradient with thermometers:
- If your thermometer spacing varies, you need to correct for that
- This is like the Ω correction for pressure gradients

Now imagine measuring gravitational pull with a spring scale:
- You just measure the force directly
- Spacing of measurements doesn't affect the reading
- This is like gravity — no correction needed

---

## 6. The Core Collapse Mechanism

### 6.1 Without Ω Correction

In the core:
1. $\rho$ is high → $h$ is small
2. Moving outward: $\partial h / \partial r > 0$
3. Spurious term: $\frac{\partial W}{\partial h} \nabla h$ points **outward**
4. This **reduces** the effective pressure gradient
5. $|\nabla P|_{\text{SPH}} < |\nabla P|_{\text{true}}$
6. Gravity wins → core contracts
7. Contraction increases $\rho$ → smaller $h$ → larger error
8. **Runaway collapse**

### 6.2 With Ω Correction

The correction factor:
$$\Omega = 1 - \frac{\partial h}{\partial \rho} \sum_j m_j \frac{\partial W}{\partial h}$$

In the core: $\Omega < 1$

The corrected pressure force:
$$\frac{P}{\rho^2 \Omega} > \frac{P}{\rho^2}$$

This **enhances** the pressure force, exactly compensating for the spurious reduction.

### 6.3 Why Gravity Softening Doesn't Help

Better gravity softening (Wendland C4 vs Hernquist-Katz):
- Makes gravity **smoother** at small r
- Doesn't change the **direction** or **scaling** of gravity
- The gravity is already correct — it's not the problem!

The collapse happens because **pressure is too weak**, not because gravity is too strong.

---

## 7. Quantitative Analysis

### 7.1 Force Balance Error

Define the equilibrium error:
$$\epsilon = \frac{|\mathbf{a}^P + \mathbf{a}^g|}{|\mathbf{a}^g|}$$

**Without Ω:**
$$\epsilon \sim \frac{|\delta \mathbf{a}^P|}{|\mathbf{a}^g|} \sim \frac{h}{R}$$

For typical resolution $h/R \sim 0.05$ to $0.1$: **5-10% force imbalance**

**With Ω:**
$$\epsilon \sim O(h^2 / R^2) \sim 0.25\% \text{ to } 1\%$$

---

### 7.2 Detailed Derivation of the Error Scaling

#### Step 1: The Spurious Pressure Acceleration

The pressure acceleration in SPH is:
$$\mathbf{a}^P_i = -\sum_j m_j \frac{P_j}{\rho_j^2} \nabla_i W(r_{ij}, h_j)$$

When $h = h(\rho)$, the kernel gradient has two parts:
$$\nabla_i W(r_{ij}, h_j) = \underbrace{\frac{\partial W}{\partial r} \hat{\mathbf{r}}_{ij}}_{\text{physical}} + \underbrace{\frac{\partial W}{\partial h} \nabla_i h_j}_{\text{spurious}}$$

The spurious acceleration is:
$$\delta \mathbf{a}^P_i = -\sum_j m_j \frac{P_j}{\rho_j^2} \frac{\partial W}{\partial h} \nabla h_j$$

#### Step 2: Estimate $\nabla h$

From the resolution constraint $h = \eta (m/\rho)^{1/d}$ in $d$ dimensions:
$$\nabla h = \frac{\partial h}{\partial \rho} \nabla \rho = -\frac{h}{d\rho} \nabla \rho$$

Therefore:
$$|\nabla h| \sim \frac{h}{\rho} |\nabla \rho|$$

#### Step 3: Estimate $\partial W / \partial h$

For a normalized kernel $W(r, h) = h^{-d} w(r/h)$:
$$\frac{\partial W}{\partial h} = -\frac{d}{h} W - \frac{r}{h^2} W'$$

Near the kernel center where $r \sim h$ and $W \sim h^{-d}$:
$$\left| \frac{\partial W}{\partial h} \right| \sim \frac{W}{h} \sim \frac{1}{h^{d+1}}$$

#### Step 4: Estimate the Spurious Acceleration

$$|\delta \mathbf{a}^P| \sim \sum_j m_j \frac{P}{\rho^2} \cdot \frac{1}{h^{d+1}} \cdot \frac{h}{\rho} |\nabla \rho|$$

The sum over neighbors gives approximately $\rho$ (since $\sum_j m_j W \approx \rho$):
$$\sim \rho \cdot \frac{1}{h^d} \cdot \frac{P}{\rho^2} \cdot \frac{1}{h^{d+1}} \cdot \frac{h}{\rho} |\nabla \rho| \cdot h^d$$

Simplifying (the $h^d$ factors from summation and kernel normalization cancel):
$$|\delta \mathbf{a}^P| \sim \frac{P}{\rho^2} \cdot \frac{1}{h} \cdot \frac{h}{\rho} |\nabla \rho| = \frac{P}{\rho^3} |\nabla \rho|$$

#### Step 5: Estimate the Gravitational Acceleration

For a self-gravitating sphere of mass $M$ and radius $R$:
$$|\mathbf{a}^g| \sim \frac{GM}{R^2} \sim G\rho R$$

using $M \sim \rho R^3$.

#### Step 6: Estimate the Density Gradient

For a polytropic sphere (Lane-Emden solution), the density varies from center to edge over scale $R$:
$$|\nabla \rho| \sim \frac{\rho}{R}$$

#### Step 7: Estimate the Pressure

In hydrostatic equilibrium:
$$\nabla P \sim \rho g \sim \rho \cdot G\rho R = G\rho^2 R$$

Therefore:
$$P \sim G\rho^2 R^2$$

(This is just the virial estimate: $P \sim GM\rho/R \sim G\rho^2 R^2$)

#### Step 8: Compute the Error Ratio via Dimensional Analysis

The cleanest approach is to compare the spurious and physical parts of the kernel gradient directly.

The physical pressure acceleration scales as:
$$|\mathbf{a}^P| \sim \frac{|\nabla P|}{\rho} \sim \frac{P}{\rho R}$$

The spurious part comes from the $\partial W / \partial h \cdot \nabla h$ term. The ratio is:
$$\frac{|\delta \mathbf{a}^P|}{|\mathbf{a}^P|} \sim \frac{|\partial W / \partial h| \cdot |\nabla h|}{|\partial W / \partial r|}$$

Now:
- $|\partial W / \partial r| \sim W/h \sim 1/h^{d+1}$
- $|\partial W / \partial h| \sim W/h \sim 1/h^{d+1}$
- $|\nabla h| \sim (h/\rho) |\nabla \rho| \sim (h/\rho) \cdot (\rho/R) = h/R$

Therefore:
$$\frac{|\delta \mathbf{a}^P|}{|\mathbf{a}^P|} \sim \frac{(1/h^{d+1}) \cdot (h/R)}{1/h^{d+1}} = \frac{h}{R}$$

Since in equilibrium $|\mathbf{a}^P| = |\mathbf{a}^g|$:
$$\boxed{\epsilon = \frac{|\delta \mathbf{a}^P|}{|\mathbf{a}^g|} \sim \frac{h}{R}}$$

#### Step 9: Physical Interpretation

The error $\epsilon \sim h/R$ has a beautiful interpretation:

- $h$ = resolution scale (where kernel varies)
- $R$ = system scale (where physical quantities vary)
- $h/R$ = ratio of numerical to physical scales

**The error is the ratio of the resolution scale to the system scale.**

This is the fundamental limit of any gradient-based method with variable resolution.

#### Step 10: Numerical Values

For typical SPH simulations:
- $N \sim 10^4$ to $10^6$ particles
- $h/R \sim N^{-1/3}$ (in 3D, h scales with mean particle spacing)
- $N = 10^4$: $h/R \sim 0.046$ → **4.6% error**
- $N = 10^5$: $h/R \sim 0.021$ → **2.1% error**
- $N = 10^6$: $h/R \sim 0.01$ → **1% error**

But even 1% systematic error accumulates over dynamical times!

---

### 7.3 Why Ω Correction Reduces Error to $O(h^2/R^2)$

#### The Ω Factor Definition

$$\Omega = 1 - \frac{\partial h}{\partial \rho} \sum_j m_j \frac{\partial W_{ij}}{\partial h}$$

#### What Ω Does

The corrected pressure acceleration is:
$$\mathbf{a}^P_{\text{corrected}} = -\sum_j m_j \frac{P_j}{\rho_j^2 \Omega_j} \nabla W_{ij}$$

The factor $1/\Omega$ **exactly cancels** the leading-order spurious term from $\partial W / \partial h \cdot \nabla h$.

#### Proof of Cancellation

Expand $\Omega^{-1}$ to first order:
$$\frac{1}{\Omega} = \frac{1}{1 - \frac{\partial h}{\partial \rho} \sum_k m_k \frac{\partial W}{\partial h}} \approx 1 + \frac{\partial h}{\partial \rho} \sum_k m_k \frac{\partial W}{\partial h} + O\left(\left(\frac{h}{R}\right)^2\right)$$

The corrected acceleration becomes:
$$\mathbf{a}^P_{\text{corrected}} \approx -\sum_j m_j \frac{P_j}{\rho_j^2} \nabla W_{ij} \cdot \left(1 + \frac{\partial h}{\partial \rho} \sum_k m_k \frac{\partial W}{\partial h}\right)$$

The second term generates a correction that **exactly matches** the spurious $\partial W / \partial h \cdot \nabla h$ term, but with opposite sign!

#### Residual Error

After the Ω correction, the remaining error is:
- Second-order in $(h/R)$
- Comes from the $O((h/R)^2)$ terms in the expansion

$$\epsilon_{\text{with } \Omega} \sim \left(\frac{h}{R}\right)^2$$

For $h/R = 0.05$: error $\sim 0.25\%$ instead of $5\%$

**This is a 20× improvement in force balance accuracy!**

---

### 7.4 Collapse Timescale

The core collapse time without Ω:
$$t_{\text{collapse}} \sim \frac{R}{\sqrt{\epsilon \cdot g \cdot R}} \sim \frac{t_{\text{dyn}}}{\sqrt{\epsilon}}$$

For $\epsilon \sim 0.1$: collapse in ~3 dynamical times.

With Ω: $\epsilon \sim 0.01$: stable for ~30+ dynamical times.

---

## 8. Why This Makes Physical Sense

### 8.1 The Virial Theorem

In equilibrium:
$$2K + U = 0$$

where $K$ is kinetic (thermal) energy and $U$ is gravitational potential energy.

The pressure gradient comes from $\nabla K$ (thermal motion).
The gravitational force comes from $\nabla U$ (potential well).

**The Ω correction is about correctly computing the thermal pressure gradient**, which is the fundamental quantity that resists gravity.

### 8.2 Information Content

- **Gravity:** Only needs to know where the mass is (positions)
- **Pressure:** Needs to know the local thermodynamic state AND its gradient

Gradients are more sensitive to discretization errors than integrals. The Ω correction is essentially fixing the gradient estimation.

---

## 9. Summary

### Why Ω Hydrodynamic Correction Works:

1. **Pressure force estimates a gradient** using kernel convolution
2. When $h = h(\rho)$, the gradient estimate is **systematically biased**
3. The bias **reduces** effective pressure in high-density regions
4. Ω correction **exactly compensates** for this bias

### Why Gravity Convolution Doesn't Help:

1. **Gravity is a direct pairwise sum**, not a gradient estimate
2. Softening changes magnitude but not the structure of the force
3. The gravity calculation is **already correct**
4. The problem is pressure being too weak, not gravity being too strong

### The Bottom Line:

$$\boxed{\text{Core collapse} = \text{Pressure underestimated} \neq \text{Gravity overestimated}}$$

The Ω correction fixes the pressure underestimation.
Gravity softening doesn't address the actual problem.

---

## 10. Implications for Your Code

### Current State:
- `g_fluid_force.cpp`: Has Ω correction ✓
- `gravity_force.cpp`: Has Wendland C4 softening ✓

### Why It Works:
The Ω correction in hydrodynamics is the **essential** ingredient. The gravity softening is nice for accuracy but doesn't affect the equilibrium stability.

### Prediction:
You could use **any reasonable** gravity softening (Hernquist-Katz, Wendland, Plummer) and still maintain hydrostatic equilibrium **as long as** the Ω correction is present in the hydrodynamics.

Conversely, even "perfect" gravity softening will not prevent collapse without the Ω correction.

---

## References

1. Springel, V. & Hernquist, L. (2002). MNRAS, 333, 649. — Original grad-h formulation

2. Price, D.J. & Monaghan, J.J. (2007). MNRAS, 374, 1347. — Gravity softening analysis

3. Hopkins, P.F. (2013). MNRAS, 428, 2840. — General Lagrangian SPH methods
