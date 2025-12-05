# Lagrangian GSPH Formulation for Cloud Scattering by a Point Mass Black Hole

## Document Information
- **Created:** 2025-12-04
- **Purpose:** Detailed mathematical formulation with explicit coordinate notation
- **Related Code:** `src/gsph/g_fluid_force.cpp`, `src/external_forces/point_mass_bh.cpp`

---

## 1. Coordinate System and Notation

### 1.1 Reference Frames

We use an **inertial (Eulerian) coordinate system** with origin at some fixed point in space. All positions are measured in this global frame.

| Symbol | Description | Type |
|--------|-------------|------|
| $\mathbf{x}_i(t)$ | Position of SPH particle $i$ at time $t$ | Vector (Lagrangian trajectory) |
| $\mathbf{x}_{BH}(t)$ | Position of black hole at time $t$ | Vector (may be fixed or moving) |
| $\mathbf{x}_j(t)$ | Position of neighboring particle $j$ | Vector |

### 1.2 Key Distinction: Lagrangian vs Eulerian

**Eulerian description:** Fields defined at fixed spatial points $\mathbf{x}$
$$\rho(\mathbf{x}, t), \quad P(\mathbf{x}, t), \quad \mathbf{v}(\mathbf{x}, t)$$

**Lagrangian description:** Following individual fluid elements (particles)
$$\rho_i(t) = \rho(\mathbf{x}_i(t), t), \quad P_i(t), \quad \mathbf{v}_i(t) = \frac{d\mathbf{x}_i}{dt}$$

**SPH is fundamentally Lagrangian:** We track particles as they move through space.

### 1.3 Relative Vectors

| Symbol | Definition | Description |
|--------|------------|-------------|
| $\mathbf{r}_{ij}$ | $\mathbf{x}_i - \mathbf{x}_j$ | Vector from particle $j$ to particle $i$ |
| $r_{ij}$ | $\|\mathbf{r}_{ij}\|$ | Distance between particles $i$ and $j$ |
| $\hat{\mathbf{e}}_{ij}$ | $\mathbf{r}_{ij} / r_{ij}$ | Unit vector from $j$ toward $i$ |
| $\mathbf{R}_i$ | $\mathbf{x}_i - \mathbf{x}_{BH}$ | Vector from BH to particle $i$ |
| $R_i$ | $\|\mathbf{R}_i\|$ | Distance from particle $i$ to BH |

```
                    Global Inertial Frame
                    
        y ↑
          │
          │         ● Particle i at x_i(t)
          │        /
          │       / R_i = x_i - x_BH
          │      /
          │     ◉ Black Hole at x_BH(t)
          │      \
          │       \ R_j = x_j - x_BH  
          │        \
          │         ● Particle j at x_j(t)
          │
          └─────────────────────────→ x
          Origin (fixed in space)
          
        Between particles:
        r_ij = x_i - x_j
        ê_ij = r_ij / |r_ij|
```

---

## 2. Governing Equations in Lagrangian Form

### 2.1 Lagrangian (Material) Derivative

For any field quantity $Q$, the Lagrangian derivative following particle $i$ is:

$$\frac{DQ_i}{Dt} \equiv \frac{d}{dt}Q(\mathbf{x}_i(t), t) = \frac{\partial Q}{\partial t}\bigg|_{\mathbf{x}_i} + \mathbf{v}_i \cdot \nabla Q\bigg|_{\mathbf{x}_i}$$

### 2.2 Equations of Motion

#### Trajectory Equation:
$$\frac{d\mathbf{x}_i}{dt} = \mathbf{v}_i$$

#### Momentum Equation:
$$\frac{d\mathbf{v}_i}{dt} = -\frac{1}{\rho_i}(\nabla P)_i + \mathbf{g}_{self,i} + \mathbf{g}_{BH,i}$$

#### Energy Equation:
$$\frac{du_i}{dt} = -\frac{P_i}{\rho_i}(\nabla \cdot \mathbf{v})_i$$

#### Equation of State:
$$P_i = (\gamma - 1)\rho_i u_i$$

---

## 3. Gravity Terms with Explicit Coordinates

### 3.1 Self-Gravity (Cloud Particles on Each Other)

Each particle $i$ feels gravitational attraction from all other cloud particles $j$:

$$\mathbf{g}_{self,i} = -G \sum_{j \neq i} m_j \frac{\mathbf{x}_i - \mathbf{x}_j}{|\mathbf{x}_i - \mathbf{x}_j|^3}$$

**In terms of relative vector $\mathbf{r}_{ij}$:**
$$\mathbf{g}_{self,i} = -G \sum_{j \neq i} m_j \frac{\mathbf{r}_{ij}}{r_{ij}^3}$$

**Properties:**
- Varies significantly across the cloud (particles at center vs edge feel different forces)
- Creates the internal hydrostatic structure: $\nabla P = \rho \mathbf{g}_{self}$
- Computed via Barnes-Hut tree in code (`bhtree.cpp`)

### 3.2 External Black Hole Gravity

Each particle $i$ feels gravitational attraction from the point mass BH:

$$\mathbf{g}_{BH,i} = -\frac{G M_{BH}}{(|\mathbf{x}_i - \mathbf{x}_{BH}|^2 + \epsilon^2)^{3/2}} (\mathbf{x}_i - \mathbf{x}_{BH})$$

**In terms of relative vector $\mathbf{R}_i$:**
$$\mathbf{g}_{BH,i} = -\frac{G M_{BH}}{(R_i^2 + \epsilon^2)^{3/2}} \mathbf{R}_i$$

where $\epsilon$ is a softening length to prevent singularity.

**Properties:**
- All particles at similar distance from BH feel nearly the same acceleration
- Does NOT create internal pressure gradient in the cloud
- Computed in `point_mass_bh.cpp`

### 3.3 Why They Enter Equations Differently

Consider two neighboring particles $i$ and $j$ in the cloud:

**Self-gravity differs significantly:**
$$\mathbf{g}_{self,i} \neq \mathbf{g}_{self,j}$$
because they are at different positions within the cloud's own gravitational field.

**BH gravity is nearly identical (to first order):**
$$\mathbf{g}_{BH,j} = \mathbf{g}_{BH,i} + \underbrace{\mathcal{O}\left(\frac{r_{ij}}{R_i}\right)}_{\text{tidal correction}}$$

Since $r_{ij} \ll R_i$ (particle separation ≪ distance to BH), both particles experience essentially the same BH acceleration.

---

## 4. GSPH Discretization with Explicit Coordinates

### 4.1 SPH Density

$$\rho_i = \sum_j m_j W(|\mathbf{x}_i - \mathbf{x}_j|, h_i) = \sum_j m_j W(r_{ij}, h_i)$$

### 4.2 The Riemann Problem Interface

For particle pair $(i, j)$, the Riemann problem is solved at their **interface** (midpoint):

$$\mathbf{x}_{interface} = \frac{\mathbf{x}_i + \mathbf{x}_j}{2}$$

```
    Particle i                Interface               Particle j
    at x_i                    at x_int                at x_j
        ●─────────────────────────*─────────────────────────●
        |←────── r_ij/2 ─────────→|←─────── r_ij/2 ────────→|
        
    Properties:              Riemann problem          Properties:
    ρ_i, P_i, v_i            solved here             ρ_j, P_j, v_j
    g_self,i                 → P*, v*                g_self,j
```

**1D projection along pair axis:**

Normal velocities:
$$v_i^n = \mathbf{v}_i \cdot \hat{\mathbf{e}}_{ij}, \quad v_j^n = \mathbf{v}_j \cdot \hat{\mathbf{e}}_{ij}$$

Self-gravity projection:
$$g_{self,i}^n = \mathbf{g}_{self,i} \cdot \hat{\mathbf{e}}_{ij}, \quad g_{self,j}^n = \mathbf{g}_{self,j} \cdot \hat{\mathbf{e}}_{ij}$$

### 4.3 Well-Balanced Pressure Extrapolation

**Standard (non-well-balanced) Riemann solver:**
Uses $P_i$ and $P_j$ directly → spurious flow in equilibrium.

**Well-balanced Riemann solver:**
Extrapolate pressures from particle positions to interface, accounting for hydrostatic gradient:

$$P_i^{int} = P_i - \rho_i \, g_{self,i}^n \cdot \frac{r_{ij}}{2}$$

$$P_j^{int} = P_j + \rho_j \, g_{self,j}^n \cdot \frac{r_{ij}}{2}$$

**Physical interpretation:**
- Particle $i$ is at distance $r_{ij}/2$ from interface
- Moving toward interface (in direction $-\hat{\mathbf{e}}_{ij}$), pressure decreases if $g_{self,i}^n > 0$
- The extrapolation reconstructs what pressure "would be" at the interface

**Why only self-gravity?**
- Self-gravity creates the equilibrium: $\nabla P = \rho \mathbf{g}_{self}$
- BH gravity does NOT create internal pressure gradient
- In the local (Lagrangian) frame of particles $i$ and $j$, BH gravity is uniform

### 4.4 Sound Speeds at Interface

$$c_i^{int} = \sqrt{\frac{\gamma P_i^{int}}{\rho_i}}, \quad c_j^{int} = \sqrt{\frac{\gamma P_j^{int}}{\rho_j}}$$

### 4.5 Riemann Solver

**Left state (from particle $j$, toward increasing $\hat{\mathbf{e}}_{ij}$):**
$$(v_L, \rho_L, P_L, c_L) = (v_j^n, \rho_j, P_j^{int}, c_j^{int})$$

**Right state (from particle $i$):**
$$(v_R, \rho_R, P_R, c_R) = (v_i^n, \rho_i, P_i^{int}, c_i^{int})$$

**HLL Solver Output:**
$$v^* = \frac{P_L - P_R + \rho_L v_L (S_L - v_L) - \rho_R v_R (S_R - v_R)}{\rho_L (S_L - v_L) - \rho_R (S_R - v_R)}$$

$$P^* = P_L + \rho_L (S_L - v_L)(v^* - v_L)$$

### 4.6 GSPH Acceleration (Fluid Only)

$$\mathbf{a}_{fluid,i} = -\sum_j m_j P^*_{ij} \left( \frac{\Omega_i}{\rho_i^2} \nabla_i W(r_{ij}, h_i) + \frac{\Omega_j}{\rho_j^2} \nabla_i W(r_{ij}, h_j) \right)$$

where $\Omega_i$ are the **grad-h correction factors** (explained in detail below).

---

## 5. The Grad-h Correction: Why It's Essential

### 5.1 The Problem: Variable Smoothing Length

In SPH, the smoothing length $h_i$ adapts to local particle density:
$$h_i \propto \rho_i^{-1/D}$$
where $D$ is the number of dimensions.

This means $h_i$ is **not constant** - it varies from particle to particle and changes over time. This creates a serious problem for energy and momentum conservation.

### 5.2 The Standard SPH Density Estimate

The SPH density is:
$$\rho_i = \sum_j m_j W(|\mathbf{x}_i - \mathbf{x}_j|, h_i)$$

When we take the **Lagrangian derivative** (following the particle):
$$\frac{D\rho_i}{Dt} = \sum_j m_j \frac{\partial W}{\partial r_{ij}} \frac{Dr_{ij}}{Dt} + \sum_j m_j \frac{\partial W}{\partial h_i} \frac{Dh_i}{Dt}$$

The **second term** arises because $h_i$ changes as the particle moves! If we ignore this term, we get **inconsistent** derivatives that violate conservation laws.

### 5.3 Derivation of the Grad-h Correction Factor

The smoothing length is typically set by requiring a constant number of neighbors:
$$\frac{4\pi}{3} h_i^D \bar{n} = N_{neighbor} = \text{const}$$

where $\bar{n} = \rho/m$ is the number density. This gives:
$$h_i \propto \rho_i^{-1/D}$$

Taking the derivative:
$$\frac{\partial h_i}{\partial \rho_i} = -\frac{h_i}{D \rho_i}$$

Now, the density itself depends on $h_i$:
$$\rho_i(h_i) = \sum_j m_j W(r_{ij}, h_i)$$

So we have an **implicit equation**. Using the chain rule:
$$\frac{d\rho_i}{dh_i} = \frac{\partial \rho_i}{\partial h_i}\bigg|_{explicit} + \frac{\partial \rho_i}{\partial h_i}\bigg|_{implicit}$$

The explicit derivative is:
$$\frac{\partial \rho_i}{\partial h_i} = \sum_j m_j \frac{\partial W}{\partial h}(r_{ij}, h_i)$$

The **grad-h correction factor** $\Omega_i$ accounts for this:

$$\boxed{\Omega_i = \frac{1}{1 + \frac{h_i}{D\rho_i} \frac{\partial \rho_i}{\partial h_i}} = \frac{1}{1 + \frac{h_i}{D\rho_i} \sum_j m_j \frac{\partial W}{\partial h}(r_{ij}, h_i)}}$$

### 5.4 Why Grad-h Prevents Core Collapse

**Without grad-h correction ($\Omega_i = 1$):**

Consider a self-gravitating cloud in equilibrium. As numerical noise causes slight compression:

1. **Density increases** in the core: $\rho_i \uparrow$
2. **Smoothing length decreases**: $h_i \downarrow$ (to maintain neighbor count)
3. **Kernel becomes sharper**: $\nabla W$ increases in magnitude
4. **Spurious pressure gradient** appears because the kernel derivative changed, NOT because the physical pressure gradient changed
5. **This spurious gradient is WRONG** - it doesn't correctly represent $\nabla P$
6. The error typically **underestimates** the pressure support
7. **Core collapses** due to insufficient pressure support

**With grad-h correction ($\Omega_i$ computed correctly):**

The correction factor compensates for the changing kernel:

1. When $h_i$ decreases (compression), $\frac{\partial \rho_i}{\partial h_i} < 0$
2. This makes $\Omega_i > 1$
3. The **amplified** kernel gradient correctly represents the true pressure gradient
4. Pressure support is maintained → **equilibrium preserved**

### 5.5 Mathematical Interpretation: Consistent Derivatives

The grad-h correction ensures that the **discrete SPH operators are consistent** with the continuum equations.

**Standard SPH pressure gradient:**
$$-\frac{1}{\rho_i}(\nabla P)_i \approx -\sum_j m_j \left(\frac{P_i}{\rho_i^2} + \frac{P_j}{\rho_j^2}\right) \nabla_i W_{ij}$$

**With grad-h correction:**
$$-\frac{1}{\rho_i}(\nabla P)_i \approx -\sum_j m_j \left(\frac{P_i}{\Omega_i \rho_i^2} + \frac{P_j}{\Omega_j \rho_j^2}\right) \nabla_i W_{ij}$$

Wait - note that in our GSPH formulation, we write it as:
$$\mathbf{a}_{fluid,i} = -\sum_j m_j P^*_{ij} \left( \frac{\Omega_i}{\rho_i^2} \nabla_i W_{ij}^{(i)} + \frac{\Omega_j}{\rho_j^2} \nabla_i W_{ij}^{(j)} \right)$$

The position of $\Omega$ (numerator vs denominator) depends on the specific SPH formulation. In our code, we use the Springel & Hernquist (2002) convention where $\Omega$ appears in the **numerator**.

### 5.6 The Grad-h Term in Code

From `g_pre_interaction.cpp`:
```cpp
// Compute kernel derivative w.r.t. h
dh_dens_i += p_j.mass * kernel->dhw(r, p_i.sml);

// Grad-h correction factor: Ω_i = 1 / (1 + (h/Dρ) * dρ/dh)
p_i.gradh = 1.0 / (1.0 + p_i.sml / (DIM * dens_i) * dh_dens_i);
```

From `g_fluid_force.cpp`:
```cpp
// Grad-h correction (Springel & Hernquist 2002)
const real omega_i = p_i.gradh;
const real omega_j = p_j.gradh;

// Standard GSPH force with pstar from well-balanced Riemann solver
const vec_t f = dw_i * (p_j.mass * pstar * rho2_inv_i * omega_i) 
              + dw_j * (p_j.mass * pstar * rho2_inv_j * omega_j);
```

### 5.7 Summary: Two Essential Corrections for Equilibrium

| Correction | Problem Solved | Without It |
|------------|----------------|------------|
| **Well-balanced Riemann** | Pressure difference misinterpreted as discontinuity | Spurious flow → collapse |
| **Grad-h correction** | Variable $h$ causes inconsistent derivatives | Pressure support error → collapse |

**Both corrections are necessary** for maintaining self-gravitating equilibrium:

1. **Well-balanced Riemann**: Handles the physics (gravity creates pressure gradient)
2. **Grad-h**: Handles the numerics (variable resolution creates derivative errors)

```
                    Stable Equilibrium Requires Both
                    
    ┌─────────────────────────────────────────────────────┐
    │                                                     │
    │   Physical: ∇P = ρg     →  Well-balanced Riemann   │
    │                                                     │
    │   Numerical: h varies   →  Grad-h correction       │
    │                                                     │
    └─────────────────────────────────────────────────────┘
                              ↓
                    Cloud remains stable
```

---

## 6. Total Acceleration: Putting It All Together

$$\frac{d\mathbf{v}_i}{dt} = \underbrace{\mathbf{a}_{fluid,i}}_{\substack{\text{GSPH with} \\ \text{well-balanced } P^*}} + \underbrace{\mathbf{g}_{self,i}}_{\substack{\text{Self-gravity} \\ \text{(in Riemann \& direct)}}} + \underbrace{\mathbf{g}_{BH,i}}_{\substack{\text{External BH} \\ \text{(direct only)}}}$$

**Explicitly:**

$$\frac{d\mathbf{v}_i}{dt} = -\sum_j m_j P^*_{ij}(\mathbf{g}_{self}) \left( \frac{\Omega_i}{\rho_i^2} \nabla_i W_{ij}^{(i)} + \frac{\Omega_j}{\rho_j^2} \nabla_i W_{ij}^{(j)} \right)$$
$$- G \sum_{j \neq i} m_j \frac{\mathbf{x}_i - \mathbf{x}_j}{|\mathbf{x}_i - \mathbf{x}_j|^3}$$
$$- \frac{G M_{BH}}{(|\mathbf{x}_i - \mathbf{x}_{BH}|^2 + \epsilon^2)^{3/2}} (\mathbf{x}_i - \mathbf{x}_{BH})$$

---

## 7. Why External BH Force is NOT in Well-Balanced Correction

### 7.1 Physical Argument: Local Comoving Frame

Consider a small region of the cloud containing particles $i$ and $j$. Transform to a **local freely-falling frame** accelerating at $\mathbf{g}_{BH}$:

$$\mathbf{x}'_i = \mathbf{x}_i - \mathbf{x}_{cm}, \quad \mathbf{v}'_i = \mathbf{v}_i - \mathbf{v}_{cm}$$

where the center of mass follows:
$$\frac{d^2 \mathbf{x}_{cm}}{dt^2} = \mathbf{g}_{BH}(\mathbf{x}_{cm})$$

In this frame:
- **BH gravity vanishes** (to first order) - equivalence principle!
- **Self-gravity remains** - it's internal to the cloud
- Particles $i$ and $j$ interact via fluid forces and self-gravity only

The Riemann problem is naturally formulated in this local comoving frame.

### 7.2 Mathematical Argument: Taylor Expansion

For neighboring particles separated by $\mathbf{r}_{ij}$:

$$\mathbf{g}_{BH,j} = \mathbf{g}_{BH,i} + (\mathbf{r}_{ij} \cdot \nabla)\mathbf{g}_{BH}\big|_{\mathbf{x}_i} + \mathcal{O}(r_{ij}^2)$$

The tidal tensor is:
$$\nabla \mathbf{g}_{BH} = -\frac{GM_{BH}}{R^3}\left(\mathbf{I} - 3\hat{\mathbf{R}}\hat{\mathbf{R}}\right)$$

The difference in BH acceleration between neighbors:
$$|\mathbf{g}_{BH,i} - \mathbf{g}_{BH,j}| \sim \frac{GM_{BH}}{R^3} r_{ij}$$

**Ratio to self-gravity:**
$$\frac{|\Delta \mathbf{g}_{BH}|}{|\mathbf{g}_{self}|} \sim \frac{M_{BH}/R^3}{M_{cloud}/R_{cloud}^3} \cdot \frac{r_{ij}}{R_{cloud}} \ll 1$$

The tidal correction is negligible at the particle-pair scale.

### 7.3 Equilibrium Structure Argument

The cloud's internal pressure profile is determined by:
$$\nabla P = \rho \mathbf{g}_{self}$$

NOT by:
$$\nabla P = \rho (\mathbf{g}_{self} + \mathbf{g}_{BH}) \quad \text{(WRONG for cloud structure)}$$

The BH doesn't create internal pressure gradients - it just accelerates the whole cloud.

---

## 8. Algorithm Summary

```
For each timestep:

1. UPDATE POSITIONS (from previous velocities)
   x_i(t+dt) = x_i(t) + v_i * dt + 0.5 * a_i * dt²   (leapfrog/predictor)

2. BUILD TREE & FIND NEIGHBORS
   Using current positions x_i(t+dt)

3. COMPUTE DENSITIES
   ρ_i = Σ_j m_j W(|x_i - x_j|, h_i)

4. COMPUTE SELF-GRAVITY → store in g_self,i
   g_self,i = -G Σ_j m_j (x_i - x_j) / |x_i - x_j|³
   
5. COMPUTE GSPH FLUID FORCE (uses g_self in well-balanced correction)
   For each pair (i,j):
     r_ij = x_i - x_j
     ê_ij = r_ij / |r_ij|
     
     # Project self-gravity (NOT BH gravity!)
     g_i^n = g_self,i · ê_ij
     g_j^n = g_self,j · ê_ij
     
     # Well-balanced pressure extrapolation
     P_i^int = P_i - ρ_i * g_i^n * |r_ij|/2
     P_j^int = P_j + ρ_j * g_j^n * |r_ij|/2
     
     # Solve Riemann problem
     (P*, v*) = RiemannSolver(left_state, right_state)
     
     # Accumulate GSPH acceleration
     a_fluid,i += ...

6. ADD SELF-GRAVITY TO ACCELERATION
   a_i = a_fluid,i + g_self,i

7. ADD EXTERNAL BH GRAVITY (direct, no Riemann modification)
   R_i = x_i - x_BH
   g_BH,i = -G M_BH * R_i / (|R_i|² + ε²)^(3/2)
   a_i += g_BH,i

8. UPDATE BH POSITION (if moving)
   x_BH(t+dt) = x_BH(t) + v_BH * dt

9. UPDATE VELOCITIES
   v_i(t+dt) = v_i(t) + a_i * dt   (corrector)

10. COMPUTE ENERGY CHANGE
    du_i/dt = ... (from GSPH)
```

---

## 9. Coordinate Summary Table

| Quantity | Symbol | Definition | Frame |
|----------|--------|------------|-------|
| Particle position | $\mathbf{x}_i$ | Absolute position | Inertial |
| BH position | $\mathbf{x}_{BH}$ | Absolute position | Inertial |
| Particle-particle separation | $\mathbf{r}_{ij} = \mathbf{x}_i - \mathbf{x}_j$ | Relative | - |
| Particle-BH separation | $\mathbf{R}_i = \mathbf{x}_i - \mathbf{x}_{BH}$ | Relative | - |
| Interface position | $\mathbf{x}_{int} = (\mathbf{x}_i + \mathbf{x}_j)/2$ | Absolute | Inertial |
| Self-gravity at $i$ | $\mathbf{g}_{self,i}$ | Computed at $\mathbf{x}_i$ | Inertial |
| BH gravity at $i$ | $\mathbf{g}_{BH,i}$ | Computed at $\mathbf{x}_i$ | Inertial |
| Riemann velocity | $v^*$ | Along $\hat{\mathbf{e}}_{ij}$ | Local 1D |
| Riemann pressure | $P^*$ | At interface | Local 1D |

---

## 10. Code Mapping

| Physics | Code Location | Variable |
|---------|---------------|----------|
| Particle position $\mathbf{x}_i$ | `particle.hpp` | `p.pos` |
| BH position $\mathbf{x}_{BH}$ | `point_mass_bh.cpp` | `m_position` |
| Self-gravity $\mathbf{g}_{self,i}$ | `gravity_force.cpp`, `bhtree.cpp` | `p.grav_acc` |
| BH gravity $\mathbf{g}_{BH,i}$ | `point_mass_bh.cpp` | Added directly to `p.acc` |
| Well-balanced extrapolation | `g_fluid_force.cpp:100-160` | `p_i_interface`, `p_j_interface` |
| Riemann solver | `g_fluid_force.cpp:170+` | `pstar`, `vstar` |
| Grad-h correction | `g_pre_interaction.cpp:109` | `p.gradh` |

---

## 11. References

1. Inutsuka (2002) - GSPH formulation
2. Hopkins (2015) - GIZMO and well-balanced schemes
3. Springel & Hernquist (2002) - Grad-h corrections
4. Katz et al. (1996) - TreeSPH with self-gravity
