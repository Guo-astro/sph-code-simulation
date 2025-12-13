# First Principles Derivation: Particle Count for Lane-Emden SPH Convergence

## 1. Problem Statement

For a Lane-Emden n=1.5 (γ=5/3) hydrostatic sphere simulated with SPH, we need to determine the minimum number of particles N required to achieve a given accuracy ε in the central density.

**Observed:** With N ≈ 10,648 particles, the central density error is ~16%
**Goal:** Derive N(ε) from first principles

---

## 2. SPH Density Estimation Error

### 2.1 SPH Density Estimator

The SPH density at particle position **r**ᵢ is:

$$\rho_i = \sum_j m_j W(|\mathbf{r}_i - \mathbf{r}_j|, h_i)$$

where W is the kernel function and h is the smoothing length.

### 2.2 Smoothing Length Scaling

For a given kernel with neighbor number $N_{nb}$, the smoothing length scales with the local inter-particle spacing:

$$h \sim \left(\frac{N_{nb} \cdot m}{\rho}\right)^{1/D}$$

where D is the spatial dimension (D=3 for 3D).

For N total particles of equal mass m = M/N in a sphere of radius R with mean density ρ̄:

$$h \sim \left(\frac{M}{N \cdot \bar{\rho}}\right)^{1/3} \sim \frac{R}{N^{1/3}}$$

### 2.3 Density Error from Kernel Truncation

The SPH density estimate has a systematic error due to kernel truncation. For a smooth density field, the leading order error is:

$$\frac{\delta\rho}{\rho} \sim \mathcal{O}\left(\frac{h^2}{L_\rho^2}\right)$$

where $L_\rho = \rho / |\nabla\rho|$ is the local density scale length.

---

## 3. Lane-Emden Specific Analysis

### 3.1 Lane-Emden Solution for n=1.5

The Lane-Emden equation for polytropic index n:

$$\frac{1}{\xi^2}\frac{d}{d\xi}\left(\xi^2 \frac{d\theta}{d\xi}\right) + \theta^n = 0$$

with boundary conditions θ(0) = 1, θ'(0) = 0.

For n = 1.5 (γ = 5/3):
- Surface radius: ξ₁ = 3.65375
- Central to mean density ratio: ρc/ρ̄ = 5.99071

### 3.2 Density Profile and Gradient

The density profile is:
$$\rho(r) = \rho_c \theta^n(\xi) = \rho_c \theta^{3/2}(\xi)$$

where ξ = r/α and α = R/ξ₁.

Near the center (ξ → 0), the Taylor expansion gives:
$$\theta(\xi) \approx 1 - \frac{\xi^2}{6} + \mathcal{O}(\xi^4)$$

The density gradient:
$$\frac{d\rho}{dr} = \frac{\rho_c}{\alpha} \cdot \frac{3}{2} \theta^{1/2} \frac{d\theta}{d\xi}$$

### 3.3 Characteristic Density Scale Length

At intermediate radii (where gradients are largest), the density scale length is:

$$L_\rho = \frac{\rho}{|d\rho/dr|} \sim \frac{\alpha}{\left|\frac{d\ln\theta}{d\xi}\right|} \sim \frac{R}{\xi_1} \cdot \mathcal{O}(1)$$

For the Lane-Emden n=1.5, the minimum density scale length occurs around ξ ≈ 2-3:

$$L_{\rho,\min} \sim \frac{R}{\xi_1} \approx 0.27 R$$

---

## 4. Error Scaling Derivation

### 4.1 Smoothing Length vs Density Scale Length

The smoothing length at the center is:

$$h_c = \eta \left(\frac{N_{nb} \cdot m}{\rho_c}\right)^{1/3}$$

where η ≈ 1.2-2.0 is the kernel-dependent factor.

For equal-mass particles with M = 1, m = 1/N:

$$h_c = \eta \left(\frac{N_{nb}}{N \cdot \rho_c}\right)^{1/3}$$

### 4.2 Resolution Criterion

For accurate density estimation, we require:
$$h \ll L_\rho$$

or equivalently:
$$\frac{h}{L_\rho} \ll 1$$

### 4.3 Error vs Particle Number

The fractional density error scales as:

$$\epsilon = \frac{|\rho_{SPH} - \rho_{exact}|}{\rho_{exact}} \sim C \left(\frac{h}{L_\rho}\right)^2$$

where C is an order-unity constant depending on the kernel and density profile.

Substituting h ~ N^(-1/3):

$$\boxed{\epsilon \sim C \cdot N^{-2/3}}$$

This is the **fundamental SPH convergence rate in 3D**.

---

## 5. Particle Count Formula

### 5.1 General Formula

From the error scaling, we can solve for the required particle count:

$$N(\epsilon) = N_0 \left(\frac{\epsilon_0}{\epsilon}\right)^{3/2}$$

where N₀ is a reference particle count with known error ε₀.

### 5.2 Calibration from Simulation

From our simulation with N₀ = 10,648 particles:
- Measured central density: ρc,sim = 1.66
- Expected central density: ρc,exact = 1.43
- Error: ε₀ = |1.66 - 1.43| / 1.43 ≈ **16%**

### 5.3 Required Particle Counts

| Target Error ε | N Required | Calculation |
|---------------|------------|-------------|
| 10% | ~21,000 | 10,648 × (0.16/0.10)^1.5 |
| 5% | ~61,000 | 10,648 × (0.16/0.05)^1.5 |
| 2% | ~240,000 | 10,648 × (0.16/0.02)^1.5 |
| 1% | ~680,000 | 10,648 × (0.16/0.01)^1.5 |

---

## 6. Alternative Derivation: Resolution Elements

### 6.1 Number of Resolution Elements

Another way to estimate particle requirements is to count the number of "resolution elements" needed across the sphere.

If we need Nres elements across the characteristic scale Lρ:
$$N_{elements} = \left(\frac{R}{L_\rho / N_{res}}\right)^3 = \left(\frac{R \cdot N_{res}}{L_\rho}\right)^3$$

With Lρ ~ R/ξ₁ and ξ₁ ≈ 3.65:
$$N_{elements} = (N_{res} \cdot \xi_1)^3 \approx 49 \cdot N_{res}^3$$

### 6.2 Resolution Requirements

| Nres (per Lρ) | Total Particles | Expected Error |
|---------------|-----------------|----------------|
| 2 | ~390 | ~50% (very poor) |
| 5 | ~6,100 | ~20% |
| 10 | ~49,000 | ~5% |
| 20 | ~390,000 | ~1% |

---

## 7. Practical Considerations

### 7.1 Kernel Choice Effects

Different kernels have different error constants C:
- **Cubic Spline:** C ≈ 1.0 (baseline)
- **Wendland C4:** C ≈ 0.8 (slightly better)
- **Quintic Spline:** C ≈ 0.7 (better but more expensive)

### 7.2 Neighbor Number Effects

Higher neighbor numbers reduce statistical noise but increase computational cost:
- $N_{nb}$ = 32: Higher noise, less accurate gradients
- $N_{nb}$ = 50: Good balance (typical choice)
- $N_{nb}$ = 100: Lower noise, more expensive

The error prefactor scales weakly with Nnb:
$$C \propto N_{nb}^{-1/3}$$

### 7.3 Initial Condition Quality

The **isentropic** requirement (K = const) is crucial:
- Without constant K, the equilibrium configuration differs from Lane-Emden
- This was the primary source of the ~16% error in our test
- With fixed K (using analytic density for internal energy), lower particle counts may achieve the same accuracy

---

## 8. Summary and Recommendations

### 8.1 Key Results

1. **Fundamental scaling:** Error ε ~ N^(-2/3) in 3D SPH

2. **Particle requirements for Lane-Emden n=1.5:**
   | Accuracy | Particles |
   |----------|-----------|
   | 10% | ~20,000 |
   | 5% | ~60,000 |
   | 2% | ~240,000 |
   | 1% | ~700,000 |

3. **Critical factor:** Initial K must be constant (isentropic) using analytic density

### 8.2 Recommendations

1. **For quick tests:** Use N ≈ 20,000-50,000 particles (5-10% accuracy)

2. **For production runs:** Use N ≈ 100,000-200,000 particles (2-3% accuracy)

3. **For convergence studies:** Run with N = 10k, 50k, 200k, 500k and verify N^(-2/3) scaling

4. **Always verify:** Check that K(r) = P/ρ^γ is constant to within a few percent

---

## 9. Mathematical Appendix

### A. Lane-Emden Parameters for n=1.5

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Polytropic index | n | 1.5 |
| Adiabatic index | γ | 5/3 |
| Dimensionless surface | ξ₁ | 3.6537540101 |
| dθ/dξ at surface | θ'₁ | -0.20330 |
| Central/mean density | ρc/ρ̄ | 5.99071 |

### B. Code Units (M=1, R=1)

| Quantity | Formula | Value |
|----------|---------|-------|
| Scaling factor | α = R/ξ₁ | 0.2737 |
| Central density | ρc = 5.99 × ρ̄ | 1.4302 |
| Mean density | ρ̄ = 3M/(4πR³) | 0.2387 |
| Polytropic K | 4πGα²ρc^(1/3)/2.5 | 0.4242 |

### C. Error Derivation Details

The SPH density estimator has the following properties:
1. **Zeroth moment:** ∫W dV = 1 (normalization)
2. **First moment:** ∫r W dV = 0 (symmetry)
3. **Second moment:** ∫r² W dV = σ²h² (kernel width)

Taylor expanding the density field:
$$\rho(\mathbf{r}') = \rho(\mathbf{r}) + (\mathbf{r}' - \mathbf{r}) \cdot \nabla\rho + \frac{1}{2}(\mathbf{r}' - \mathbf{r})^T \mathbf{H} (\mathbf{r}' - \mathbf{r}) + ...$$

The SPH estimate becomes:
$$\rho_{SPH} = \rho + \mathcal{O}(h^2 \nabla^2\rho)$$

The relative error:
$$\frac{\delta\rho}{\rho} = \frac{h^2 \nabla^2\rho}{\rho} \sim \frac{h^2}{L_\rho^2}$$

where $L_\rho^{-2} \sim \nabla^2\rho / \rho$ is the curvature scale.

---

*Document created: 2025-12-04*
*SPH Code Simulation Project*
