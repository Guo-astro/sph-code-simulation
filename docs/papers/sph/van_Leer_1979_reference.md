# van Leer (1979) - Riemann Problem Solution Reference

## Citation

**van Leer, B. (1979). Towards the ultimate conservative difference scheme. V. A second-order sequel to Godunov's method.**  
*Journal of Computational Physics*, 32(1), 101-136.  
DOI: [10.1016/0021-9991(79)90145-1](https://doi.org/10.1016/0021-9991(79)90145-1)

**Highly Cited**: Over 9,800 citations (as of 2025)

## Overview

This seminal 1979 paper by Bram van Leer presents a detailed explanation for the solution of the one-dimensional Riemann problem within the context of developing second-order accurate conservative finite-difference schemes for hyperbolic conservation laws.

## Key Contributions

### 1. Exact Riemann Problem Solution

The paper provides a comprehensive treatment of solving the 1D Riemann problem for the Euler equations of gas dynamics:

**Governing Equations:**
```
∂ρ/∂t + ∂(ρu)/∂x = 0           (mass)
∂(ρu)/∂t + ∂(ρu² + P)/∂x = 0   (momentum)
∂E/∂t + ∂[(E + P)u]/∂x = 0     (energy)
```

where:
- ρ = density
- u = velocity
- P = pressure
- E = total energy per unit volume = ρ(e + u²/2)
- e = specific internal energy = P/[ρ(γ-1)]

### 2. Riemann Problem Structure

The Riemann problem consists of initial conditions:
```
U(x,0) = U_L  for x < 0
U(x,0) = U_R  for x > 0
```

The solution consists of three waves separating four constant states:

```
U_L | Rarefaction/Shock | U*_L | Contact | U*_R | Rarefaction/Shock | U_R
```

### 3. Wave Types and Solution Method

#### Shock Waves

For a shock wave, the Rankine-Hugoniot jump conditions apply:
```
ρ_R(u_R - S) = ρ_L(u_L - S)
ρ_R(u_R - S)u_R + P_R = ρ_L(u_L - S)u_L + P_L
ρ_R(u_R - S)[(e_R + u_R²/2) + P_R/ρ_R] = ρ_L(u_L - S)[(e_L + u_L²/2) + P_L/ρ_L]
```

where S is the shock speed.

#### Rarefaction Waves

For a rarefaction fan, the solution is self-similar (depends on ξ = x/t):
```
ρ = ρ_0 * [1 + (γ-1)/(2c_0) * (u_0 - ξ)]^(2/(γ-1))
u = 2/(γ+1) * [c_0 + (γ-1)/2 * u_0 + ξ]
P = P_0 * [1 + (γ-1)/(2c_0) * (u_0 - ξ)]^(2γ/(γ-1))
```

where c = √(γP/ρ) is the sound speed.

#### Contact Discontinuity

Across the contact discontinuity:
- Pressure is continuous: P*_L = P*_R = P*
- Velocity is continuous: u*_L = u*_R = u*
- Density is discontinuous: ρ*_L ≠ ρ*_R

### 4. Iterative Solution Procedure

The Riemann problem is solved by finding P* and u* such that:

1. **Pressure function**: f(P*) = 0, where f relates the star states to initial states
2. **Newton-Raphson iteration**: P*^(n+1) = P*^n - f(P*^n)/f'(P*^n)
3. **Convergence**: Typically converges in 2-5 iterations

### 5. Special Cases

#### Vacuum Formation

Vacuum forms when:
```
2/(γ-1) * (c_L + c_R) ≤ (u_R - u_L)
```

In this case, the solution consists of two rarefaction waves with a vacuum region between them.

#### Strong Shock Limit

For strong shocks (P*/P → ∞):
```
ρ*/ρ → (γ+1)/(γ-1)
```

This is the maximum density ratio across a shock.

## Application to SPH

In SPH (Smoothed Particle Hydrodynamics), the Riemann problem solution is used for:

1. **Approximate Riemann Solvers**: Computing fluxes at particle interfaces
   - HLLC (Harten-Lax-van Leer-Contact) solver
   - Roe solver
   - Exact Riemann solver

2. **Test Cases**: Validating SPH implementations
   - Shock tube (Sod problem)
   - Vacuum test
   - Strong shock tests

3. **Boundary Conditions**: Handling inflow/outflow and wall boundaries

## van Leer's Second-Order Method

The paper introduces the MUSCL (Monotonic Upstream-centered Scheme for Conservation Laws) approach:

1. **Piecewise linear reconstruction** of states within cells
2. **Slope limiting** to prevent oscillations near discontinuities
3. **Riemann solver** at cell interfaces
4. **Conservative update** of cell averages

This approach became the foundation for modern high-resolution shock-capturing schemes.

## Key Formulas from van Leer (1979)

### Pressure in Star Region

For two rarefaction waves:
```
P* = [(c_L + c_R - (γ-1)/2 * (u_R - u_L)) / 
      (c_L/P_L^((γ-1)/(2γ)) + c_R/P_R^((γ-1)/(2γ)))]^(2γ/(γ-1))
```

For shock-shock case, use iterative solver.

### Velocity in Star Region

```
u* = 1/2 * (u_L + u_R) + 1/2 * (f_R(P*) - f_L(P*))
```

where f_L and f_R are the shock/rarefaction functions.

### Density in Star Region

For rarefaction:
```
ρ* = ρ * (P*/P)^(1/γ)
```

For shock:
```
ρ* = ρ * [(γ+1)*P* + (γ-1)*P] / [(γ-1)*P* + (γ+1)*P]
```

## Comparison with Other Methods

| Method | Accuracy | Computational Cost | Robustness |
|--------|----------|-------------------|------------|
| Exact Riemann | Exact | High (iterative) | Excellent |
| HLLC | Very Good | Low | Very Good |
| HLL | Good | Very Low | Excellent |
| Roe | Very Good | Medium | Good (entropy fix) |

## Implementation Notes

### Numerical Considerations

1. **Pressure positivity**: Ensure P > 0 in all states
2. **Vacuum detection**: Check vacuum formation condition
3. **Iteration tolerance**: Use relative error ~ 10^-6 for P*
4. **Entropy fix**: Apply correction for transonic rarefactions

### Performance Tips

1. **Initial guess**: Use P* = (P_L + P_R)/2 for Newton iteration
2. **Caching**: Store c_L, c_R to avoid repeated sqrt operations
3. **Simplified solver**: For smooth flows, use linearized Riemann solver

## Related Papers in the van Leer Series

The "Towards the Ultimate Conservative Difference Scheme" series consists of:

- **Part I** (1973): Introduction and basic concepts
- **Part II** (1974): Monotonicity and conservation
- **Part III** (1977): Upstream-centered finite-difference schemes
- **Part IV** (1977): New approach to numerical convection
- **Part V** (1979): **[THIS PAPER]** Second-order sequel to Godunov's method

## Code Implementation Example (Pseudocode)

```cpp
// Riemann solver for Euler equations (simplified)
void solveRiemann(State left, State right, State& star_left, State& star_right) {
    // 1. Compute sound speeds
    double c_L = sqrt(gamma * left.P / left.rho);
    double c_R = sqrt(gamma * right.P / right.rho);
    
    // 2. Check for vacuum
    if (2.0/(gamma-1) * (c_L + c_R) <= (right.u - left.u)) {
        solveVacuumCase(left, right, star_left, star_right);
        return;
    }
    
    // 3. Iterate for P*
    double P_star = 0.5 * (left.P + right.P);  // Initial guess
    for (int iter = 0; iter < maxIter; ++iter) {
        double f = pressureFunction(P_star, left, right);
        double df = pressureFunctionDerivative(P_star, left, right);
        P_star = P_star - f / df;
        if (abs(f) < tolerance) break;
    }
    
    // 4. Compute u*
    double u_star = computeStarVelocity(P_star, left, right);
    
    // 5. Compute star states
    star_left = computeStarState(P_star, u_star, left, "left");
    star_right = computeStarState(P_star, u_star, right, "right");
}
```

## Additional Resources

### Textbooks

1. **Toro, E. F. (2009). Riemann Solvers and Numerical Methods for Fluid Dynamics.**
   - Chapter 4: Complete treatment of exact Riemann solver
   - Includes implementation details and test cases

2. **LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems.**
   - Chapter 14: Riemann problems for Euler equations
   - Excellent pedagogical treatment

### Online Resources

1. **University of Michigan Digital Library**
   - URL: https://deepblue.lib.umich.edu/
   - Collection of van Leer's papers

2. **Clawpack Riemann Solvers**
   - URL: https://www.clawpack.org/riemann.html
   - Open-source implementations

## Verification Test Cases

### 1. Sod's Shock Tube (1978)

Initial conditions:
```
Left:  ρ=1.0, u=0.0, P=1.0
Right: ρ=0.125, u=0.0, P=0.1
γ = 1.4, t = 0.2
```

Expected structure: Left rarefaction, contact, right shock

### 2. Vacuum Test

Initial conditions:
```
Left:  ρ=1.0, u=-2.0, P=0.4
Right: ρ=1.0, u=+2.0, P=0.4
γ = 1.4, t = 0.15
```

Expected structure: Two rarefactions with vacuum in between

### 3. Strong Shock

Initial conditions:
```
Left:  ρ=1.0, u=0.0, P=1000.0
Right: ρ=1.0, u=0.0, P=0.01
γ = 1.4, t = 0.012
```

Expected structure: Left rarefaction, contact, strong right shock

## Historical Context

van Leer's 1979 paper built upon Godunov's seminal 1959 work on conservative finite-difference schemes. The key innovation was achieving second-order accuracy while maintaining monotonicity and conservation, which was previously thought to be impossible (Godunov's barrier theorem).

The MUSCL approach introduced in this paper became the foundation for:
- Modern CFD codes
- SPH approximate Riemann solvers
- High-resolution shock-capturing schemes
- Adaptive mesh refinement methods

## Citation in Research

When citing this work for SPH or Riemann solver implementations, use:

```bibtex
@article{vanLeer1979,
  title={Towards the ultimate conservative difference scheme. V. A second-order sequel to Godunov's method},
  author={van Leer, Bram},
  journal={Journal of computational Physics},
  volume={32},
  number={1},
  pages={101--136},
  year={1979},
  publisher={Elsevier},
  doi={10.1016/0021-9991(79)90145-1}
}
```

## Notes for This Repository

The vacuum test case in `sample/vacuum/` implements the analytical solution based on the Riemann problem formulation described in van Leer (1979). The rarefaction wave formulas and vacuum detection criteria follow directly from this paper's treatment.

Key implementation files:
- `sample/vacuum/scripts/vacuum_analytical.py`: Analytical Riemann solver
- `src/sample/vacuum.cpp`: Initial conditions setup
- `sample/vacuum/Makefile.vacuum`: Test automation

---

**Document prepared by**: SPH Code Development Team  
**Last updated**: November 12, 2025  
**Status**: Reference material for vacuum test implementation

## Access Information

**Official Paper Access**:
- **Publisher**: Elsevier (Science Direct)
- **DOI**: https://doi.org/10.1016/0021-9991(79)90145-1
- **University Library**: Available through institutional access
- **ResearchGate**: Author's profile may have preprint

**Alternative Access**:
- Many universities have institutional subscriptions to Journal of Computational Physics
- Author's personal website may have preprints
- Google Scholar provides citation tracking and may link to open-access versions

**Note**: This is a copyrighted work. Please access through legitimate channels such as university libraries or purchase from the publisher.
