# Riemann Solver Implementation Gap Analysis

## Current Implementation: HLL Solver

The current implementation in `g_fluid_force.cpp::hll_solver()` uses a **non-iterative HLL (Harten-Lax-van Leer) solver**:

```cpp
void FluidForce::hll_solver()
{
    m_solver = [&](const real left[], const real right[], real & pstar, real & vstar) {
        const real u_l   = left[0];   // velocity
        const real rho_l = left[1];   // density
        const real p_l   = left[2];   // pressure
        const real c_l   = left[3];   // sound speed

        const real u_r   = right[0];
        const real rho_r = right[1];
        const real p_r   = right[2];
        const real c_r   = right[3];

        // Roe averaging
        const real roe_l = std::sqrt(rho_l);
        const real roe_r = std::sqrt(rho_r);
        const real roe_inv = 1.0 / (roe_l + roe_r);
        const real u_t = (roe_l * u_l + roe_r * u_r) * roe_inv;
        const real c_t = (roe_l * c_l + roe_r * c_r) * roe_inv;
        
        // Wave speeds
        const real s_l = std::min(u_l - c_l, u_t - c_t);
        const real s_r = std::max(u_r + c_r, u_t + c_t);

        // HLL state
        const real c1 = rho_l * (s_l - u_l);
        const real c2 = rho_r * (s_r - u_r);
        const real c3 = 1.0 / (c1 - c2);
        const real c4 = p_l - u_l * c1;
        const real c5 = p_r - u_r * c2;
        
        vstar = (c5 - c4) * c3;
        pstar = (c1 * c5 - c2 * c4) * c3;
    };
}
```

### Key Characteristics of Current HLL Implementation:
1. **Non-iterative**: Single-pass calculation
2. **Roe averaging**: Uses density-weighted averaging for intermediate states
3. **Wave speed estimation**: Uses minimum/maximum of characteristic speeds
4. **Direct formula**: Computes P* and v* directly without iteration

---

## Paper's Iterative Riemann Solver (van Leer 1997)

From **"Implementations and tests of Godunov-type particle hydrodynamics"** (Cha & Whitworth 2003), Section 3.2:

### Algorithm Overview:

The paper describes an **iterative Lagrangian Riemann solver** with these components:

#### 1. Lagrangian Shock Speed (W)
For each state (a or b):
```
W = √(W² + C²)  when P* ≥ P
W = C(1 - P*/P) / (1 - (P*/P)^((γ-1)/(2γ)))  otherwise (rarefaction)

where:
- C = √(γP/ρ) = Lagrangian sound speed
- γ = ratio of specific heats
```

#### 2. Tangential Slope (Z = dP*/dv*)
```
Z = ρW  when P* ≥ P (shock)

Z = ρC / √(1 + ((γ+1)/(2γ)) * (P*-P)/P)  when P* < P (rarefaction)
```

#### 3. Initial Pressure Estimate (Equation 15)
```
P*^(1) = (Ca*Pb + Cb*Pa - Ca*Cb*(va - vb)) / (Ca + Cb)

where Ca, Cb are sound speeds at states a and b
```

#### 4. Iteration (Equation 19)
```
P*^(n+1) = P*^(n) - (Zb^(n) * Za^(n) * (va*^(n) - vb*^(n))) / (Zb^(n) + Za^(n))
```

**Convergence criterion**: Iterate until P*^(n+1) differs by less than **1.5%** from P*^(n)

**Iteration trigger**: Only iterate if P*^(1) differs by more than **1%** from both Pa and Pb

#### 5. Velocity Updates at Each Iteration (Equations 20-21)
```
va*^(n) = va + (P*^(n) - Pa) / Wa
vb*^(n) = vb - (P*^(n) - Pb) / Wb
```

#### 6. Final Velocity (Equation 22)
```
v* = (Zb*vb* + Za*va*) / (Zb + Za)
```

#### 7. Pressure Protection
The paper mentions using **pressure protection** (van Leer 1997) to prevent divergence of the iteration, especially in strong rarefaction waves that can create fluid cavities (very low density regions).

---

## Implementation Gaps

### 1. **Iteration Loop** ❌
- **Paper**: Iterative refinement with convergence check (1.5% tolerance)
- **Current**: Single-pass calculation
- **Impact**: Lower accuracy for complex wave interactions

### 2. **Shock vs Rarefaction Detection** ❌
- **Paper**: Different formulas for W and Z depending on whether P* ≥ P (shock) or P* < P (rarefaction)
- **Current**: No distinction between shocks and rarefactions
- **Impact**: May not handle rarefaction waves correctly

### 3. **Lagrangian Shock Speed (W)** ❌
- **Paper**: Computes W based on pressure ratio and γ
- **Current**: Uses characteristic wave speeds (s_l, s_r) instead
- **Impact**: Different physical approach to wave propagation

### 4. **Tangential Slope (Z)** ❌
- **Paper**: Z = ρW or ρC/√(...) depending on shock/rarefaction
- **Current**: Not explicitly computed
- **Impact**: Missing physical quantity used in iteration

### 5. **Initial Estimate** ⚠️ Different approach
- **Paper**: P*^(1) = (Ca*Pb + Cb*Pa - Ca*Cb*(va - vb)) / (Ca + Cb)
- **Current**: Uses HLL averaging with different coefficients
- **Similarity**: Both use linear combination of left/right states

### 6. **Velocity Calculation** ⚠️ Different approach
- **Paper**: v* = (Zb*vb* + Za*va*) / (Zb + Za) with updated va*, vb*
- **Current**: vstar = (c5 - c4) * c3 (HLL formula)
- **Impact**: Different weighting scheme

### 7. **Pressure Protection** ❌
- **Paper**: Explicitly mentions pressure protection to prevent iteration divergence
- **Current**: No pressure protection mechanism
- **Impact**: May fail in extreme rarefaction scenarios (fluid cavities)

### 8. **Iteration Control** ❌
- **Paper**: 
  - Only iterate if |P*^(1) - Pa| > 1% and |P*^(1) - Pb| > 1%
  - Convergence when |P*^(n+1) - P*^(n)| / P*^(n) < 1.5%
- **Current**: No iteration control
- **Impact**: Missing adaptive accuracy mechanism

### 9. **Second-Order MUSCL** ✅ Implemented
- **Paper**: Describes MUSCL with van Leer limiter for gradient reconstruction
- **Current**: Implemented in `calculation()` with `m_is_2nd_order` flag
- **Status**: ✅ This part is correctly implemented

---

## Summary Table

| Feature | Paper (van Leer 1997) | Current HLL | Status |
|---------|----------------------|-------------|--------|
| Solver Type | Iterative Lagrangian | Non-iterative Eulerian | ❌ Different |
| Shock/Rarefaction | Separate formulas | Unified approach | ❌ Missing |
| Lagrangian W | Computed from P*, γ | Not used | ❌ Missing |
| Tangential slope Z | Computed explicitly | Not used | ❌ Missing |
| Iteration | Yes (1.5% tolerance) | No | ❌ Missing |
| Pressure protection | Yes | No | ❌ Missing |
| Initial estimate | Equation 15 | HLL averaging | ⚠️ Different |
| Velocity solution | Equation 22 | HLL formula | ⚠️ Different |
| MUSCL reconstruction | van Leer limiter | van Leer limiter | ✅ Implemented |
| Gradient calculation | Yes | Yes | ✅ Implemented |

---

## Recommendations

### Priority 1: Core Algorithm Replacement
Replace the current HLL solver with the iterative van Leer solver to match the paper's approach.

### Priority 2: Add Pressure Protection
Implement pressure protection mechanism to handle extreme rarefaction waves.

### Priority 3: Add Iteration Control
Implement convergence checking and iteration triggering logic.

### Priority 4: Performance Optimization
The iterative solver will be slower than HLL. Consider:
- Making it optional via configuration
- Providing both solvers as options
- Optimizing iteration convergence

### Priority 5: Testing
Test with problems mentioned in the paper:
- Shock tube tests
- Strong rarefaction waves
- Sjögreen test (Section 5.5)
- Fluid cavity formation scenarios

---

## Code Structure Suggestion

```cpp
void FluidForce::van_leer_iterative_solver()
{
    m_solver = [&](const real left[], const real right[], real & pstar, real & vstar) {
        // Extract states
        const real u_a = left[0], rho_a = left[1], p_a = left[2], c_a = left[3];
        const real u_b = right[0], rho_b = right[1], p_b = right[2], c_b = right[3];
        
        // Initial estimate (Eq 15)
        real p_star = (c_a * p_b + c_b * p_a - c_a * c_b * (u_a - u_b)) / (c_a + c_b);
        
        // Check if iteration is needed (1% threshold)
        const real threshold = 0.01;
        if (std::abs(p_star - p_a) / p_a < threshold && 
            std::abs(p_star - p_b) / p_b < threshold) {
            // Skip iteration, use initial estimate
            vstar = compute_final_velocity(p_star, ...);
            pstar = p_star;
            return;
        }
        
        // Iterate until convergence (1.5% tolerance)
        const real conv_tol = 0.015;
        const int max_iter = 100;
        
        for (int iter = 0; iter < max_iter; ++iter) {
            // Compute W for both sides (shock vs rarefaction)
            real w_a = compute_lagrangian_shock_speed(p_star, p_a, rho_a, c_a);
            real w_b = compute_lagrangian_shock_speed(p_star, p_b, rho_b, c_b);
            
            // Compute Z (tangential slope)
            real z_a = compute_tangential_slope(p_star, p_a, rho_a, c_a, w_a);
            real z_b = compute_tangential_slope(p_star, p_b, rho_b, c_b, w_b);
            
            // Update velocities (Eq 20-21)
            real u_a_star = u_a + (p_star - p_a) / w_a;
            real u_b_star = u_b - (p_star - p_b) / w_b;
            
            // Update pressure (Eq 19)
            real p_star_new = p_star - (z_b * z_a * (u_a_star - u_b_star)) / (z_b + z_a);
            
            // Apply pressure protection
            p_star_new = apply_pressure_protection(p_star_new, p_a, p_b);
            
            // Check convergence
            if (std::abs(p_star_new - p_star) / p_star < conv_tol) {
                p_star = p_star_new;
                break;
            }
            
            p_star = p_star_new;
        }
        
        // Compute final velocity (Eq 22)
        vstar = compute_final_velocity(p_star, u_a, u_b, ...);
        pstar = p_star;
    };
}
```

---

## References

- **Paper**: Cha, S.-H., & Whitworth, A. P. (2003). "Implementations and tests of Godunov-type particle hydrodynamics." *Monthly Notices of the Royal Astronomical Society*, 340(1), 73-90.
- **van Leer solver**: van Leer, B. (1997) - Referenced in the paper for the iterative Riemann solver algorithm
- **Current implementation**: `/Users/guo/Downloads/sphcode/src/gsph/g_fluid_force.cpp`
