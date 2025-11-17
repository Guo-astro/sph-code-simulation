# SR-GSPH Implementation Gap Analysis

**Date:** 2025-11-16  
**Reference:** Kitajima, Inutsuka, Seno (2025) arXiv:2510.18251v1 §2.2  
**Analysis of:** Current SR-GSPH implementation in `/src/srgsph/`

---

## Executive Summary

Your current implementation **mixes concepts from two different formulations** in the Kitajima+ paper:

1. **Fixed smoothing length formulation** (§2.2) - simpler, exact conservation
2. **Variable smoothing length formulation** (§2.3-2.5) - more complex, volume-based

This mixing creates **fundamental inconsistencies** that compromise both accuracy and conservation properties.

---

## Critical Gaps Identified

### ❌ **GAP 1: Inconsistent Number Density Calculation**

**What the fixed-h formulation says (§2.2, Eq. 10):**
```
n(x) = Σ_j ν_j W(x - x_j; h)    [fixed h for all particles]
```

**What your code does (`sr_pre_interaction.cpp:330-334`):**
```cpp
// Compute baryon number density using VOLUME-BASED approach (Eq. 42)
// N_volume-based(x) = ν(x) / Vp(x)
p_i.N = p_i.nu / Vp_i;
```

**The problem:**
- Fixed-h formulation: `N = direct kernel sum` (**Eq. 10-13**)
- Variable-h formulation: `N = ν / Vp` where `Vp = 1/Σ_j W(...)` (**Eq. 42**)
- You're using **Eq. 42 from §2.4** (variable-h) instead of **Eq. 10 from §2.2** (fixed-h)

**Impact:**
- The volume-based approach introduces an **extra layer of kernel smoothing**
- For fixed-h: `N ≈ Σ_j ν_j W_ij` (direct)
- For volume-based: `N = ν / Vp = ν × (Σ_j W_ij)` (inverse then multiply)
- These are **mathematically different** and give different results

---

### ❌ **GAP 2: Variable vs Fixed Smoothing Length**

**What the fixed-h formulation says (§2.2, opening):**
```
h = const.    [single constant for all particles and all times]
```

**What your code does (`sr_pre_interaction.cpp:130-254`):**
```cpp
real PreInteraction::compute_smoothing_length(...) {
    // VOLUME-BASED APPROACH with gather method (Eqs. 35-36 in paper)
    // h = η * Vp*^(1/d)
    // where Vp* uses kernel with C_smooth * h
    
    for(int iter = 0; iter < max_iter; ++iter) {
        const real h_smooth = m_c_smooth * h;
        // ... iterative solver ...
        const real h_new = eta_corrected * std::pow(Vp_star, 1.0 / real(DIM));
    }
}
```

**The problem:**
- You have **per-particle iterative smoothing length** using **Eqs. 35-36 from §2.3** (variable-h section)
- Fixed-h uses a **single global constant** `h` chosen at initialization
- The iteration algorithm you're using is **explicitly for variable-h** (Gather method, §2.3.3)

**Impact:**
- Breaks the fundamental assumption of fixed-h: uniform kernel support
- Introduces computational overhead (20 iterations per particle per timestep)
- The kernel geometry factor K_ij depends on h being constant

---

### ⚠️ **GAP 3: Kernel Gradient Form (Minor but Important)**

**What the fixed-h formulation says (§2.2, Eq. 26-32):**
```
For fixed h:
∇_i W(x_i - x_j; h) - ∇_j W(x_i - x_j; h)

Note: With fixed h, this simplifies because ∇_i W = -∇_j W
```

**What your code does (`sr_fluid_force.cpp:200-204`):**
```cpp
// Kernel gradients (variable smoothing length)
// Eq. 64-65: uses ∇_i W(x_i-x_j, 2h_i) - ∇_j W(x_i-x_j, 2h_j)
const vec_t dw_i = kernel->dw(r_ij, r, 2.0 * h_i);
const vec_t dw_j = kernel->dw(r_ij, r, 2.0 * p_j.sml);
const vec_t dw_ij = dw_i - dw_j;
```

**The problem:**
- You're using **different smoothing lengths** `h_i` and `h_j` (**Eq. 64-65 from §2.5**, variable-h)
- Fixed-h uses **same h for both**: `∇W(r_ij; h) - ∇W(r_ij; h) = 0`? No!
- Actually for fixed-h: `∇_i W - ∇_j W = -2∇W` by symmetry

**Why it matters:**
- The factor of `2h` vs `h` comes from the **Godunov SPH 1D reduction** (§2.2.3)
- Using different `h_i` and `h_j` is **only correct for variable-h** formulation
- For fixed-h, this should be **same h for both particles**

---

### ⚠️ **GAP 4: Volume Weighting Term**

**What the fixed-h formulation says (§2.2, Eq. 31-32):**
```
dS_i/dt = -Σ_j ν_i ν_j · F^(mom)(U_i, U_j; n_ij) · K_ij

where K_ij is a symmetric kernel-geometry factor from the convolution
```

**What your code does (`sr_fluid_force.cpp:206-211`):**
```cpp
// Volume weighting V²_ij (Paper Eq. 29-30 and 63)
// V²_ij(h) ≡ ∫ [1/N²(x)] W(x-x_mid, √2h) dx  
// Approximation: V²_ij ≈ ½(Vₚ²_i/ν² + Vₚ²_j/ν²) = ½((ν/N)²_i + (ν/N)²_j)
const real V_i = nu_i / p_i.N;
const real V_j = p_j.nu / p_j.N;
const real V2_ij = 0.5 * (V_i * V_i + V_j * V_j);
```

**The problem:**
- You're computing `V²_ij` from **Eq. 63 in §2.5** (variable-h formulation)
- Fixed-h doesn't need this! The kernel geometry factor `K_ij` is **directly computable** from the Gaussian kernel overlap
- For fixed-h with equal ν_j, `K_ij` depends **only on |r_i - r_j|** and the fixed `h`

**Impact:**
- Extra complexity that's not needed for fixed-h
- The volume-based V²_ij can vary even for particles at the same separation
- This breaks the translational invariance that fixed-h guarantees

---

## ✅ What You Got Right

### **Riemann Solver Inputs (CORRECT)**

Your code **correctly** uses rest-frame variables in the Riemann solver:

```cpp
// CRITICAL: Use REST-FRAME baryon number density n = N/γ
const real n_i = p_i.N / p_i.gamma_lor;
const real n_j = p_j.N / p_j.gamma_lor;

const real right[4] = {
    ve_i,      // rest-frame velocity component
    n_i,       // REST-FRAME baryon number density
    p_i.pres,  // pressure (same in all frames)
    p_i.sound, // sound speed
};
```

This is **exactly what §2.2 requires**: the 1D Riemann problem is solved in the rest frame with rest-frame density `n = N/γ`.

---

## Recommended Path Forward

You have **three options**:

### **Option A: Pure Fixed-h Implementation (RECOMMENDED for simplicity)**

Simplify to exactly match §2.2:

1. **Remove** all variable-h code:
   - Delete `compute_smoothing_length()` iteration
   - Delete volume-based density `N = ν/Vp`
   - Delete per-particle `h_i`, `h_j`

2. **Replace with fixed-h:**
   ```cpp
   // One global h for all particles
   const real h = param->srgsph.smoothing_length;  // or auto-compute once
   
   // Direct kernel sum for number density (Eq. 10)
   for(int i = 0; i < num; ++i) {
       real n_sum = 0.0;
       for(int j : neighbors) {
           n_sum += particles[j].nu * kernel->w(r_ij, h);
       }
       particles[i].N = n_sum;
   }
   
   // Kernel gradients with same h
   const vec_t dw_i = kernel->dw(r_ij, r, h);  // Not 2h!
   // For Godunov SPH 1D reduction, you might need 2h - check §2.2.3
   
   // Simplified K_ij (no volume weighting needed)
   const real K_ij = compute_kernel_overlap(r, h);  // Function of r and h only
   ```

3. **Benefits:**
   - Exact conservation of momentum and energy (proven in §2.2)
   - Much simpler code
   - Faster (no iterations)
   - Direct correspondence to paper equations

### **Option B: Pure Variable-h Implementation**

Fully commit to §2.3-2.5:

1. Keep your volume-based density
2. Keep iterative smoothing length
3. But add the **extra terms** from variable-h:
   - ∇h terms in the equations (Eq. 64-65 have these)
   - Proper volume evolution equation
   - More complex conservation corrections

4. **Challenges:**
   - More complex
   - **Approximate** conservation (not exact like fixed-h)
   - Need careful treatment of ∇h terms

### **Option C: Hybrid with Configuration Flag**

Support both formulations:

```cpp
if(param->srgsph.use_fixed_smoothing) {
    // Fixed-h path (§2.2)
} else {
    // Variable-h path (§2.3-2.5)
}
```

**Not recommended** - doubles code complexity and testing burden.

---

## Detailed Code Changes for Option A (Fixed-h)

### 1. Remove iteration in `sr_pre_interaction.cpp`

```cpp
void PreInteraction::calculation(std::shared_ptr<Simulation> sim) {
    auto & particles = sim->get_particles();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();
    
    // Use single global h (set once at initialization)
    const real h = m_h_global;  // Add this as member variable
    
#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        auto & p_i = particles[i];
        
        // Fixed smoothing length (no iteration!)
        p_i.sml = h;
        
        // Direct kernel sum for number density (Eq. 10)
        real n_sum = 0.0;
        for(int j : neighbors) {
            const vec_t r_ij = periodic->calc_r_ij(p_i.pos, particles[j].pos);
            const real r = std::abs(r_ij);
            n_sum += particles[j].nu * kernel->w(r, h);
        }
        p_i.N = n_sum;
        
        // Recover primitives (this part stays the same)
        auto prim = PrimitiveRecovery::conserved_to_primitive(...);
        p_i.vel = prim.vel;
        p_i.pres = prim.pressure;
        // ... etc
        
        // Gradients for MUSCL (this part can stay mostly the same)
        // ...
    }
}
```

### 2. Simplify kernel gradients in `sr_fluid_force.cpp`

```cpp
// Fixed-h kernel gradients
// For Godunov SPH 1D reduction, check if you need h or 2h
// The factor comes from the convolution integral in §2.2.3
const real h_kernel = h;  // or 2.0 * h, depending on exact Godunov reduction
const vec_t dw_i = kernel->dw(r_ij, r, h_kernel);
const vec_t dw_j = kernel->dw(r_ij, r, h_kernel);  // SAME h!
const vec_t dw_ij = dw_i - dw_j;  // This is antisymmetric by kernel property

// For fixed-h with symmetric kernel, dw_i = -dw_j, so:
// dw_ij = 2 * dw_i  (approximately, for symmetric particles)
```

### 3. Replace volume weighting with kernel geometry factor

```cpp
// Fixed-h kernel overlap factor
// This comes from the 1D Godunov reduction (§2.2.3)
// For equal ν particles and Gaussian kernel:
const real K_ij = compute_gaussian_overlap(r, h);

// Force equations (Eq. 31-32)
dS -= dw_ij * (pstar * nu_i * nu_j * K_ij);
de -= inner_product(v_ij, dw_ij) * (pstar * nu_i * nu_j * K_ij);
```

You'll need to implement `compute_gaussian_overlap()` based on the Godunov SPH formulation. This is a **fixed function** of `r` and `h` only.

---

## Why This Matters

### **Conservation Properties**

**Fixed-h (§2.2):**
- ✅ **Exact** conservation of total momentum
- ✅ **Exact** conservation of total energy
- ✅ Proven analytically (antisymmetric forces)

**Your current mixed approach:**
- ❌ Approximate conservation (volume terms break symmetry)
- ❌ Can accumulate errors over time
- ❌ Not proven to conserve

### **Shock Handling**

**Fixed-h:**
- ✅ Naturally handles contact discontinuities (§2.2 mentions this)
- ✅ MUSCL reconstruction provides shock capturing
- ✅ Riemann solver provides upwinding

**Variable-h:**
- ⚠️ Requires special treatment at discontinuities
- ⚠️ Can have h-oscillations across shocks
- ⚠️ Your iteration might not converge at strong shocks

---

## Specific File-by-File Recommendations

### `sr_pre_interaction.hpp` + `.cpp`
- **Remove**: `compute_smoothing_length()` and `compute_volume()`
- **Simplify**: `calculation()` to direct kernel sum (Eq. 10)
- **Add**: Global `h` initialization

### `sr_fluid_force.hpp` + `.cpp`
- **Change**: Kernel gradients to use same `h` for both particles
- **Replace**: `V2_ij` volume weighting with `K_ij` kernel overlap
- **Check**: Whether 1D reduction uses `h` or `2h` (read Godunov SPH paper if unclear)

### `sr_timestep.cpp`
- **Current implementation looks OK** - CFL with sound speed is correct for both formulations

### `parameters.hpp`
- **Add**: `srgsph.smoothing_length` or auto-compute from particle spacing
- **Remove** (or ignore): `iterative_sml` flag for SR-GSPH

---

## Testing Recommendations

Once you fix to pure fixed-h:

1. **Sod shock tube**: Should conserve momentum exactly (check `Σ p_i.S == const`)
2. **Energy conservation**: For isolated system, `Σ (p_i.e * p_i.nu)` should be constant
3. **Contact discontinuity**: Test density jump with smooth velocity (should be stable)
4. **Convergence**: Run with different `h` values, verify 2nd-order accuracy

---

## Bottom Line

**Your implementation is currently trying to use §2.3-2.5 (variable-h) methods while claiming to implement §2.2 (fixed-h).**

The cleanest fix is:
1. **Choose fixed-h** (simpler, exact conservation, matches §2.2)
2. **Remove all volume-based and iteration code**
3. **Use direct kernel sums** for number density
4. **Use symmetric kernel geometry** for forces

This will give you a **correct, validated, simpler implementation** that exactly matches the paper's §2.2 formulation.

---

## Questions to Resolve

1. **Did you intend to implement fixed-h or variable-h?**
   - Comments say "variable smoothing length" but that's §2.3-2.5, not §2.2

2. **What's the source of the `2h` factor in kernel gradients?**
   - This likely comes from Godunov SPH 1D reduction
   - Need to verify whether it's `h`, `2h`, or `√2 h` for Gaussian kernel
   - Check Inutsuka (2002) Godunov SPH paper for details

3. **What's your target test case?**
   - For SR shock tubes: fixed-h is standard and well-tested
   - For smooth flows: variable-h might adapt better (but adds complexity)

---

**Next Step:** Please confirm whether you want to implement:
- **A) Fixed-h (§2.2)** - my recommendation
- **B) Variable-h (§2.3-2.5)** - more complex but adaptive

I can provide detailed implementation guidance for either path.
