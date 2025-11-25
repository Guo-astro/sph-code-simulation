# Riemann Solver Comparison: Fortran MM94 vs C++ Implementation

## Critical Differences Found

### 1. **SHOCK VELOCITY FORMULA - MAJOR DIFFERENCE**

#### Fortran MM94 (lines 673-675):
```fortran
A      = J**2+(RHOA*WA)**2
B      = -VELA*RHOA**2*WA**2
VSHOCK = (-B+SIGN*J**2*DSQRT(1.D0+RHOA**2/J**2))/A
```

Simplified:
```
Vs = (v*D² ± j*√(j² + D²)) / (j² + D²)
```
where `D = ρ*W` and `W = 1/√(1-v²)`

#### C++ Pons et al. (sr_exact_riemann.cpp lines 67-76):
```cpp
const real term = j * j + Da * Da * (1.0 - va * va);
const real sqrt_term = std::sqrt(term);
const real denom_vs = Da * Da + j * j;

Vs = (Da * Da * va ± j * sqrt_term) / denom_vs;
```

Simplified:
```
Vs = (v*D² ± j*√(j² + D²*(1-v²))) / (j² + D²)
```

#### KEY DIFFERENCE:
- **MM94**: `√(j² + D²)`
- **Pons**: `√(j² + D²*(1-v²))`

When `W² = 1/(1-v²)`, this becomes:
- **MM94**: `√(j² + ρ²/(1-v²))`
- **Pons**: `√(j² + ρ²)`

**For 1D normal shocks (vt=0), these give VERY DIFFERENT shock speeds!**

---

### 2. **LORENTZ FACTOR COMPUTATION**

#### Fortran MM94 (lines 104-105):
```fortran
WL  = 1.D0/DSQRT(1.D0-VELL**2)
```
Uses **NORMAL velocity only** (VELL is 1D velocity from line 76)

#### C++ (sr_exact_riemann.cpp lines 61-62):
```cpp
const real v2a_total = va * va + vta * vta;  // Total velocity squared
const real Wa = 1.0 / std::sqrt(1.0 - v2a_total);
```
Uses **TOTAL velocity** (normal + tangential)

**ISSUE**: For pure 1D shocks (vt=0), these are equivalent. But the formulas themselves differ in philosophy:
- MM94: W from normal velocity → D = ρ*W_normal
- Pons: W from total velocity → D = ρ*W_total

---

### 3. **POST-SHOCK VELOCITY FORMULA**

#### Fortran MM94 (lines 683-686):
```fortran
A = WSHOCK*(P-PA)/J+HA*WA*VELA
B = HA*WA+(P-PA)*(WSHOCK*VELA/J+1.D0/RHOA/WA)
VEL = A/B
```

#### C++ (sr_exact_riemann.cpp lines 87-92):
```cpp
const real num = ha * Wa * va + Ws * (pb - pa) / j_signed;
const real den = ha * Wa + (pb - pa) * (Ws * va / j_signed + 1.0 / Da);
vxb = num / den;
```

**These match!** ✓

---

### 4. **ENTHALPY CALCULATION**

#### Fortran MM94 (lines 628-645):
```fortran
A  = 1.D0+(GAMMA-1.D0)*(PA-P)/GAMMA/P
B  = 1.D0-A
C  = HA*(PA-P)/RHOA-HA**2
H = (-B+DSQRT(B**2-4.D0*A*C))/2.D0/A
```

Let me expand:
- `A = 1 + (γ-1)(P_a-P_b)/(γ*P_b)`
- `B = 1 - A = -(γ-1)(P_a-P_b)/(γ*P_b)`
- `C = h_a(P_a-P_b)/ρ_a - h_a²`

Quadratic: `A*h² + B*h + C = 0`

#### C++ (sr_exact_riemann.cpp lines 29-39):
```cpp
const real A = (gamma - 1.0) * (pa - pb) / (gamma * pb);
const real B = ha * (pa - pb) / rhoa;

const real qa = 1.0 + A;
const real qb = -A;
const real qc = B - ha * ha;

hb = (-qb + std::sqrt(delta)) / (2.0 * qa);
```

Quadratic coefficients:
- `qa = 1 + A` ✓
- `qb = -A`
- `qc = B - ha²` ✓

**Wait, qb is different!**
- Fortran: `B = 1 - A`
- C++: `qb = -A`

But `B = 1 - A = 1 - A`, while `qb = -A`.
So Fortran B ≠ C++ qb!

Let me recalculate. Fortran quadratic is `A*h² + B*h + C = 0` where:
- `A = 1 + (γ-1)(P_a-P_b)/(γ*P_b)`
- `B = 1 - [1 + (γ-1)(P_a-P_b)/(γ*P_b)] = -(γ-1)(P_a-P_b)/(γ*P_b)`
- `C = h_a(P_a-P_b)/ρ_a - h_a²`

C++ quadratic is `qa*h² + qb*h + qc = 0` where:
- `qa = 1 + (γ-1)(P_a-P_b)/(γ*P_b)` ✓
- `qb = -(γ-1)(P_a-P_b)/(γ*P_b)` ✓  (This equals Fortran B!)
- `qc = h_a(P_a-P_b)/ρ_a - h_a²` ✓

**These match!** ✓ (The variable naming is different but the math is the same)

---

## SUMMARY OF POTENTIAL BUGS

### ⚠️ BUG #1: Inconsistent Shock Speed Formula

**Location**: `sr_exact_riemann.cpp:67`

**Issue**: Using Pons et al. (2000) formula with `√(j² + D²*(1-v²))` factor, but Fortran MM94 uses `√(j² + D²)`.

**Impact**:
- For 1D Sod shock (v ≈ 0.3-0.5), the difference is significant!
- MM94: `√(j² + ρ²/(1-v²))` grows as v increases
- Pons: `√(j² + ρ²)` is constant
- This affects shock propagation speed and can cause errors in shock position

**For the 1D Sod test (no tangential velocity):**
- Should use MM94 formula: `√(j² + D²)` where `D = ρ*W` and `W = 1/√(1-v²)`
- Currently using: `√(j² + D²*(1-v²))` which gives `√(j² + ρ²)`

**Recommendation**:
- For `vt = 0` (pure normal shock), use MM94 formula
- For `vt ≠ 0` (shock with tangential component), use Pons et al. formula

---

### ⚠️ BUG #2: Possible Confusion Between Normal and Total Velocity

**Location**: `sr_exact_riemann.cpp:61-62`

**Issue**: Computing Lorentz factor from total velocity `(vx² + vt²)`, but then using in formula that expects W from normal velocity.

**Fortran approach** (1D only):
```
W = 1/√(1-v_normal²)
D = ρ*W
shock_speed_term = √(j² + D²)
```

**C++ approach** (with tangent):
```cpp
W = 1/√(1-(vx² + vt²))  // W from TOTAL velocity
D = ρ*W
shock_speed_term = √(j² + D²*(1-vx²))  // But uses vx only in (1-vx²) factor!
```

**This is INCONSISTENT!** The Pons formula likely expects:
- Either: `W` from normal velocity AND `√(j² + D²*(1-v²))`
- Or: `W` from total velocity AND different formula

**Need to verify**: Which Lorentz factor should be used in the Pons formulation?

---

## RECOMMENDATIONS

1. **For 1D normal shocks (vt=0)**:
   - Use MM94 formulation exactly as in Fortran
   - Shock speed: `Vs = (v*D² ± j*√(j² + D²)) / (j² + D²)`

2. **For shocks with tangential velocity (vt≠0)**:
   - Verify correct Pons et al. formulation from original paper
   - Ensure consistency between W calculation and shock speed formula

3. **Test**: Run 1D Sod test with MM94 formula and compare accuracy

---

## CODE CHANGES NEEDED

```cpp
// In solve_shock(), line 59-76:

// OPTION 1: Use MM94 formula when vt=0
if (vta < 1e-10) {
    // Pure normal shock - use MM94 formula
    const real Wa = 1.0 / std::sqrt(1.0 - va * va);  // W from normal velocity only
    const real Da = rhoa * Wa;

    const real term = j * j + Da * Da;  // MM94: no (1-v²) factor
    const real sqrt_term = std::sqrt(term);
    const real denom_vs = Da * Da + j * j;

    if (is_left_wave) {
        Vs = (Da * Da * va - j * sqrt_term) / denom_vs;
    } else {
        Vs = (Da * Da * va + j * sqrt_term) / denom_vs;
    }
} else {
    // Shock with tangential component - use Pons formula
    // ... current code ...
}

// OPTION 2: Always use MM94 formula for better accuracy
const real Wa = 1.0 / std::sqrt(1.0 - va * va);  // W from normal velocity
const real Da = rhoa * Wa;
const real term = j * j + Da * Da;  // MM94 formula
// ... rest of shock speed calculation ...
```
