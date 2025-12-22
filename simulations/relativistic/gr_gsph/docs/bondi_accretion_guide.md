# Bondi Accretion in Schwarzschild Spacetime

## A Beginner's Guide to GR-GSPH Benchmarks

This guide explains the physics and mathematics behind the Bondi accretion benchmarks, based on:
- **Liptai & Price (2019)** MNRAS 485, 819 - "General relativistic smoothed particle hydrodynamics"
- **Hawley, Smarr & Wilson (1984)** - Original analytic solutions

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [The Schwarzschild Metric](#2-the-schwarzschild-metric)
3. [Geodesic Flow (Pressure-less)](#3-geodesic-flow-pressure-less)
4. [Sonic Point Flow (With Pressure)](#4-sonic-point-flow-with-pressure)
5. [Running the Simulations](#5-running-the-simulations)
6. [Understanding the Results](#6-understanding-the-results)

---

## 1. Introduction

### What is Bondi Accretion?

Bondi accretion describes gas falling radially onto a compact object (like a black hole). It's one of the simplest models of accretion and serves as an excellent test for numerical relativity codes.

```
       Gas flowing inward
            ↓  ↓  ↓
    ══════════════════════
         ↓  ↓  ↓  ↓
           ●────────→ r
        Black Hole
        (mass M)
```

### Why is it Important?

1. **Analytic solutions exist** - We can compare simulations to exact answers
2. **Tests GR effects** - Strong gravity near the black hole
3. **Tests hydrodynamics** - Pressure, density, velocity evolution
4. **Benchmark standard** - Used by many GR codes for validation

### Two Test Cases

| Test | Physics | Key Feature |
|------|---------|-------------|
| **Geodesic Flow** | No pressure (u ≪ 1) | Tests gravity + advection |
| **Sonic Point Flow** | With pressure | Tests full hydrodynamics |

---

## 2. The Schwarzschild Metric

### What is a Metric?

A metric tells us how to measure distances in spacetime. In flat space:
```
ds² = -dt² + dx² + dy² + dz²    (Minkowski metric, c=1)
```

Near a black hole, spacetime is **curved**, and distances are measured differently.

### The Schwarzschild Solution

For a non-rotating black hole of mass M, the metric is:

```
ds² = -(1 - 2M/r)dt² + (1 - 2M/r)⁻¹dr² + r²dΩ²
```

where:
- `ds²` = spacetime interval (proper distance squared)
- `dt` = coordinate time interval
- `dr` = radial coordinate interval
- `dΩ² = dθ² + sin²θ dφ²` = angular part
- `r` = radial coordinate (NOT physical distance!)
- `M` = black hole mass (in geometric units where G=c=1)

### Key Quantity: The Lapse Function

The **lapse function** α measures how time flows:

```
α² = 1 - 2M/r
```

- At r → ∞: α → 1 (normal time flow)
- At r = 2M: α = 0 (event horizon - time "stops" for distant observer)
- For r < 2M: Inside the black hole (not simulated)

### The Event Horizon

The **event horizon** is at r = 2M (called the Schwarzschild radius):

```
r_s = 2M = 2GM/c²    (in physical units)
```

For M = 1 (our simulations): r_s = 2

**Important:** Our simulation domain is r ∈ [2.5, 18], safely outside the horizon.

---

## 3. Geodesic Flow (Pressure-less)

### Physical Setup

Imagine dust (pressure-less gas) falling radially into a black hole:
- No pressure forces (P = 0)
- Gas follows geodesics (free-fall paths)
- Internal energy is negligible (u ≪ 1)

### Step-by-Step Derivation

#### Step 1: Conservation of Energy

For a particle falling from rest at infinity, energy conservation gives:

```
E = (1 - 2M/r) × dt/dτ = 1
```

where τ is proper time (time measured by the falling particle).

Solving for dt/dτ:
```
dt/dτ = 1/(1 - 2M/r) = 1/α²
```

#### Step 2: Radial Velocity

From the metric with ds² = -dτ² (for a massive particle):

```
-1 = -(1 - 2M/r)(dt/dτ)² + (1 - 2M/r)⁻¹(dr/dτ)²
```

Substituting dt/dτ:
```
-1 = -1/(1 - 2M/r) + (dr/dτ)²/(1 - 2M/r)
```

Solving for dr/dτ (the **proper velocity**):
```
(dr/dτ)² = 2M/r
dr/dτ = -√(2M/r)    (negative = falling inward)
```

#### Step 3: Coordinate Velocity

The **coordinate velocity** (what a distant observer measures) is:

```
v^r = dr/dt = (dr/dτ)/(dt/dτ) = -√(2M/r) × (1 - 2M/r)
```

**This is Equation (941) in Liptai & Price (2019):**
```
┌─────────────────────────────────────┐
│  v^r = -√(2M/r) × (1 - 2M/r)       │
└─────────────────────────────────────┘
```

#### Step 4: Density Profile

Mass conservation (continuity equation) in spherical symmetry:

```
∂ρ/∂t + (1/r²) ∂(r² ρ v^r)/∂r = 0
```

For steady-state (∂/∂t = 0):
```
r² ρ v^r = constant
```

Substituting v^r and solving:
```
┌─────────────────────────────────────┐
│  ρ = ρ₀/r² × √(r/2M)               │
└─────────────────────────────────────┘
```

This is Equation (942) - density increases as gas falls inward.

#### Step 5: Thermal Energy (Adiabatic)

For adiabatic compression with γ = 5/3:

```
┌─────────────────────────────────────┐
│  u = u₀ × (r² √(2M/r))^(1-γ)       │
└─────────────────────────────────────┘
```

Since γ = 5/3 > 1, the exponent (1-γ) = -2/3 is negative, so:
- u **increases** as r decreases (compression heats the gas)
- We use u₀ = 10⁻⁹ to keep u ≪ 1 (pressure-less assumption)

### Summary: Geodesic Flow Equations

| Quantity | Formula | At r=2.5M | At r=18M |
|----------|---------|-----------|----------|
| Velocity v^r | -√(2M/r)(1-2M/r) | -0.18 | -0.30 |
| Density ρ | ρ₀/r² √(r/2M) | 0.18 | 0.008 |
| Energy u | u₀(r²√(2M/r))^(-2/3) | 3×10⁻¹⁰ | 2×10⁻¹¹ |

---

## 4. Sonic Point Flow (With Pressure)

### Physical Setup

Now include pressure effects:
- Gas has significant thermal energy
- Pressure gradient opposes gravity
- Flow passes through a **sonic point** where v = c_s

### What is a Sonic Point?

The **sonic point** is where the flow velocity equals the sound speed:

```
|v| = c_s    at r = r_c (critical radius)
```

- For r > r_c: **Subsonic** (|v| < c_s) - pressure dominates
- For r < r_c: **Supersonic** (|v| > c_s) - gravity dominates

For γ = 5/3: **r_c = 8M**

### Step-by-Step Derivation

#### Step 1: Define Temperature Variable

Instead of internal energy u, we use a "temperature" variable:

```
T ≡ P/ρ = (γ-1)u
```

This is NOT the physical temperature, but a convenient variable.

#### Step 2: Polytropic Index

Define:
```
n = 1/(γ-1)
```

For γ = 5/3: n = 3/2 = 1.5

#### Step 3: Four-Velocity

The four-velocity U^μ satisfies:
- U^μ U_μ = -1 (normalization)
- U^r = dr/dτ (radial component)

The solution has the form:

```
┌─────────────────────────────────────┐
│  U^r = C₁ / (r² T^n)               │
└─────────────────────────────────────┘
```

where C₁ is a constant determined by the sonic point.

#### Step 4: Density and Energy

From the equations of state:

```
┌─────────────────────────────────────┐
│  ρ = K₀ T^n                        │
│  u = n × T                          │
└─────────────────────────────────────┘
```

where K₀ = 1 is the entropy normalization.

#### Step 5: The Implicit Equation for T(r)

T(r) is found by solving:

```
┌───────────────────────────────────────────────────────┐
│  C₂ = [1 + (n+1)T]² × [1 - 2M/r + C₁²/(r⁴ T^(2n))]  │
└───────────────────────────────────────────────────────┘
```

This is an **implicit equation** - T appears on both sides!

#### Step 6: Critical Point Conditions

At r = r_c, the flow is sonic. This gives us C₁ and C₂:

**Critical velocity:**
```
U_c = √(M / 2r_c)
```

**Critical coordinate velocity:**
```
v_c = U_c / √(1 - 3U_c²)
```

**Critical temperature:**
```
T_c = v_c² n / [(1+n)(1 - n v_c²)]
```

**Constants:**
```
C₁ = U_c × r_c² × T_c^n
C₂ = [1 + (1+n)T_c]² × [1 - 2M/r_c + C₁²/(r_c⁴ T_c^(2n))]
```

#### Step 7: Convert to Coordinate Velocity

From U^r, compute the coordinate velocity v^r:

```
U⁰ = √(1 - 2M/r + (U^r)²) / (1 - 2M/r)

v^r = -U^r / U⁰    (negative for infall)
```

### Numerical Values for r_c = 8M

| Quantity | Value |
|----------|-------|
| U_c | 0.25 |
| v_c | 0.277 |
| T_c | 0.0522 |
| C₁ | 0.191 |
| C₂ | 1.038 |

### Solving the Implicit Equation

The implicit equation must be solved numerically:

1. **For r > r_c (subsonic):** T < T_c, solution exists in [0, 2T_c]
2. **For r < r_c (supersonic):** May need to extrapolate with T = T_c

```python
# Python approach (bisection/brentq)
from scipy.optimize import brentq

def equation(T, r, M, n, C1, C2):
    term1 = (1 + (n+1)*T)**2
    term2 = 1 - 2*M/r + C1**2 / (r**4 * T**(2*n))
    return term1 * term2 - C2

T = brentq(equation, 1e-10, 2*T_c, args=(r, M, n, C1, C2))
```

### Summary: Sonic Point Equations

| Quantity | Formula |
|----------|---------|
| Temperature | Solve implicit equation |
| Density | ρ = K₀ T^n |
| 4-velocity | U^r = C₁/(r² T^n) |
| Coord. velocity | v^r = -U^r/U⁰ |
| Internal energy | u = n T |

---

## 5. Running the Simulations

### Prerequisites

```bash
# Build the SPH code
cd /path/to/sphcode
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
```

### Run Geodesic Flow

```bash
./build/sph simulations/relativistic/gr_gsph/config/presets/gr_bondi_geodesic.json
```

Output: `simulations/relativistic/gr_gsph/results/bondi_geodesic/`

### Run Sonic Point Flow

```bash
./build/sph simulations/relativistic/gr_gsph/config/presets/gr_bondi_sonic.json
```

Output: `simulations/relativistic/gr_gsph/results/bondi_sonic/`

### Generate Benchmark Plots

```bash
python3 simulations/relativistic/gr_gsph/scripts/bondi_accretion_benchmark.py
```

Output: `simulations/relativistic/gr_gsph/results/bondi_*.png`

### Configuration Parameters

Key parameters in the JSON config files:

| Parameter | Geodesic | Sonic | Description |
|-----------|----------|-------|-------------|
| `testType` | "geodesic" | "sonic" | Which flow to simulate |
| `blackHoleMass` | 1.0 | 1.0 | M in geometric units |
| `rIn` | 2.5 | 2.5 | Inner boundary (> 2M) |
| `rOut` | 18.0 | 18.0 | Outer boundary |
| `N` | 500 | 500 | Number of particles |
| `gamma` | 5/3 | 5/3 | Adiabatic index |
| `avAlpha` | 0.0 | 0.0 | Artificial viscosity |

---

## 6. Understanding the Results

### Reading the Output Plots

The benchmark generates 6-panel plots:

```
┌─────────────┬─────────────┬─────────────┐
│  Geodesic   │  Geodesic   │  Geodesic   │
│  Density    │  Velocity   │  Energy     │
├─────────────┼─────────────┼─────────────┤
│ Sonic Point │ Sonic Point │ Sonic Point │
│  Density    │  Velocity   │  Energy     │
└─────────────┴─────────────┴─────────────┘
```

- **Green lines:** Analytic solution
- **Blue dots:** C++ GR-GSPH simulation
- **Red dashed line (sonic only):** Sonic point at r = 8M

### What Good Agreement Looks Like

**Geodesic Flow:**
- Blue dots should follow green curves closely
- Small deviations at inner boundary (fewer neighbors)
- Smooth profiles throughout domain

**Sonic Point Flow:**
- Blue dots follow green curves in subsonic region (r > 8M)
- Some oscillations near sonic point (r ≈ 8M)
- Transition region may show numerical diffusion

### Common Issues and Causes

| Issue | Possible Cause |
|-------|----------------|
| Density drop at inner boundary | Fewer SPH neighbors |
| Velocity noise | Need more particles or smaller timestep |
| Energy spike | Artificial viscosity heating |
| Mismatch at large r | Boundary condition effects |

### Physical Interpretation

**Geodesic Flow (Top Row):**
- Density increases toward black hole (ρ ∝ r^(-3/2))
- Velocity has maximum at r ≈ 6M, then decreases near horizon
- Energy increases due to adiabatic compression

**Sonic Point Flow (Bottom Row):**
- Density roughly constant (because T also varies)
- Velocity increases smoothly through sonic point
- Energy peaks near sonic point, then decreases

---

## References

1. **Liptai & Price (2019)** - "General relativistic smoothed particle hydrodynamics"
   MNRAS 485, 819-842. [arXiv:1901.08064](https://arxiv.org/abs/1901.08064)

2. **Hawley, Smarr & Wilson (1984)** - "A numerical study of nonspherical black hole accretion"
   ApJ 277, 296-311.

3. **Michel (1972)** - "Accretion of matter by condensed objects"
   Astrophys. Space Sci. 15, 153-160.

---

## Appendix: Unit System

Our simulations use **geometric units** where:

| Quantity | Code Value | Physical Meaning |
|----------|------------|------------------|
| G | 1 | Gravitational constant |
| c | 1 | Speed of light |
| M | 1 | Black hole mass |

To convert to physical units:
- Length: L = GM/c²
- Time: T = GM/c³
- Velocity: V = c

For a solar mass black hole (M = M_☉):
- L ≈ 1.5 km
- T ≈ 5 μs
- r_s = 2M ≈ 3 km

---

*Last updated: December 2024*
*Part of the SPHCode GR-GSPH benchmarks*
