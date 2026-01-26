# SPH/GSPH Relaxation Implementation Report
Generated: 2026-01-24

## Executive Summary

The shamrock codebase implements **TRUE Lane-Emden relaxation** for SPH initial conditions. The method subtracts analytical equilibrium forces from SPH-computed forces to drive particles toward hydrostatic equilibrium. Three implementations exist:

1. **Lane-Emden (n=1.5)** - Polytropic spheres
2. **Isothermal Bonnor-Ebert** - True isothermal self-gravitating spheres
3. **Koyama-Inutsuka (2000)** - Barotropic EOS with K&I cooling tables

---

## 1. Core Relaxation Algorithm

### Method: TRUE LANE-EMDEN Relaxation

**Mathematical Foundation:**

In hydrostatic equilibrium:
```
dP/dr = -ρ GM(r)/r²
```

This gives:
```
a_hydro = -(1/ρ) dP/dr     [pressure gradient, OUTWARD]
a_grav  = -GM(r)/r²        [gravity, INWARD]
a_total = 0                [equilibrium]
```

**Relaxation Strategy** (? INFERRED from code patterns):
```cpp
// Net acceleration drives particles toward equilibrium
a_net = a_SPH_total - a_analytical_pressure
      = (a_grav + a_SPH_pressure) - a_analytical_pressure
      = a_grav + (a_SPH_pressure - a_analytical_pressure)
```

**When it works:** SPH density ≈ analytical → small corrections → convergence  
**When it fails:** SPH density ≫ analytical → huge corrections → divergence

### Perturbative Nature (✓ VERIFIED)

From `/docs/RELAXATION_FORCE_REVIEW.md` (line 138):
> "The relaxation method uses a **perturbative approach**:
> Correction = (SPH_value - Analytical_value)
> This works ONLY if |Correction| << Analytical_value."

**Critical requirement:** Good initial conditions (within ~30% of target density)

**Solution:** Glass-making initialization
- Creates uniform sphere → radius mapping to match cumulative mass
- Initial density error: 1000% → 27% (40× improvement!)
- Relaxation then converges successfully

---

## 2. Relaxation Loop Implementation

**Location:** `/src/solver.cpp`, lines 1970-2025 (✓ VERIFIED)

```cpp
for(int step = start_step; step < target_step; ++step) {
    // 1. Rebuild tree (particles move during relaxation)
    tree->resize(num_p);
    tree->make(particles, num_p);
    
    // 2. Calculate SPH forces (pressure + gravity)
    m_pre->calculation(m_sim);          // Density, pressure
    m_fforce->calculation(m_sim);       // SPH pressure gradient
    if(m_param->gravity.is_valid) {
        m_gforce->calculation(m_sim);   // Gravitational force
    }
    
    // 3. Apply relaxation (subtract analytical forces)
    m_lane_emden_relax->apply_relaxation(m_sim, 0.0);
    
    // 4. Calculate adaptive timestep (CFL condition)
    m_timestep->calculation(m_sim);
    real dt_relax = m_sim->get_dt();
    
    // 5. Integrate positions (QUASI-STATIC method)
    #pragma omp parallel for
    for(int i = 0; i < num_p; ++i) {
        if(particles[i].is_ghost) continue;
        
        // Zero velocities (quasi-static approximation)
        particles[i].vel = 0.0;
        
        // Kinematic integration: Δx = ½at²
        const real half_dt2 = 0.5 * dt_relax * dt_relax;
        particles[i].pos += particles[i].acc * half_dt2;
        
        periodic->apply(particles[i].pos);
    }
    
    // 6. Remove escaping particles
    m_lane_emden_relax->remove_escaping_particles(m_sim, 1.1);
    
    accumulated_time += dt_relax;
}
```

### Damping Strategy (✓ VERIFIED)

**Velocity damping:** Velocities zeroed EVERY step (line 2007-2011)
- Method: Quasi-static relaxation
- NOT gradual damping (0.9 factor unused)
- Particles move purely in direction of net force

**Position update:** `Δx = ½at²` (line 2014-2018)
- Standard kinematic integration with v₀=0
- NOT steepest descent (`Δx = a·dt` mentioned in comment but not used)

---

## 3. Force Computation During Relaxation

### Lane-Emden n=1.5 Implementation

**Location:** `/src/relaxation/lane_emden_relaxation.cpp`, lines 35-120 (✓ VERIFIED)

```cpp
vec_t LaneEmdenRelaxation::compute_relaxation_force(const SPHParticle& p) const
{
    // Spherical radius
    const real r = std::abs(p.pos);
    
    // Dimensionless coordinate: ξ = r / α
    const real xi = r / m_params.alpha_scaling;
    
    // Get Lane-Emden solution from pre-computed table
    const real theta = m_data.get_theta(xi);      // θ(ξ)
    const real dtheta = m_data.dtheta_dxi(xi);    // dθ/dξ
    
    // Analytical pressure gradient acceleration
    // For polytrope P = K ρ^γ with ρ = ρ_c θ^n:
    // a_r = -(1/ρ) dP/dr = -K γ ρ_c^(γ-1) (n/α) dθ/dξ
    //
    // For n=1.5, γ=5/3: exponent θ^(nγ-n-1) = θ^0 = 1
    // Therefore: a_r = -K γ ρ_c^(γ-1) (n/α) dθ/dξ
    
    const real n = 1.5;
    const real prefactor = m_params.K * m_params.gamma 
                         * pow(m_params.rho_center, m_params.gamma - 1.0) 
                         * n / m_params.alpha_scaling;
    
    // Since dθ/dξ < 0 in interior, a_r > 0 (OUTWARD pressure support)
    const real a_r = -prefactor * dtheta;
    
    // Convert to Cartesian vector
    const real r_inv = 1.0 / r;
    vec_t force;
    force[0] = a_r * p.pos[0] * r_inv;
    force[1] = a_r * p.pos[1] * r_inv;
    force[2] = a_r * p.pos[2] * r_inv;  // (if DIM==3)
    
    return force;
}
```

**Application** (line 136-148):
```cpp
void LaneEmdenRelaxation::apply_relaxation(std::shared_ptr<Simulation> sim, real damping_factor)
{
    auto& particles = sim->get_particles();
    
    #pragma omp parallel for
    for(int i = 0; i < num_p; ++i) {
        vec_t relax_acc = compute_relaxation_force(particles[i]);
        
        // SUBTRACT analytical pressure gradient from SPH acceleration
        particles[i].acc -= relax_acc;
    }
    
    // NOTE: Velocities zeroed in solver loop, NOT here
}
```

### Sign Verification (✓ VERIFIED)

From `/docs/RELAXATION_FORCE_REVIEW.md`:
- **Interior:** `dθ/dξ < 0` → `a_r = -prefactor × (negative) = POSITIVE` → **OUTWARD** ✓
- **Mathematical derivation:** Lines 12-68 show rigorous proof
- **Code implementation:** Lines 53-91 confirm sign is CORRECT

---

## 4. Isothermal Bonnor-Ebert Relaxation

**Location:** `/src/relaxation/isothermal_relaxation.cpp` (✓ VERIFIED)

### Isothermal Lane-Emden Equation

```
d²ψ/dξ² + (2/ξ)dψ/dξ = exp(-ψ)
```

where:
- `ξ = r / r_0` (dimensionless radius)
- `r_0 = c_s / sqrt(4πGρ_c)` (scale length)
- `ρ(ξ) = ρ_c × exp(-ψ(ξ))` (TRUE BE density profile)

**Solver:** RK4 integration (lines 37-86)

```cpp
// Initial conditions (near origin): ψ ~ ξ²/6
xi[0] = 1e-6;
psi[0] = xi[0] * xi[0] / 6.0;
dpsi[0] = xi[0] / 3.0;

// RK4 integration of coupled system
for (int i = 1; i < n_points; ++i) {
    // f1 = dψ/dξ = z
    // f2 = dz/dξ = exp(-ψ) - 2z/ξ
    
    // ... RK4 coefficients k1, k2, k3, k4 ...
    
    psi[i] = y + (k1_y + 2*k2_y + 2*k3_y + k4_y) / 6.0;
    dpsi[i] = z + (k1_z + 2*k2_z + 2*k3_z + k4_z) / 6.0;
}
```

### Equilibrium Force (lines 248-287)

```cpp
vec_t IsothermalRelaxation::compute_relaxation_force(const SPHParticle& p) const
{
    // From hydrostatic equilibrium: dP/dr = -ρ g
    // For isothermal: P = ρ c_s², so c_s² d(ρ)/dr = -ρ g
    // Therefore: a_pressure = -c_s² d(ln ρ)/dr
    //
    // Since ρ = ρ_c exp(-ψ), we have ln(ρ) = ln(ρ_c) - ψ
    // So d(ln ρ)/dr = -dψ/dr = -(dψ/dξ)(dξ/dr) = -(dψ/dξ)/r_0
    //
    // Therefore: a_pressure = c_s² (dψ/dξ) / r_0
    
    real r = std::abs(p.pos);
    real xi = r / m_r_0;
    
    // Get dψ/dξ from Lane-Emden solution
    real dpsi_dxi = interpolate_dpsi(xi);
    
    // Equilibrium acceleration: a_eq = c_s² (dψ/dξ) / r_0
    real a_r_mag = m_c_s_sq * dpsi_dxi / m_r_0;
    
    vec_t r_hat = p.pos / r;
    return r_hat * a_r_mag;
}
```

**Application:** Same subtraction pattern (lines 289-328)

```cpp
void IsothermalRelaxation::apply_relaxation(std::shared_ptr<Simulation> sim, real damping_factor)
{
    auto& particles = sim->get_particles();
    
    #pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        if (p_i.is_ghost) continue;
        
        vec_t relax_acc = compute_relaxation_force(p_i);
        
        // SUBTRACT analytical pressure gradient from SPH acceleration
        p_i.acc -= relax_acc;
    }
    
    // Velocities zeroed in solver loop (quasi-static)
}
```

---

## 5. Bonnor-Ebert Initialization

**Location:** `/src/sample/isothermal_bonnor_ebert.cpp` (✓ VERIFIED)

### Physical Model (lines 1-20)

```
Isothermal Lane-Emden: (1/ξ²) d/dξ(ξ² dψ/dξ) = exp(-ψ)
where:
  ξ = r / r_0          (dimensionless radius)
  r_0 = c_s / sqrt(4πGρ_c)  (scale length)
  ρ = ρ_c exp(-ψ)      (TRUE hydrostatic density)

Critical radius: ξ_crit ≈ 6.45
  ξ < ξ_crit: stable equilibrium
  ξ > ξ_crit: gravitationally unstable (collapse)
```

### Particle Placement Strategy (lines 277-300, ? INFERRED)

**Method:** Cumulative mass mapping

1. **Solve Lane-Emden equation** → get `ρ(ξ)`, `dψ/dξ(ξ)`
2. **Compute enclosed mass:**
   ```cpp
   M(r) = 4π ρ_c r_0³ ξ² (dψ/dξ)
   ```
3. **Build cumulative mass table:**
   ```cpp
   M_cumulative[i] = M_enc[i] / M_total
   ```
4. **Glass-making initialization** (likely):
   - Create uniform sphere with equal particle spacing
   - Map radii to match analytical cumulative mass
   - Ensures SPH density ≈ analytical density from start

### Key Parameters (lines 132-165)

```cpp
// Input
int N;                    // Grid resolution (N³ particles)
real M_cloud;             // Cloud mass [M_☉]
real T_cloud;             // Temperature [K]
real mu;                  // Mean molecular weight (1.27 atomic, 2.33 molecular)
real xi_s;                // Truncation radius (default 6.0, critical 6.45)

// Derived
real c_s = sqrt(k_B T / (mu m_H));            // Sound speed [km/s]
real r_0 = c_s / sqrt(4π G ρ_c);              // Scale length [pc]
real R_cloud = xi_s × r_0;                    // Physical radius [pc]
real rho_edge = rho_c × exp(-ψ(xi_s));        // Edge density
real P_ext = n_edge × T_cloud;                // External pressure [K cm⁻³]
```

---

## 6. Koyama-Inutsuka (2000) Relaxation

**Location:** `/src/relaxation/koyama_inutsuka_relaxation.cpp` (✓ VERIFIED)

### Barotropic EOS (lines 69-108)

```cpp
real KoyamaInutsukaRelaxation::get_c_eff_squared(real rho) const
{
    // Effective sound speed: c_eff² = dP/dρ
    // For barotropic P = n k_B T_eq(n):
    // c_eff² = (k_B T / m_n) × [1 + d ln T / d ln n]
    
    real n_H = rho * m_params.density_to_n;
    real T_eq = m_ki_cooling.equilibrium_temperature(n_H, m_params.N_H);
    
    // Numerical derivative: d ln T / d ln n
    const real eps = 0.02;
    real n_lo = n_H * (1.0 - eps);
    real n_hi = n_H * (1.0 + eps);
    real T_lo = m_ki_cooling.equilibrium_temperature(n_lo, m_params.N_H);
    real T_hi = m_ki_cooling.equilibrium_temperature(n_hi, m_params.N_H);
    
    real dlnT_dlnn = (log(T_hi) - log(T_lo)) / (log(n_hi) - log(n_lo));
    
    real c_eff_sq_cgs = (k_B * T_eq / m_n) * (1.0 + dlnT_dlnn);
    
    // Convert to code units: [cm/s]² → [km/s]²
    return c_eff_sq_cgs / 1.0e10;
}
```

### Hydrostatic Integration (lines 110-253)

```cpp
void KoyamaInutsukaRelaxation::compute_equilibrium_profile()
{
    // Integrate hydrostatic equilibrium ODEs:
    //   dρ/dr = -ρ G M(r) / (r² c_eff²)
    //   dM/dr = 4π r² ρ
    // Until: P(ρ) ≤ P_ext (pressure truncation)
    
    real r = 1e-6 * R_cloud;           // Start near center
    real rho = rho_center;
    real M = (4/3)π ρ r³;
    
    for (int i = 0; i < n_points_max && r < r_max; ++i) {
        // Current state
        real T = get_T_eq_from_rho(rho);
        real P_cgs = n_H × k_B × T;
        
        // Store profile point
        m_r_table.push_back(r);
        m_rho_table.push_back(rho);
        m_P_table.push_back(P_cgs / k_B);
        m_M_table.push_back(M);
        m_T_table.push_back(T);
        
        // Check truncation
        if (P_cgs <= P_ext_cgs) {
            truncated = true;
            break;
        }
        
        // Get effective sound speed
        real c_eff_sq = get_c_eff_squared(rho);
        
        // Hydrostatic gradient
        real drho_dr = -rho × G × M / (r² × c_eff_sq);
        
        // RK4 integration step
        // ... k1, k2, k3, k4 coefficients ...
        
        rho = rho + (dr/6) × (k1_rho + 2×k2_rho + 2×k3_rho + k4_rho);
        M = M + (dr/6) × (k1_M + 2×k2_M + 2×k3_M + k4_M);
        r = r + dr;
    }
}
```

### Relaxation Force (lines 462-511)

```cpp
vec_t KoyamaInutsukaRelaxation::compute_relaxation_force(const SPHParticle& p) const
{
    real r = std::abs(p.pos);
    
    // Get equilibrium profile at this radius
    real rho_eq = get_rho_eq(r);
    real drho_dr = get_drho_dr(r);
    real c_eff_sq = get_c_eff_squared(rho_eq);
    
    // Analytical pressure gradient acceleration:
    // a_r = -(1/ρ) dP/dr
    // For barotropic: dP/dr = c_eff² dρ/dr
    // Therefore: a_r = -c_eff² (1/ρ) dρ/dr
    
    real a_r = -c_eff_sq * drho_dr / rho_eq;
    
    const real r_inv = 1.0 / r;
    vec_t force;
    force[0] = a_r * p.pos[0] * r_inv;
    force[1] = a_r * p.pos[1] * r_inv;
    force[2] = a_r * p.pos[2] * r_inv;
    
    return force;
}
```

**Application:** Same subtraction pattern (lines 513-546)

---

## 7. Convergence Parameters

### Relaxation Control (? INFERRED from solver parameters)

**From `/src/solver.cpp`:**

```cpp
m_use_relaxation = input.get<bool>("useRelaxation", false);
m_relaxation_steps = input.get<int>("relaxationSteps", 0);
m_relaxation_output_freq = input.get<int>("relaxationOutputFreq", 10);
m_relaxation_only = input.get<bool>("relaxationOnly", false);

// GLASS pre-relaxation (optional)
m_use_glass_relaxation = input.get<bool>("useGlassRelaxation", false);
m_glass_relaxation_steps = input.get<int>("glassRelaxationSteps", 5000);
```

### Convergence Monitoring (✓ VERIFIED from solver loop)

**From `/src/solver.cpp`, lines 2057-2062:**

```cpp
// Calculate max acceleration (convergence metric)
real max_acc = 0.0;
for(const auto& pi : particles) {
    real a = sqrt(pi.acc·pi.acc);
    max_acc = max(max_acc, a);
}
```

**Expected behavior:**
- `max_acc` should decrease as particles settle
- No explicit convergence threshold (runs fixed number of steps)
- User monitors `max_acc` output to decide when to stop

### Particle Removal (✓ VERIFIED)

```cpp
int remove_escaping_particles(Simulation* sim, real tolerance_factor = 1.1)
{
    // Remove particles beyond tolerance_factor × R_cloud
    // Default: 1.1 (10% buffer zone)
    
    auto new_end = std::remove_if(particles.begin(), particles.end(),
        [R_max](const SPHParticle& p) {
            return std::abs(p.pos) > R_max;
        });
    
    particles.erase(new_end, particles.end());
    return removed_count;
}
```

---

## 8. Key Algorithms Summary

### Relaxation Method

| Component | Algorithm | Implementation |
|-----------|-----------|----------------|
| **Force computation** | Analytical pressure gradient from Lane-Emden solution | `compute_relaxation_force()` |
| **Force application** | Subtract analytical from SPH: `a_net = a_SPH - a_analytical` | `apply_relaxation()` |
| **Damping** | Zero velocities every step (quasi-static) | Solver loop line 2007 |
| **Integration** | Kinematic: `Δx = ½at²` with `v₀=0` | Solver loop line 2014 |
| **Timestep** | Adaptive CFL condition from SPH | `m_timestep->calculation()` |
| **Convergence** | Monitor `max_acc` (no auto-stop) | Solver loop line 2057 |

### Isothermal Profile Computation

| Step | Method | Lines |
|------|--------|-------|
| 1. Solve ODE | RK4 integration of `d²ψ/dξ² + (2/ξ)dψ/dξ = exp(-ψ)` | isothermal_relaxation.cpp:37-86 |
| 2. Interpolate | Binary search + linear interpolation | isothermal_relaxation.cpp:175-225 |
| 3. Compute density | `ρ(r) = ρ_c exp(-ψ(r/r_0))` | isothermal_relaxation.cpp:227-240 |
| 4. Compute force | `a_eq = (c_s²/r_0)(dψ/dξ)r̂` | isothermal_relaxation.cpp:248-287 |

### K&I Profile Computation

| Step | Method | Lines |
|------|--------|-------|
| 1. Get T_eq(n) | Lookup from K&I 2000 tables | koyama_inutsuka_relaxation.cpp:58-67 |
| 2. Compute c_eff² | `c_eff² = (k_B T/m_n)(1 + d ln T/d ln n)` | koyama_inutsuka_relaxation.cpp:69-108 |
| 3. Integrate ODE | RK4: `dρ/dr = -ρGM/(r²c_eff²)` | koyama_inutsuka_relaxation.cpp:110-253 |
| 4. Interpolate | Binary search + linear | koyama_inutsuka_relaxation.cpp:410-435 |
| 5. Compute force | `a_r = -c_eff²(1/ρ)dρ/dr` | koyama_inutsuka_relaxation.cpp:462-511 |

---

## 9. Critical Insights

### Why "TRUE LANE-EMDEN"? (✓ VERIFIED)

From file headers and comments:
- **TRUE** means using exact analytical solution, NOT approximations
- Solves full Lane-Emden ODE numerically (RK4)
- Subtracts analytical forces directly (not fitting parameters)

### Why It Can Fail (✓ VERIFIED)

From `/docs/RELAXATION_FORCE_REVIEW.md`:

**Equal-mass shells initialization:**
- Initial density error: 1000% (ρ_SPH = 14.36, ρ_analytic = 1.43)
- Initial pressure error: 5400% (P ∝ ρ^(5/3))
- Huge outward force → particles blown away → complete failure

**Glass-making initialization:**
- Initial density error: 27% (ρ_SPH ≈ 1.17, ρ_analytic ≈ 0.92)
- Initial pressure error: 75%
- Small corrections → successful convergence

### Quasi-Static Assumption (✓ VERIFIED)

Velocities zeroed every step means:
- Particles have NO kinetic energy
- System evolves through sequence of static configurations
- Valid when relaxation timescale ≫ sound crossing time
- Faster convergence than full dynamic simulation

---

## 10. File Locations

### Core Relaxation Code

| File | Purpose | Lines | Status |
|------|---------|-------|--------|
| `/src/relaxation/lane_emden_relaxation.cpp` | Lane-Emden n=1.5 relaxation | 180 | ✓ VERIFIED |
| `/src/relaxation/isothermal_relaxation.cpp` | Isothermal BE relaxation | 366 | ✓ VERIFIED |
| `/src/relaxation/koyama_inutsuka_relaxation.cpp` | K&I 2000 barotropic relaxation | 658 | ✓ VERIFIED |
| `/include/relaxation/*.hpp` | Header interfaces | - | ✓ VERIFIED |

### Sample Generators

| File | Purpose | Lines | Status |
|------|---------|-------|--------|
| `/src/sample/lane_emden.cpp` | Lane-Emden IC generator | ~20k | ? INFERRED |
| `/src/sample/isothermal_bonnor_ebert.cpp` | Isothermal BE IC | 32447 | ✓ VERIFIED (partial) |
| `/src/sample/bonnor_ebert_ki2000.cpp` | K&I BE IC | 15955 | ✓ VERIFIED (partial) |

### Solver Integration

| File | Purpose | Lines | Status |
|------|---------|-------|--------|
| `/src/solver.cpp` | Main relaxation loop | 1970-2025 | ✓ VERIFIED |
| `/docs/RELAXATION_FORCE_REVIEW.md` | Algorithm review & validation | 214 | ✓ VERIFIED |

---

## 11. Open Questions

1. **Glass-making details** (? INFERRED):
   - Exact algorithm for uniform sphere generation
   - Radius mapping procedure
   - Located in sample generators (not yet read fully)

2. **Convergence criteria** (? INFERRED):
   - No automatic stopping condition found
   - User must monitor `max_acc` and decide when converged
   - Typical values for "converged" state?

3. **GLASS pre-relaxation** (? INFERRED from parameters):
   - `useGlassRelaxation`, `glassRelaxationSteps` mentioned
   - Implementation not yet traced
   - Likely in sample generators

4. **Polytropic slab relaxation** (not explored):
   - `/src/relaxation/polytropic_slab*.cpp` exist
   - 1D/2D variants for planar geometry

---

## 12. Recommended Next Steps

1. **Read glass-making code** in sample generators
2. **Trace GLASS pre-relaxation** implementation
3. **Check polytropic slab** relaxation (1D/2D variants)
4. **Find typical convergence thresholds** from example configs
5. **Understand profile update** mechanism (K&I has `update_profile_from_sph()`)

---

## Verification Markers

- ✓ VERIFIED: Read source code directly, traced logic
- ? INFERRED: Based on code patterns, comments, or indirect evidence
- ✗ UNCERTAIN: Not yet investigated

All file paths are absolute from project root:
`/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/`
