# Technical Specifications - SPH Code Simulation

## 1. System Architecture

### 1.1 Overview
The SPH simulation code is a high-performance computational fluid dynamics application designed for astrophysical simulations. It uses Smoothed Particle Hydrodynamics (SPH) to model fluid behavior through discrete particles.

### 1.2 Core Components

```
┌─────────────────────────────────────────────────────────────┐
│                     Main Simulation Loop                     │
├─────────────────────────────────────────────────────────────┤
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────────────┐  │
│  │   Config    │  │  Particle   │  │    Time Integrator  │  │
│  │   Parser    │  │   Manager   │  │   (Leapfrog/RK2)    │  │
│  └─────────────┘  └─────────────┘  └─────────────────────┘  │
├─────────────────────────────────────────────────────────────┤
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────────────┐  │
│  │ Barnes-Hut  │  │    SPH      │  │      Physics        │  │
│  │    Tree     │  │   Kernels   │  │      Modules        │  │
│  └─────────────┘  └─────────────┘  └─────────────────────┘  │
├─────────────────────────────────────────────────────────────┤
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────────────┐  │
│  │    I/O      │  │   Logging   │  │    Diagnostics      │  │
│  │   (HDF5)    │  │             │  │                     │  │
│  └─────────────┘  └─────────────┘  └─────────────────────┘  │
└─────────────────────────────────────────────────────────────┘
```

### 1.3 Dimension Handling
- Compile-time dimension selection via `-DSPH_DIM=1|2|3`
- Template-based Vector<DIM> class for position, velocity, acceleration
- Dimension-specific tree structure (binary tree/quadtree/octree)

## 2. Data Models and Structures

### 2.1 Particle Structure
```cpp
template<int DIM>
struct Particle {
    // Identification
    size_t id;

    // Kinematics
    Vector<DIM> position;
    Vector<DIM> velocity;
    Vector<DIM> acceleration;

    // Thermodynamics
    double mass;
    double density;
    double pressure;
    double internal_energy;
    double entropy;  // For entropy-based formulation

    // SPH Properties
    double smoothing_length;  // h
    int neighbor_count;

    // Time stepping
    double dt;  // Individual timestep (optional)

    // Flags
    uint32_t flags;  // Boundary, active, etc.
};
```

### 2.2 Vector Template
```cpp
template<int DIM>
struct Vector {
    double data[DIM];

    // Operators: +, -, *, /, dot, cross (3D only), magnitude, normalize
    // Element access: operator[]
};
```

### 2.3 Barnes-Hut Tree Node
```cpp
template<int DIM>
struct TreeNode {
    Vector<DIM> center;
    double size;

    // Monopole moment (for gravity)
    double mass;
    Vector<DIM> center_of_mass;

    // Children: 2^DIM children (2 for 1D, 4 for 2D, 8 for 3D)
    std::array<TreeNode*, (1 << DIM)> children;

    // Particles in leaf node
    std::vector<size_t> particle_indices;

    bool is_leaf;
};
```

### 2.4 Simulation State
```cpp
template<int DIM>
struct SimulationState {
    std::vector<Particle<DIM>> particles;
    double time;
    double dt;
    size_t step;

    // Domain
    Vector<DIM> box_min, box_max;
    bool periodic[DIM];

    // Physics parameters
    double gamma;  // Adiabatic index
    double G;      // Gravitational constant
};
```

## 3. SPH Kernel Specifications

### 3.1 Supported Kernels

#### Cubic Spline (M4)
```
W(r,h) = σ/h^d * {
    1 - 3/2*q^2 + 3/4*q^3,           0 ≤ q < 1
    1/4*(2-q)^3,                      1 ≤ q < 2
    0,                                q ≥ 2
}
where q = r/h, σ = {2/3, 10/(7π), 1/π} for d = {1,2,3}
```

#### Wendland C2
```
W(r,h) = σ/h^d * (1-q/2)^4 * (2q+1),  0 ≤ q < 2
where σ = {5/8, 7/(4π), 21/(16π)} for d = {1,2,3}
```

#### Wendland C4
```
W(r,h) = σ/h^d * (1-q/2)^6 * (35/12*q^2 + 3q + 1),  0 ≤ q < 2
where σ = {3/2, 9/(4π), 495/(256π)} for d = {1,2,3}
```

### 3.2 Kernel Requirements
- Normalization: ∫W(r,h)dV = 1
- Symmetry: W(r,h) = W(-r,h)
- Compact support: W(r,h) = 0 for r > κh (typically κ=2)
- Positive definite: W(r,h) ≥ 0

## 4. Physics Modules

### 4.1 Equation of State (EOS)
**Ideal Gas:**
```
P = (γ - 1) * ρ * u
c_s = sqrt(γ * P / ρ)
```
where:
- P = pressure
- γ = adiabatic index (default 5/3 for monatomic, 7/5 for diatomic)
- ρ = density
- u = specific internal energy
- c_s = sound speed

### 4.2 SPH Density Summation
```
ρ_i = Σ_j m_j W(r_ij, h_i)
```

### 4.3 Momentum Equation (Standard SPH)
```
dv_i/dt = -Σ_j m_j (P_i/ρ_i^2 + P_j/ρ_j^2 + Π_ij) ∇W_ij
```
where Π_ij is artificial viscosity term.

### 4.4 Artificial Viscosity (Monaghan-Balsara)
```
Π_ij = {
    (-α*c̄_ij*μ_ij + β*μ_ij^2) / ρ̄_ij,  v_ij · r_ij < 0
    0,                                    otherwise
}
μ_ij = h̄_ij * v_ij · r_ij / (r_ij^2 + η^2)
```
Default: α=1, β=2, η=0.01*h

### 4.5 Energy Equation
```
du_i/dt = 1/2 * Σ_j m_j (P_i/ρ_i^2 + P_j/ρ_j^2 + Π_ij) v_ij · ∇W_ij
```

### 4.6 Self-Gravity (Barnes-Hut)
Opening angle criterion: θ = s/d < θ_crit (default θ_crit = 0.5)
```
a_grav = -G * Σ M_node * r / |r|^3  (monopole approximation)
```

### 4.7 ISM Cooling (Koyama-Inutsuka 2000)
```
Λ(T) = 2×10^-19 * exp(-1.184×10^5 / (T + 1000)) + 2.8×10^-28 * sqrt(T) * exp(-92/T)
Γ = 2×10^-26  (constant heating rate)
du/dt = -(Λ(T)*n - Γ) / ρ
```
where n = number density, T = temperature

## 5. Godunov SPH Specifications

### 5.1 GSPH Formulation
Replace pressure terms with Riemann solution:
```
dv_i/dt = -2 * Σ_j m_j P*_ij / (ρ_i * ρ_j) * ∇W_ij
du_i/dt = 2 * Σ_j m_j P*_ij * v*_ij / (ρ_i * ρ_j) * ∇W_ij
```
where P*_ij, v*_ij are solutions to the Riemann problem.

### 5.2 Riemann Solvers

#### HLL Solver
```
F_HLL = (S_R * F_L - S_L * F_R + S_L * S_R * (U_R - U_L)) / (S_R - S_L)
```
Wave speed estimates: S_L = min(v_L - c_L, v_R - c_R), S_R = max(v_L + c_L, v_R + c_R)

#### HLLC Solver
Includes contact wave S_* for improved contact discontinuity resolution.

## 6. Time Integration

### 6.1 Leapfrog (Kick-Drift-Kick)
```
v^(n+1/2) = v^n + (dt/2) * a^n
x^(n+1) = x^n + dt * v^(n+1/2)
compute a^(n+1)
v^(n+1) = v^(n+1/2) + (dt/2) * a^(n+1)
```

### 6.2 Timestep Criteria
```
dt = C_CFL * min(dt_hydro, dt_force, dt_grav)
dt_hydro = h / (c_s + |v|)
dt_force = sqrt(h / |a|)
dt_grav = sqrt(ε / |a_grav|)  (ε = softening length)
```
Default C_CFL = 0.3

## 7. Configuration System

### 7.1 JSON Configuration Schema
```json
{
    "simulation": {
        "name": "string",
        "dimensions": 2,
        "end_time": 1.0,
        "output_interval": 0.01
    },
    "physics": {
        "gamma": 1.6667,
        "gravity": true,
        "G": 1.0,
        "cooling": false
    },
    "sph": {
        "kernel": "wendland_c2",
        "neighbor_count": 50,
        "alpha_av": 1.0,
        "beta_av": 2.0
    },
    "initial_conditions": {
        "type": "preset",
        "preset": "sod_shock"
    },
    "output": {
        "format": "hdf5",
        "directory": "./output",
        "fields": ["position", "velocity", "density", "pressure"]
    }
}
```

### 7.2 Preset Initial Conditions
- `sod_shock`: 1D Sod shock tube
- `sedov_blast`: 2D/3D Sedov-Taylor blast wave
- `kelvin_helmholtz`: 2D KH instability
- `evrard_collapse`: 3D gravitational collapse
- `imbh_cloud`: IMBH tidal disruption setup

## 8. I/O Specifications

### 8.1 HDF5 Output Format
```
/Header
    /NumParticles
    /Time
    /Dimension
    /BoxSize
/PartType0
    /Coordinates    [N, DIM]
    /Velocities     [N, DIM]
    /Masses         [N]
    /Density        [N]
    /Pressure       [N]
    /InternalEnergy [N]
    /SmoothingLength [N]
    /ParticleIDs    [N]
```

### 8.2 Checkpoint Format
Full simulation state for restart capability.

## 9. Performance Requirements

### 9.1 Target Performance
- 10^5 particles: < 1 minute per 100 timesteps
- 10^6 particles: < 10 minutes per 100 timesteps
- Linear scaling with OpenMP threads (up to 16 cores)

### 9.2 Memory Requirements
- Per particle: ~200 bytes
- Tree overhead: ~50% of particle memory
- 10^6 particles: ~300 MB RAM

### 9.3 Parallelization Strategy
- OpenMP for shared-memory parallelism
- Parallel regions: neighbor search, force calculation, time integration
- Critical sections: tree construction, output

## 10. Build System Requirements

### 10.1 CMake Configuration
```cmake
cmake_minimum_required(VERSION 3.13)
project(sph VERSION 1.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dimension selection
set(SPH_DIM "2" CACHE STRING "Simulation dimension (1, 2, or 3)")
add_compile_definitions(SPH_DIM=${SPH_DIM})

# Dependencies
find_package(OpenMP REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(Boost REQUIRED)
find_package(nlohmann_json REQUIRED)

# Executable
add_executable(sph src/main.cpp ...)
target_link_libraries(sph PRIVATE OpenMP::OpenMP_CXX HDF5::HDF5 ...)
```

### 10.2 Compiler Flags
- Debug: `-g -O0 -DDEBUG`
- Release: `-O3 -march=native -DNDEBUG`
- Profile: `-O2 -g -pg`

## 11. Testing Requirements

### 11.1 Unit Tests (Google Test)
- Vector operations
- Kernel normalization and gradients
- Tree construction and traversal
- EOS calculations
- Riemann solver accuracy

### 11.2 Integration Tests
- Sod shock tube (compare to analytical solution)
- Sedov blast wave (self-similar solution)
- Evrard collapse (energy conservation)

### 11.3 Acceptance Criteria
- Sod shock: position error < 1% at contact
- Sedov: radius error < 5%
- Energy conservation: ΔE/E < 10^-4 per dynamical time

## 12. Relativistic Extensions (Advanced)

### 12.1 Special Relativistic Hydrodynamics
Conserved variables: D = ρW, S = ρhW^2v, E = ρhW^2 - P
where W = 1/sqrt(1 - v^2) is Lorentz factor, h = 1 + ε + P/ρ is enthalpy

### 12.2 SRMHD
Additional magnetic field evolution with relativistic MHD equations.

### 12.3 SR-GSPH
Godunov SPH with relativistic Riemann solvers.
