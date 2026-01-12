# SPH Code Simulation - Project Requirements

## Overview
High-performance Smoothed Particle Hydrodynamics (SPH) simulation code for computational fluid dynamics and astrophysical simulations.

## Core Features

### Multi-dimensional Support
- 1D, 2D, and 3D simulations via compile-time dimension selection (`-DSPH_DIM=1/2/3`)
- Dimension-agnostic particle data structures

### SPH Variants
- Standard SPH (SSPH) with various kernels
- Godunov SPH (GSPH) with Riemann solvers
- Gradient-h SPH for variable smoothing length
- Special Relativistic MHD (SRMHD)
- Special Relativistic GSPH (SR-GSPH)

### Physics Modules
- Hydrodynamics with ideal gas EOS
- Self-gravity via Barnes-Hut tree (O(N log N))
- Thermal physics with ISM cooling (Koyama-Inutsuka 2000)
- External forces (point mass gravity, uniform fields)

### Performance
- OpenMP parallelization for multi-core systems
- Barnes-Hut tree for efficient neighbor search and gravity
- HDF5 output for large-scale data

### Configuration
- JSON-based preset system for simulation setup
- Flexible initial condition generators
- Comprehensive logging and diagnostics

## Technical Stack
- Language: C++14
- Build: CMake 3.13+
- Dependencies: OpenMP, Boost, HDF5, nlohmann/json
- Testing: Google Test

## Directory Structure
```
src/           - C++ implementation
include/       - Header files
simulations/   - Preset configs and results
scripts/       - Python analysis tools
tools/sph-viz/ - Web visualization
docs/          - Technical documentation
tests/         - Unit tests
```

## Build Commands
```bash
cd build && cmake .. -DSPH_DIM=2 && make -j8
./sph <config.json>
ctest --output-on-failure
```

## Current Development Areas
- IMBH-cloud tidal disruption simulations
- ISM thermal instability with K&I cooling
- Relativistic MHD validation
- Performance optimization
