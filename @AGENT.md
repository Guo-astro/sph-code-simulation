# Agent Build Instructions - SPH Code Simulation

## Project Setup

### Prerequisites
```bash
# macOS (Homebrew)
brew install cmake boost hdf5 libomp

# Ubuntu/Debian
sudo apt install cmake libboost-all-dev libhdf5-dev
```

### Build Commands
```bash
# Standard 2D build (most common)
cd build && cmake .. -DSPH_DIM=2 && make -j8

# 1D build (for shock tube tests)
cd build && cmake .. -DSPH_DIM=1 && make -j8

# 3D build (for astrophysical simulations)
cd build && cmake .. -DSPH_DIM=3 && make -j8

# Debug build with sanitizers
cmake .. -DCMAKE_BUILD_TYPE=Debug -DENABLE_SANITIZERS=ON && make -j8

# Clean rebuild
rm -rf build/* && cd build && cmake .. && make -j8
```

## Running Tests
```bash
# Run all tests
cd build && ctest --output-on-failure

# Run specific test
./build/test_kernel_gravity
```

## Running Simulations
```bash
# Run with config file
./build/sph <config.json>

# Example workflows (via Make)
make lane_emden_help    # Polytrope simulations
make shock_tube_help    # Shock tube tests
make sedov_help         # Sedov blast wave
make imbh_help          # IMBH-cloud interactions
```

## Visualization
```bash
# Export and view simulation results
make viz SIM=simulations/benchmarks/sedov/results/gsph_wendland

# Export only
make viz_export SIM=<simulation_path>

# Start viz server only
make viz_server
```

## Key Learnings

### Build Options Reference
| Option | Values | Description |
|--------|--------|-------------|
| `SPH_DIM` | 1, 2, 3 | Spatial dimension |
| `CMAKE_BUILD_TYPE` | Release, Debug | Build type |
| `ENABLE_NATIVE_ARCH` | ON/OFF | CPU-specific optimizations |
| `ENABLE_SANITIZERS` | ON/OFF | Address/UB sanitizers |
| `BUILD_TESTS` | ON/OFF | Build unit tests |

### Common Issues
- **Missing OpenMP on macOS**: `brew install libomp`
- **HDF5 not found**: `cmake .. -DHDF5_ROOT=/path/to/hdf5`
- **Compilation errors**: Clean rebuild with `rm -rf build/*`

## Feature Development Quality Standards

### Testing Requirements
- Unit tests use Google Test framework
- Run tests: `cd build && ctest --output-on-failure`
- Validate physics with standard benchmark problems (Sod, Sedov)

### Git Workflow
1. Commit with clear messages: `git commit -m "feat(gsph): add HLLC solver"`
2. Use conventional commits: `feat:`, `fix:`, `docs:`, `test:`, `refactor:`
3. Push regularly to remote

### Documentation Requirements
- Update relevant docs when implementation changes
- Keep @fix_plan.md current with task progress
- Document physics assumptions in code comments

### Feature Completion Checklist
- [ ] All tests pass: `ctest --output-on-failure`
- [ ] Code compiles without warnings
- [ ] Changes committed with descriptive messages
- [ ] @fix_plan.md updated
- [ ] Documentation updated if needed
