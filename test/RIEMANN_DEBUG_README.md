# Iterative Riemann Solver Debugging Suite

This directory contains tools to debug and compare the iterative relativistic Riemann solver implementation in `sr_fluid_force.cpp` with the reference Python implementation in `kitajima_solver.py`.

## Overview

The debugging suite consists of:

1. **C++ Test Program** (`test_iterative_riemann_solver.cpp`)
   - Standalone program that extracts the iterative solver from `sr_fluid_force.cpp`
   - Runs predefined test cases with detailed debug output
   - Outputs iteration history and convergence information

2. **Python Comparison Script** (`compare_riemann_solvers.py`)
   - Runs the same test cases using `kitajima_solver.py`
   - Compares results with C++ implementation
   - Generates detailed comparison reports

3. **Visualization Script** (`visualize_comparison.py`)
   - Parses output from both solvers
   - Creates comparison plots showing:
     - Convergence history
     - Iteration-by-iteration pressure evolution
     - Residual reduction
     - Final results comparison
     - Error metrics

4. **Automated Run Script** (`run_riemann_debug.sh`)
   - Runs the complete debugging workflow
   - Compiles C++, runs both solvers, generates visualizations

## Test Cases

The suite includes 5 standard test cases:

1. **SR_Sod_Shock_Tube**: Standard relativistic Sod shock tube
   - Left: P=13.33, n=10.0, v=0.0
   - Right: P=1e-8, n=1.0, v=0.0

2. **Two_Rarefactions**: Symmetric rarefaction waves
   - Left: P=1.0, n=1.0, v=-0.5
   - Right: P=1.0, n=1.0, v=0.5

3. **Two_Shocks**: Symmetric shock waves
   - Left: P=1.0, n=1.0, v=0.5
   - Right: P=1.0, n=1.0, v=-0.5

4. **Left_Shock_Right_Rarefaction**: Mixed wave pattern
   - Left: P=10.0, n=1.0, v=0.3
   - Right: P=1.0, n=1.0, v=0.0

5. **With_Tangential_Velocity**: Test tangential velocity coupling (Pons et al. 1999)
   - Left: P=1.0, n=1.0, v=0.0, v_t=0.3
   - Right: P=0.1, n=1.0, v=0.0, v_t=0.0

## Usage

### Quick Start

Run the complete debugging suite:

```bash
cd test
./run_riemann_debug.sh
```

This will:
1. Compile the C++ test program
2. Run C++ tests and save results to `test_results_cpp.txt`
3. Run Python tests and save results to `test_results_python.txt`
4. Generate comparison plots in `comparison_plots/`

### Manual Steps

If you prefer to run steps individually:

```bash
# 1. Compile C++ test
cd test
g++ -std=c++17 -O3 -Wall -Wextra test_iterative_riemann_solver.cpp -o test_solver

# 2. Run C++ test
./test_solver

# 3. Run Python comparison
python3 compare_riemann_solvers.py

# 4. Generate visualizations
python3 visualize_comparison.py
```

## Output Files

### Text Output

- **test_results_cpp.txt**: Detailed C++ solver output
  - Input states for each test
  - Iteration-by-iteration convergence data
  - Final pressure P* and velocity v*
  - Convergence status and residuals

- **test_results_python.txt**: Python solver output
  - Same test cases run through kitajima_solver.py
  - Final P* and v* values for comparison

### Visualization Output

All plots are saved in `comparison_plots/`:

- **summary_comparison.png**: Overview of all tests
  - Pressure P* relative errors across all tests
  - Velocity v* absolute errors across all tests

- **{TestName}_comparison.png**: Detailed plot for each test
  - Convergence history (pressure)
  - Convergence history (residual)
  - Initial conditions table
  - Results comparison table
  - Error metrics bar chart

## Interpreting Results

### Expected Behavior

Both solvers should produce nearly identical results:
- **Pressure P***: Relative error < 1e-10 (typically ~1e-14)
- **Velocity v***: Absolute error < 1e-10
- **Convergence**: C++ should converge in 5-15 iterations for most cases

### Common Issues

1. **Large P* errors**: Check initial guess calculation
2. **Non-convergence**: Check derivative computation or relaxation factor
3. **Velocity mismatch**: Check post-wave velocity calculation
4. **Tangential velocity issues**: Check Lorentz factor computation

### Debugging Workflow

If discrepancies are found:

1. Check `test_results_cpp.txt` for iteration history
2. Look at convergence plots to see where divergence occurs
3. Compare intermediate values (H, cs, gamma) between solvers
4. Add debug prints to specific iteration in C++ code
5. Re-run with modified test case to isolate issue

## Adding New Test Cases

To add a new test case:

1. **In C++ file** (`test_iterative_riemann_solver.cpp`):
   ```cpp
   test_cases.push_back({
       "Your_Test_Name",
       5.0/3.0,  // gamma
       1.0,      // c
       {vl, nl, Pl, 0.0, vtl},  // left state
       {vr, nr, Pr, 0.0, vtr}   // right state
   });
   ```

2. **In Python script** (`compare_riemann_solvers.py`):
   ```python
   tc = TestCase(
       "Your_Test_Name",
       5.0/3.0,
       1.0,
       [vl, nl, Pl, 0.0, vtl],
       [vr, nr, Pr, 0.0, vtr]
   )
   tc.left[3] = compute_sound_speed(tc.left[2], tc.left[1], tc.gamma, tc.c)
   tc.right[3] = compute_sound_speed(tc.right[2], tc.right[1], tc.gamma, tc.c)
   test_cases.append(tc)
   ```

3. Re-run the suite

## Dependencies

### C++ Requirements
- C++17 compatible compiler (g++, clang++)
- Standard library (no external dependencies)

### Python Requirements
- Python 3.7+
- numpy
- matplotlib
- kitajima_solver.py from relativistic-riemann package

Install Python dependencies:
```bash
pip install numpy matplotlib
```

## References

The iterative solver is based on:

1. **Pons, Martí & Müller (1999)**
   "The exact solution of the Riemann problem with non-zero tangential velocities in relativistic hydrodynamics"
   J. Fluid Mech., vol. 422, pp. 125-139

2. **Kitajima, Inutsuka, Seno (2025)**
   "Special Relativistic Hydrodynamics with Smoothed Particle Hydrodynamics"
   arXiv:2510.18251v1

## Troubleshooting

### Compilation Errors

If you get compilation errors:
- Check that you have C++17 support: `g++ --version`
- Try with clang++: `clang++ -std=c++17 -O3 test_iterative_riemann_solver.cpp -o test_solver`

### Python Import Errors

If Python can't find kitajima_solver:
- Make sure you're running from the test/ directory
- Check that relativistic-riemann/ exists in parent directory
- Verify kitajima_solver.py is present

### Visualization Errors

If matplotlib fails:
- Install: `pip install matplotlib`
- For headless systems, add: `export MPLBACKEND=Agg` before running

## Contact

For issues or questions about this debugging suite, please refer to the main repository documentation or open an issue.
