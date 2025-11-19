#!/bin/bash

# Iterative Riemann Solver Debug Suite
# This script runs the complete debugging workflow

set -e

echo "=================================================="
echo "Iterative Riemann Solver Debugging Suite"
echo "=================================================="
echo ""

# Change to test directory
cd "$(dirname "$0")"

# Step 1: Compile C++ test program
echo "Step 1: Compiling C++ test program..."
echo "--------------------------------------------------"
g++ -std=c++17 -O3 -Wall -Wextra test_iterative_riemann_solver.cpp -o test_solver
echo "✓ Compilation successful"
echo ""

# Step 2: Run C++ tests
echo "Step 2: Running C++ tests..."
echo "--------------------------------------------------"
./test_solver
echo "✓ C++ tests complete"
echo ""

# Step 3: Run Python tests
echo "Step 3: Running Python tests..."
echo "--------------------------------------------------"
python3 compare_riemann_solvers.py
echo "✓ Python tests complete"
echo ""

# Step 4: Generate visualizations
echo "Step 4: Generating comparison visualizations..."
echo "--------------------------------------------------"
python3 visualize_comparison.py
echo "✓ Visualization complete"
echo ""

# Step 5: Summary
echo "=================================================="
echo "All tests complete!"
echo "=================================================="
echo ""
echo "Generated files:"
echo "  - test_results_cpp.txt     : C++ solver detailed output"
echo "  - test_results_python.txt  : Python solver detailed output"
echo "  - comparison_plots/        : Visualization directory"
echo ""
echo "To view results:"
echo "  - Check comparison_plots/summary_comparison.png for overview"
echo "  - Check individual test plots in comparison_plots/"
echo "  - Review text files for detailed numerical data"
echo ""
