#!/bin/bash
# Compare performance with and without optimizations

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

echo "=== Performance Comparison Benchmark ==="
echo ""

# Build directory for optimized version
OPTIMIZED_DIR="$PROJECT_DIR/build_optimized"
mkdir -p "$OPTIMIZED_DIR"

# Build directory for baseline version
BASELINE_DIR="$PROJECT_DIR/build_baseline"
mkdir -p "$BASELINE_DIR"

echo "Building OPTIMIZED version..."
cd "$OPTIMIZED_DIR"
cmake "$PROJECT_DIR" -DSPH_DIM=3 -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-DUSE_MORTON_ORDERING -DUSE_ITERATIVE_TRAVERSAL -DSPH_USE_THREAD_LOCAL_NEIGHBOR_LIST -DSPH_USE_DISTANCE_CACHING -DSPH_USE_DYNAMIC_SCHEDULING" \
    > /dev/null 2>&1
make benchmark_performance -j8 > /dev/null 2>&1

echo "Building BASELINE version (no optimizations)..."
cd "$BASELINE_DIR"
cmake "$PROJECT_DIR" -DSPH_DIM=3 -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-UUSE_MORTON_ORDERING -UUSE_ITERATIVE_TRAVERSAL -USPH_USE_THREAD_LOCAL_NEIGHBOR_LIST -USPH_USE_DISTANCE_CACHING -USPH_USE_DYNAMIC_SCHEDULING" \
    > /dev/null 2>&1

# Remove optimization defines from defines.hpp for baseline build
sed -i.bak 's/#define USE_MORTON_ORDERING/\/\/ #define USE_MORTON_ORDERING/' "$PROJECT_DIR/include/defines.hpp"
sed -i.bak 's/#define USE_ITERATIVE_TRAVERSAL/\/\/ #define USE_ITERATIVE_TRAVERSAL/' "$PROJECT_DIR/include/defines.hpp"
sed -i.bak 's/#define SPH_USE_THREAD_LOCAL_NEIGHBOR_LIST/\/\/ #define SPH_USE_THREAD_LOCAL_NEIGHBOR_LIST/' "$PROJECT_DIR/include/defines.hpp"
sed -i.bak 's/#define SPH_USE_DISTANCE_CACHING/\/\/ #define SPH_USE_DISTANCE_CACHING/' "$PROJECT_DIR/include/defines.hpp"
sed -i.bak 's/#define SPH_USE_DYNAMIC_SCHEDULING/\/\/ #define SPH_USE_DYNAMIC_SCHEDULING/' "$PROJECT_DIR/include/defines.hpp"

make benchmark_performance -j8 > /dev/null 2>&1

# Restore optimization defines
mv "$PROJECT_DIR/include/defines.hpp.bak" "$PROJECT_DIR/include/defines.hpp"

echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║                    BASELINE RESULTS                          ║"
echo "╚══════════════════════════════════════════════════════════════╝"
cd "$BASELINE_DIR"
./benchmark_performance --gtest_filter=PerformanceBenchmark.TreeConstruction 2>&1 | grep -A 10 "=== Tree"
./benchmark_performance --gtest_filter=PerformanceBenchmark.FullNeighborSearchPass 2>&1 | grep -A 10 "=== Full"

echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║                   OPTIMIZED RESULTS                          ║"
echo "╚══════════════════════════════════════════════════════════════╝"
cd "$OPTIMIZED_DIR"
./benchmark_performance --gtest_filter=PerformanceBenchmark.TreeConstruction 2>&1 | grep -A 10 "=== Tree"
./benchmark_performance --gtest_filter=PerformanceBenchmark.FullNeighborSearchPass 2>&1 | grep -A 10 "=== Full"

echo ""
echo "Benchmark comparison complete!"
