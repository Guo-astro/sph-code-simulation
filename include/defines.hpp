#pragma once

#include <cmath>

// DIM must be injected via CMake: cmake -DSPH_DIM=1/2/3 ..
// This ensures dimension is explicitly set per build configuration
#ifndef DIM
#error "DIM not defined! Please configure CMake with: cmake -DSPH_DIM=1 .. (or 2, or 3)"
#endif

typedef double real;

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971693993751
#endif

inline real inner_product(const real (&v1)[DIM], const real (&v2)[DIM])
{
#if DIM == 1
    return v1[0] * v2[0];
#elif DIM == 2
    return v1[0] * v2[0] + v1[1] * v2[1];
#elif DIM == 3
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
#endif
}

inline real sqr(real x) { return x * x; }
inline real pow3(real x) { return x * x * x; }
inline real pow4(real x) { return x * x * x * x; }
inline real pow5(real x) { return x * x * x * x * x; }
inline real pow6(real x) { return x * x * x * x * x * x; }

// Lane-Emden has steep density gradient (~10^7:1), needs large neighbor list
// SR Sod also has sharp discontinuity requiring large neighbor lists
// For 1M+ particles, need larger buffer during initial smoothing iterations
constexpr int neighbor_list_size = 20000;


// =============================================================================
// Tree Optimization Switches
// =============================================================================

// Morton code particle reordering for cache-friendly tree traversal
// Reorders particles along Z-order space-filling curve to improve cache locality
// Expected speedup: 2-3x due to reduced cache misses
#define USE_MORTON_ORDERING

// Iterative tree traversal instead of recursive
// Uses explicit stack to eliminate function call overhead
// Expected speedup: 10-20%
#define USE_ITERATIVE_TRAVERSAL

// =============================================================================
// Performance Optimization Switches
// =============================================================================

// Thread-local neighbor list storage to avoid heap allocation in hot loops
// Moves neighbor_list outside the loop and uses thread-local storage
// Expected speedup: 5-15% by reducing allocator contention
#define SPH_USE_THREAD_LOCAL_NEIGHBOR_LIST

// Distance caching before sorting neighbors
// Caches distances to avoid recomputing during sort comparison
// Expected speedup: 10-30% for neighbor-heavy operations
#define SPH_USE_DISTANCE_CACHING

// OpenMP dynamic scheduling for load balancing
// Uses schedule(dynamic, 64) for better load distribution with variable density
// Expected speedup: 5-20% on non-uniform particle distributions
#define SPH_USE_DYNAMIC_SCHEDULING

// Template-based kernel dispatch to eliminate virtual function overhead
// Uses CRTP pattern for compile-time polymorphism instead of virtual calls
// Expected speedup: 10-15% in kernel-heavy loops
#define SPH_USE_TEMPLATE_KERNEL

