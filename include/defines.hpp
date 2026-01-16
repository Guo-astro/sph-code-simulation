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
constexpr int neighbor_list_size = 5000;

// for debug - DISABLED: tree search validated against exhaustive search
// Uncomment to enable O(N^2) exhaustive search for debugging
// #define EXHAUSTIVE_SEARCH
