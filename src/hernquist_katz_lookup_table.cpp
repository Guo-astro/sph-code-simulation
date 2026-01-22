/**
 * @file hernquist_katz_lookup_table_stub.cpp
 * @brief Implementation of Hernquist-Katz lookup table
 *
 * Precomputed lookup table for Hernquist & Katz (1989) gravitational
 * potential (f) and force kernel (g) functions using linear interpolation.
 */

#include "hernquist_katz_lookup_table.hpp"
#include <cmath>

namespace sph {

// Static member definitions for constexpr (C++14 requirement)
constexpr int HernquistKatzLookupTable::TABLE_SIZE;
constexpr double HernquistKatzLookupTable::U_MAX;
constexpr double HernquistKatzLookupTable::DU;
constexpr double HernquistKatzLookupTable::INV_DU;

/**
 * @brief Compute f(u) polynomial - potential kernel (normalized)
 * @param u Normalized distance r/e where e = h/2
 * @return Dimensionless potential value
 */
static double compute_f_polynomial(double u) {
    if (u >= 2.0) {
        return 2.0 / u;
    }

    if (u < 1.0) {
        // f = (-0.5*u^2*(1/3 - 3/20*u^2 + u^3/20) + 1.4) / e
        // Normalized: f * e = (-0.5*u^2*(1/3 - 3/20*u^2 + u^3/20) + 1.4)
        double u2 = u * u;
        double u3 = u2 * u;
        return -0.5 * u2 * (1.0/3.0 - 3.0/20.0 * u2 + u3/20.0) + 1.4;
    } else {
        // 1 <= u < 2
        double u2 = u * u;
        double u3 = u2 * u;
        return -1.0/(15.0*u) + (-u2 * (4.0/3.0 - u + 0.3*u2 - u3/30.0) + 1.6);
    }
}

/**
 * @brief Compute g(u) polynomial - force kernel (normalized)
 * @param u Normalized distance r/e where e = h/2
 * @return Dimensionless force value
 */
static double compute_g_polynomial(double u) {
    if (u >= 2.0) {
        return 1.0 / (u * u * u);
    }

    if (u < 1e-10) {
        // At center: g = 4/3
        return 4.0 / 3.0;
    }

    if (u < 1.0) {
        // g = (4/3 - 1.2*u^2 + 0.5*u^3) / e^3
        // Normalized: g * e^3 = 4/3 - 1.2*u^2 + 0.5*u^3
        double u2 = u * u;
        double u3 = u2 * u;
        return 4.0/3.0 - 1.2*u2 + 0.5*u3;
    } else {
        // 1 <= u < 2
        double u2 = u * u;
        double u3 = u2 * u;
        double u4 = u2 * u2;
        double u5 = u4 * u;
        double u6 = u3 * u3;
        return (-1.0/15.0 + 8.0/3.0*u3 - 3.0*u4 + 1.2*u5 - u6/6.0) / u3;
    }
}

HernquistKatzLookupTable::HernquistKatzLookupTable() {
    // Initialize tables with polynomial values at grid points
    for (int i = 0; i <= TABLE_SIZE; ++i) {
        double u = static_cast<double>(i) * DU;
        f_table_[i] = compute_f_polynomial(u);
        g_table_[i] = compute_g_polynomial(u);
    }

    // Precompute slopes for fast linear interpolation
    for (int i = 0; i < TABLE_SIZE; ++i) {
        f_slope_[i] = (f_table_[i + 1] - f_table_[i]) * INV_DU;
        g_slope_[i] = (g_table_[i + 1] - g_table_[i]) * INV_DU;
    }
}

const HernquistKatzLookupTable& HernquistKatzLookupTable::get_instance() {
    // Thread-safe initialization via Meyers' singleton (C++11 static local)
    static HernquistKatzLookupTable instance;
    return instance;
}

double HernquistKatzLookupTable::f(double u) const {
    // Clamp u to valid range [0, 2]
    if (u <= 0.0) {
        return f_table_[0];
    }
    if (u >= U_MAX) {
        return 2.0 / u;  // Point mass: 1/r with normalization
    }

    // Check if we're near a discontinuity boundary
    // Boundaries are at u=1 (index 2048) and u=2 (index 4096)
    // Use direct polynomial evaluation near boundaries to avoid interpolation artifacts
    constexpr int INNER_BOUNDARY_IDX = TABLE_SIZE / 2;  // u = 1.0
    constexpr int OUTER_BOUNDARY_IDX = TABLE_SIZE - 1;  // last index before u = 2.0

    double idx_f = u * INV_DU;
    int idx = static_cast<int>(idx_f);

    // If interpolating across u=1 or u=2 boundary, use polynomial directly
    if (idx == INNER_BOUNDARY_IDX - 1 || idx == INNER_BOUNDARY_IDX ||
        idx >= OUTER_BOUNDARY_IDX) {
        return compute_f_polynomial(u);
    }

    double frac = idx_f - idx;

    return f_table_[idx] + frac * (f_table_[idx + 1] - f_table_[idx]);
}

double HernquistKatzLookupTable::g(double u) const {
    // Clamp u to valid range [0, 2]
    if (u <= 0.0) {
        return g_table_[0];
    }
    if (u >= U_MAX) {
        return 1.0 / (u * u * u);  // Point mass: 1/r^3
    }

    // Check if we're near a discontinuity boundary
    // Boundaries are at u=1 (index 2048) and u=2 (index 4096)
    constexpr int INNER_BOUNDARY_IDX = TABLE_SIZE / 2;  // u = 1.0
    constexpr int OUTER_BOUNDARY_IDX = TABLE_SIZE - 1;  // last index before u = 2.0

    double idx_f = u * INV_DU;
    int idx = static_cast<int>(idx_f);

    // If interpolating across u=1 or u=2 boundary, use polynomial directly
    if (idx == INNER_BOUNDARY_IDX - 1 || idx == INNER_BOUNDARY_IDX ||
        idx >= OUTER_BOUNDARY_IDX) {
        return compute_g_polynomial(u);
    }

    double frac = idx_f - idx;

    return g_table_[idx] + frac * (g_table_[idx + 1] - g_table_[idx]);
}

double HernquistKatzLookupTable::f_full(double r, double h) const {
    // e = h/2, u = r/e = 2r/h
    double e = h * 0.5;
    double u = r / e;

    if (u >= 2.0) {
        // Point mass regime: 1/r
        return (r > 1e-30) ? 1.0 / r : 0.0;
    }

    // Inside kernel: f(u) / e
    return f(u) / e;
}

double HernquistKatzLookupTable::g_full(double r, double h) const {
    // e = h/2, u = r/e = 2r/h
    double e = h * 0.5;
    double u = r / e;

    if (u >= 2.0) {
        // Point mass regime: 1/r^3
        return (r > 1e-30) ? 1.0 / (r * r * r) : 0.0;
    }

    // Inside kernel: g(u) / e^3
    double e3 = e * e * e;
    return g(u) / e3;
}

} // namespace sph
