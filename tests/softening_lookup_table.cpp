/**
 * @file softening_lookup_table_stub.cpp
 * @brief Implementation of Wendland C4 softening lookup table
 *
 * Provides fast evaluation of gravitational potential (phi) and force kernel (g)
 * using linear interpolation from precomputed tables.
 */

#include "softening_lookup_table.hpp"
#include <cmath>

namespace sph {

// =============================================================================
// Polynomial coefficients from gravity_force.cpp
// =============================================================================

namespace {

// Potential coefficients: phi(q) = a0 + a1*q + a2*q^2 + ... + a9*q^9
constexpr double a0 =  3.4374743761;
constexpr double a1 = -0.0031873250;
constexpr double a2 = -10.2154807743;
constexpr double a3 = -1.1577720555;
constexpr double a4 =  36.1013669755;
constexpr double a5 = -26.3399094060;
constexpr double a6 = -44.1079372114;
constexpr double a7 =  82.6543766683;
constexpr double a8 = -50.5921624056;
constexpr double a9 =  11.2232565249;

// Derivative coefficients: dphi/dq = b0 + b1*q + ... + b8*q^8
constexpr double b0 = -0.0031873250;
constexpr double b1 = -20.4309615486;
constexpr double b2 = -3.4733161665;
constexpr double b3 = 144.4054679020;
constexpr double b4 = -131.6995470300;
constexpr double b5 = -264.6476232684;
constexpr double b6 = 578.5806366781;
constexpr double b7 = -404.7372992448;
constexpr double b8 = 101.0093087241;

/**
 * @brief Evaluate phi polynomial at q (dimensionless, without 1/h factor)
 */
double eval_phi_poly(double q) {
    const double q2 = q * q;
    const double q3 = q2 * q;
    const double q4 = q2 * q2;
    const double q5 = q4 * q;
    const double q6 = q3 * q3;
    const double q7 = q6 * q;
    const double q8 = q4 * q4;
    const double q9 = q8 * q;

    return a0 + a1*q + a2*q2 + a3*q3 + a4*q4 + a5*q5 + a6*q6 + a7*q7 + a8*q8 + a9*q9;
}

/**
 * @brief Evaluate g(q) = -dphi_dq / q (dimensionless, without 1/h^3 factor)
 */
double eval_g_poly(double q) {
    if (q < 1e-10) {
        return 0.0;  // Force is zero at center
    }

    const double q2 = q * q;
    const double q3 = q2 * q;
    const double q4 = q2 * q2;
    const double q5 = q4 * q;
    const double q6 = q3 * q3;
    const double q7 = q6 * q;

    const double dphi_dq = b0 + b1*q + b2*q2 + b3*q3 + b4*q4 + b5*q5 + b6*q6 + b7*q7 + b8*q*q7;

    return -dphi_dq / q;
}

} // anonymous namespace

// =============================================================================
// WendlandC4LookupTable Implementation
// =============================================================================

// Static member definitions for constexpr (C++14 requirement)
constexpr int WendlandC4LookupTable::TABLE_SIZE;
constexpr double WendlandC4LookupTable::Q_MAX;
constexpr double WendlandC4LookupTable::DQ;
constexpr double WendlandC4LookupTable::INV_DQ;

WendlandC4LookupTable::WendlandC4LookupTable() {
    // Initialize tables with polynomial values at each grid point
    for (int i = 0; i <= TABLE_SIZE; ++i) {
        double q = static_cast<double>(i) * DQ;
        phi_table_[i] = eval_phi_poly(q);
        g_table_[i] = eval_g_poly(q);
    }

    // Precompute slopes for fast linear interpolation
    for (int i = 0; i < TABLE_SIZE; ++i) {
        phi_slope_[i] = phi_table_[i + 1] - phi_table_[i];
        g_slope_[i] = g_table_[i + 1] - g_table_[i];
    }
}

const WendlandC4LookupTable& WendlandC4LookupTable::get_instance() {
    // Thread-safe initialization via Meyers' singleton (C++11 static local)
    static WendlandC4LookupTable instance;
    return instance;
}

// Note: phi() is defined inline in the header for performance

double WendlandC4LookupTable::g(double q) const {
    // Handle q >= 1: point mass
    if (q >= 1.0) {
        return 1.0 / (q * q * q);
    }

    // Handle q = 0: force is zero at center
    if (q <= 0.0) {
        return 0.0;
    }

    // For small q, linear interpolation has high error due to steep gradient.
    // Use polynomial evaluation directly for q < threshold where threshold
    // is chosen so linear interpolation error < 1e-6.
    // With TABLE_SIZE=2048, DQ=1/2048, threshold of 40*DQ covers q~0.02.
    constexpr double SMALL_Q_THRESHOLD = 40.0 * DQ;  // ~0.02
    if (q < SMALL_Q_THRESHOLD) {
        return eval_g_poly(q);
    }

    // Find table index and fractional part
    double idx_float = q * INV_DQ;
    int idx = static_cast<int>(idx_float);
    double frac = idx_float - idx;

    // Clamp index to valid range
    if (idx >= TABLE_SIZE) {
        idx = TABLE_SIZE - 1;
        frac = 1.0;
    }

    // Linear interpolation
    return g_table_[idx] + frac * (g_table_[idx + 1] - g_table_[idx]);
}

double WendlandC4LookupTable::phi_full(double r, double h) const {
    double q = r / h;

    if (q >= 1.0) {
        // Outside kernel: point mass potential phi = 1/r
        return 1.0 / r;
    }

    // Inside kernel: phi(r,h) = phi(q) / h
    return phi(q) / h;
}

double WendlandC4LookupTable::g_full(double r, double h) const {
    double q = r / h;

    if (q >= 1.0) {
        // Outside kernel: point mass g = 1/r^3
        return 1.0 / (r * r * r);
    }

    // Handle r = 0
    if (r < 1e-30) {
        return 0.0;
    }

    // Inside kernel: g(r,h) = -dphi_dq / (h^3 * q)
    // The lookup returns g(q) = -dphi_dq / q (dimensionless)
    // So we need to divide by h^3
    double h3 = h * h * h;
    return g(q) / h3;
}

// =============================================================================
// SofteningLookupTable Implementation (configurable size)
// =============================================================================

SofteningLookupTable::SofteningLookupTable(int table_size)
    : table_size_(table_size)
    , dq_(WendlandC4LookupTable::Q_MAX / table_size)
    , inv_dq_(static_cast<double>(table_size) / WendlandC4LookupTable::Q_MAX)
    , phi_table_(table_size + 1)
    , g_table_(table_size + 1)
{
    // Initialize tables with polynomial values
    for (int i = 0; i <= table_size_; ++i) {
        double q = static_cast<double>(i) * dq_;
        phi_table_[i] = eval_phi_poly(q);
        g_table_[i] = eval_g_poly(q);
    }
}

double SofteningLookupTable::phi(double q) const {
    // Handle q >= 1: point mass
    if (q >= 1.0) {
        return 1.0 / q;
    }

    // Clamp negative q to zero
    if (q < 0.0) {
        q = 0.0;
    }

    // Find table index and fractional part
    double idx_float = q * inv_dq_;
    int idx = static_cast<int>(idx_float);
    double frac = idx_float - idx;

    // Clamp index to valid range
    if (idx >= table_size_) {
        idx = table_size_ - 1;
        frac = 1.0;
    }

    // Linear interpolation
    return phi_table_[idx] + frac * (phi_table_[idx + 1] - phi_table_[idx]);
}

double SofteningLookupTable::g(double q) const {
    // Handle q >= 1: point mass
    if (q >= 1.0) {
        return 1.0 / (q * q * q);
    }

    // Handle q = 0: force is zero at center
    if (q <= 0.0) {
        return 0.0;
    }

    // For small q, linear interpolation has high error due to steep gradient.
    // Use polynomial evaluation directly for q < threshold.
    double small_q_threshold = 40.0 * dq_;
    if (q < small_q_threshold) {
        return eval_g_poly(q);
    }

    // Find table index and fractional part
    double idx_float = q * inv_dq_;
    int idx = static_cast<int>(idx_float);
    double frac = idx_float - idx;

    // Clamp index to valid range
    if (idx >= table_size_) {
        idx = table_size_ - 1;
        frac = 1.0;
    }

    // Linear interpolation
    return g_table_[idx] + frac * (g_table_[idx + 1] - g_table_[idx]);
}

double SofteningLookupTable::phi_full(double r, double h) const {
    double q = r / h;

    if (q >= 1.0) {
        return 1.0 / r;
    }

    return phi(q) / h;
}

double SofteningLookupTable::g_full(double r, double h) const {
    double q = r / h;

    if (q >= 1.0) {
        return 1.0 / (r * r * r);
    }

    if (r < 1e-30) {
        return 0.0;
    }

    double h3 = h * h * h;
    return g(q) / h3;
}

} // namespace sph
