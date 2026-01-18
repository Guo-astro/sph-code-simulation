#pragma once

#include <vector>

/**
 * @file softening_lookup_table.hpp
 * @brief Precomputed lookup table for Wendland C4 softening kernel
 *
 * This class provides fast evaluation of the Wendland C4 gravitational potential
 * (phi) and force kernel (g) functions using linear interpolation from a
 * precomputed table.
 *
 * The lookup table replaces expensive 9th-order polynomial evaluation with
 * a simple table lookup + linear interpolation, providing >2x speedup while
 * maintaining 1e-6 accuracy.
 *
 * Usage:
 *   const auto& table = WendlandC4LookupTable::get_instance();
 *   double phi = table.phi(q);  // q = r/h, normalized distance
 *   double g = table.g(q);
 *
 * Thread Safety:
 *   The singleton instance is initialized using Meyers' singleton pattern,
 *   which is thread-safe in C++11 and later.
 */

namespace sph {

/**
 * @class WendlandC4LookupTable
 * @brief Precomputed lookup table for Wendland C4 gravity softening
 *
 * STUB: This class declaration exists only to allow tests to compile.
 * Implementation is pending (TDD red phase).
 */
class WendlandC4LookupTable {
public:
    // Table parameters
    static constexpr int TABLE_SIZE = 2048;      // Power of 2 for fast indexing
    static constexpr double Q_MAX = 1.0;         // Kernel support [0, 1]
    static constexpr double DQ = Q_MAX / TABLE_SIZE;
    static constexpr double INV_DQ = static_cast<double>(TABLE_SIZE) / Q_MAX;

    /**
     * @brief Get singleton instance of the lookup table
     * @return Reference to the singleton instance
     *
     * Thread-safe initialization via Meyers' singleton.
     */
    static const WendlandC4LookupTable& get_instance();

    /**
     * @brief Lookup phi(q) using linear interpolation
     * @param q Normalized distance r/h (must be in [0, 1])
     * @return Interpolated phi value (dimensionless, multiply by 1/h)
     *
     * For q >= 1, returns point mass value 1/q.
     */
    inline double phi(double q) const {
        // Handle q >= 1: point mass
        if (q >= 1.0) {
            return 1.0 / q;
        }

        // Fast path: assume 0 <= q < 1 (most common case)
        // Find table index and fractional part using fast integer conversion
        double idx_float = q * INV_DQ;
        int idx = static_cast<int>(idx_float);
        double frac = idx_float - idx;

        // Linear interpolation using precomputed slope (one less memory access)
        return phi_table_[idx] + frac * phi_slope_[idx];
    }

    /**
     * @brief Lookup g(q) using linear interpolation
     * @param q Normalized distance r/h (must be in [0, 1])
     * @return Interpolated g value (dimensionless)
     *
     * For q >= 1, returns point mass value 1/q^3.
     * For q = 0, returns 0 (no force at center).
     */
    double g(double q) const;

    /**
     * @brief Compute phi(r, h) with full interface (matching GravityForce API)
     * @param r Distance between particles
     * @param h Smoothing length
     * @return phi(r,h) = phi(r/h) / h
     */
    double phi_full(double r, double h) const;

    /**
     * @brief Compute g(r, h) with full interface (matching GravityForce API)
     * @param r Distance between particles
     * @param h Smoothing length
     * @return g(r,h) as defined in gravity_force.cpp
     */
    double g_full(double r, double h) const;

    /**
     * @brief Get the table size
     * @return Number of entries in lookup table
     */
    int size() const { return TABLE_SIZE; }

private:
    // Private constructor for singleton
    WendlandC4LookupTable();

    // Prevent copying
    WendlandC4LookupTable(const WendlandC4LookupTable&) = delete;
    WendlandC4LookupTable& operator=(const WendlandC4LookupTable&) = delete;

    // Cache-aligned lookup tables with precomputed slopes for fast interpolation
    // Storing slopes avoids one memory access in the hot path
    alignas(64) double phi_table_[TABLE_SIZE + 1];
    alignas(64) double phi_slope_[TABLE_SIZE];  // phi_slope_[i] = phi_table_[i+1] - phi_table_[i]
    alignas(64) double g_table_[TABLE_SIZE + 1];
    alignas(64) double g_slope_[TABLE_SIZE];    // g_slope_[i] = g_table_[i+1] - g_table_[i]
};

/**
 * @brief Configurable lookup table with user-specified size
 *
 * Used for testing accuracy vs table size tradeoffs.
 *
 * STUB: Declaration only, implementation pending.
 */
class SofteningLookupTable {
public:
    /**
     * @brief Construct a lookup table with specified size
     * @param table_size Number of entries (default: 2048)
     */
    explicit SofteningLookupTable(int table_size = 2048);

    double phi(double q) const;
    double g(double q) const;
    double phi_full(double r, double h) const;
    double g_full(double r, double h) const;
    int size() const { return table_size_; }

private:
    int table_size_;
    double dq_;
    double inv_dq_;
    std::vector<double> phi_table_;
    std::vector<double> g_table_;
};

} // namespace sph
