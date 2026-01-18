#pragma once

/**
 * @file hernquist_katz_lookup_table.hpp
 * @brief Precomputed lookup table for Hernquist-Katz (cubic spline) softening kernel
 *
 * This class provides fast evaluation of the Hernquist & Katz (1989) gravitational
 * potential (f) and force kernel (g) functions using linear interpolation from a
 * precomputed table.
 *
 * The Hernquist-Katz kernel uses:
 * - Softening length: e = h/2
 * - Normalized distance: u = r/e
 * - Support: u in [0, 2] (vs [0, 1] for Wendland C4)
 *
 * TDD STUB: This class declaration exists only to allow tests to compile.
 * Implementation is pending (TDD red phase).
 */

namespace sph {

/**
 * @class HernquistKatzLookupTable
 * @brief Precomputed lookup table for Hernquist-Katz gravity softening
 *
 * STUB: Returns placeholder values that will cause tests to FAIL.
 * This is intentional for TDD red phase - tests written first, implementation later.
 */
class HernquistKatzLookupTable {
public:
    // Table parameters
    // Larger than Wendland C4 (2048) due to 2x support radius [0, 2] vs [0, 1]
    static constexpr int TABLE_SIZE = 4096;
    static constexpr double U_MAX = 2.0;         // Kernel support [0, 2]
    static constexpr double DU = U_MAX / TABLE_SIZE;
    static constexpr double INV_DU = static_cast<double>(TABLE_SIZE) / U_MAX;

    /**
     * @brief Get singleton instance of the lookup table
     * @return Reference to the singleton instance
     *
     * Thread-safe initialization via Meyers' singleton.
     */
    static const HernquistKatzLookupTable& get_instance();

    /**
     * @brief Lookup f(u) - potential kernel (normalized)
     * @param u Normalized distance r/e where e = h/2
     * @return f value (dimensionless, multiply by 1/e for physical value)
     *
     * STUB: Returns 0.0 to cause test failures
     */
    double f(double u) const;

    /**
     * @brief Lookup g(u) - force kernel (normalized)
     * @param u Normalized distance r/e where e = h/2
     * @return g value (dimensionless, multiply by 1/e^3 for physical value)
     *
     * STUB: Returns 0.0 to cause test failures
     */
    double g(double u) const;

    /**
     * @brief Compute f(r, h) with full interface (matching gravity_force.cpp API)
     * @param r Distance between particles
     * @param h Smoothing length
     * @return f(r,h) as defined in gravity_force.cpp
     *
     * STUB: Returns 0.0 to cause test failures
     */
    double f_full(double r, double h) const;

    /**
     * @brief Compute g(r, h) with full interface (matching gravity_force.cpp API)
     * @param r Distance between particles
     * @param h Smoothing length
     * @return g(r,h) as defined in gravity_force.cpp
     *
     * STUB: Returns 0.0 to cause test failures
     */
    double g_full(double r, double h) const;

    /**
     * @brief Get the table size
     * @return Number of entries in lookup table
     */
    int size() const { return TABLE_SIZE; }

private:
    // Private constructor for singleton
    HernquistKatzLookupTable();

    // Prevent copying
    HernquistKatzLookupTable(const HernquistKatzLookupTable&) = delete;
    HernquistKatzLookupTable& operator=(const HernquistKatzLookupTable&) = delete;

    // Cache-aligned lookup tables with precomputed slopes for fast interpolation
    alignas(64) double f_table_[TABLE_SIZE + 1];
    alignas(64) double f_slope_[TABLE_SIZE];
    alignas(64) double g_table_[TABLE_SIZE + 1];
    alignas(64) double g_slope_[TABLE_SIZE];
};

} // namespace sph
