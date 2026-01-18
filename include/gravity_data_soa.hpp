#pragma once

#include <vector>
#include <cstddef>
#include "defines.hpp"
#include "vector_type.hpp"
#include "particle.hpp"

/**
 * @file gravity_data_soa.hpp
 * @brief Structure-of-Arrays (SoA) layout for gravity computation
 *
 * This class provides cache-efficient data layout for gravity calculations.
 * Instead of Array-of-Structures (AoS) like SPHParticle[], we use Structure-of-Arrays:
 *
 * AoS (current): particles[i].pos, particles[i].mass (scattered in memory)
 * SoA (new):     pos_x[], pos_y[], pos_z[], mass[]  (contiguous same-type data)
 *
 * Benefits:
 * - Better cache utilization (contiguous access patterns)
 * - Enables SIMD vectorization
 * - Better memory prefetching
 *
 * TDD STUB: These declarations exist only to allow tests to compile.
 * Implementation is pending (TDD red phase).
 */

namespace sph {

// Forward declaration for softening type
enum class GravitySofteningType;

/**
 * @class AlignedVector
 * @brief std::vector with 64-byte (cache line) alignment
 *
 * Wrapper to ensure arrays are cache-line aligned for optimal performance.
 */
template<typename T>
class AlignedVector : public std::vector<T> {
public:
    using std::vector<T>::vector;

    // Ensure data is aligned (platform-dependent)
    T* data() { return std::vector<T>::data(); }
    const T* data() const { return std::vector<T>::data(); }
};

/**
 * @class GravityDataSoA
 * @brief Structure-of-Arrays layout for gravity computation data
 *
 * Provides:
 * - Cache-aligned contiguous arrays for positions, masses, and smoothing lengths
 * - Output arrays for gravity acceleration and potential
 * - Conversion utilities between AoS and SoA layouts
 *
 * TDD STUB: Returns placeholder values that will cause tests to FAIL.
 */
class GravityDataSoA {
public:
    // Input data arrays (cache-aligned)
    AlignedVector<real> pos_x;   // X positions
    AlignedVector<real> pos_y;   // Y positions (DIM >= 2)
    AlignedVector<real> pos_z;   // Z positions (DIM == 3)
    AlignedVector<real> mass;    // Particle masses
    AlignedVector<real> sml;     // Smoothing lengths

    // Output data arrays (cache-aligned)
    AlignedVector<real> grav_acc_x;  // X acceleration
    AlignedVector<real> grav_acc_y;  // Y acceleration (DIM >= 2)
    AlignedVector<real> grav_acc_z;  // Z acceleration (DIM == 3)
    AlignedVector<real> phi;         // Gravitational potential

    /**
     * @brief Construct SoA with given number of particles
     * @param n Number of particles
     */
    explicit GravityDataSoA(size_t n = 0);

    /**
     * @brief Create SoA from AoS particle array
     * @param particles Vector of SPH particles
     * @return GravityDataSoA containing same data in SoA layout
     *
     * TDD STUB: Returns structure with empty arrays
     */
    static GravityDataSoA from_aos(const std::vector<SPHParticle>& particles);

    /**
     * @brief Convert back to AoS layout
     * @return Vector of SPH particles
     *
     * TDD STUB: Returns empty vector
     */
    std::vector<SPHParticle> to_aos() const;

    /**
     * @brief Copy computed gravity results back to AoS particles
     * @param particles Target particle array (modified in place)
     *
     * TDD STUB: Does nothing
     */
    void copy_results_to_aos(std::vector<SPHParticle>& particles) const;

    /**
     * @brief Get number of particles
     * @return Number of particles in data structure
     */
    size_t size() const { return pos_x.size(); }

    /**
     * @brief Resize all arrays
     * @param n New size
     */
    void resize(size_t n);

    /**
     * @brief Clear all arrays
     */
    void clear();

private:
    // Allocate aligned memory for arrays
    void allocate_aligned(size_t n);
};

/**
 * @brief Compute gravity for a single particle using SoA data
 * @param target_idx Index of target particle
 * @param data SoA data containing all particles
 * @param softening Softening length
 * @return Gravitational acceleration vector
 *
 * TDD STUB: Returns zero vector
 */
vec_t compute_gravity_soa_single(int target_idx,
                                  const GravityDataSoA& data,
                                  real softening);

/**
 * @brief Compute gravity at arbitrary position using SoA data
 * @param pos Test position
 * @param data SoA data containing source particles
 * @param softening Softening length
 * @return Gravitational acceleration vector
 *
 * TDD STUB: Returns zero vector
 */
vec_t compute_gravity_soa_at_position(const vec_t& pos,
                                       const GravityDataSoA& data,
                                       real softening);

/**
 * @brief Compute gravity for all particles using SoA layout
 * @param data SoA data (modified in place with grav_acc and phi)
 * @param G Gravitational constant
 * @param type Softening kernel type
 * @param softening Softening length
 *
 * This is the main entry point for SoA-based gravity computation.
 * Results are stored in data.grav_acc_x/y/z and data.phi.
 *
 * TDD STUB: Does nothing
 */
void compute_gravity_soa(GravityDataSoA& data, real G,
                         GravitySofteningType type, real softening);

} // namespace sph
