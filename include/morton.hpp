#pragma once

#include <cstdint>
#include <vector>
#include <algorithm>
#include <numeric>

#include "defines.hpp"
#include "particle.hpp"

namespace sph {
namespace morton {

// Number of bits per dimension for Morton code
// 21 bits per dimension for 3D (63 bits total in uint64_t)
// 32 bits per dimension for 2D (64 bits total)
// 64 bits for 1D
#if DIM == 3
constexpr int MORTON_BITS_PER_DIM = 21;
#elif DIM == 2
constexpr int MORTON_BITS_PER_DIM = 32;
#else
constexpr int MORTON_BITS_PER_DIM = 64;
#endif

// Bit interleaving helper functions for Morton code encoding

#if DIM == 3
// Spread 21 bits into 63 bits (every 3rd bit)
// Input:  00000000 00000000 00000XXX XXXXXXXX XXXXXXXX (21 bits)
// Output: 00X00X00 X00X00X0 0X00X00X 00X00X00 X00X00X0 0X00X00X 00X00X00 X (63 bits spread)
inline uint64_t spread_bits_3d(uint64_t x) {
    // Mask to keep only 21 bits
    x &= 0x1FFFFF;

    // Magic bit-spreading using parallel bit deposit pattern
    // This spreads 21 bits across 63 positions (every 3rd bit)
    x = (x | (x << 32)) & 0x1F00000000FFFFULL;
    x = (x | (x << 16)) & 0x1F0000FF0000FFULL;
    x = (x | (x << 8))  & 0x100F00F00F00F00FULL;
    x = (x | (x << 4))  & 0x10C30C30C30C30C3ULL;
    x = (x | (x << 2))  & 0x1249249249249249ULL;

    return x;
}

// Encode 3D position into Morton code (Z-order curve)
inline uint64_t encode_3d(uint64_t x, uint64_t y, uint64_t z) {
    return spread_bits_3d(x) | (spread_bits_3d(y) << 1) | (spread_bits_3d(z) << 2);
}
#endif

#if DIM == 2
// Spread 32 bits into 64 bits (every 2nd bit)
// Input:  XXXXXXXX XXXXXXXX XXXXXXXX XXXXXXXX (32 bits)
// Output: 0X0X0X0X 0X0X0X0X 0X0X0X0X 0X0X0X0X 0X0X0X0X 0X0X0X0X 0X0X0X0X 0X0X0X0X (64 bits)
inline uint64_t spread_bits_2d(uint64_t x) {
    x &= 0xFFFFFFFF;
    x = (x | (x << 16)) & 0x0000FFFF0000FFFFULL;
    x = (x | (x << 8))  & 0x00FF00FF00FF00FFULL;
    x = (x | (x << 4))  & 0x0F0F0F0F0F0F0F0FULL;
    x = (x | (x << 2))  & 0x3333333333333333ULL;
    x = (x | (x << 1))  & 0x5555555555555555ULL;
    return x;
}

// Encode 2D position into Morton code
inline uint64_t encode_2d(uint64_t x, uint64_t y) {
    return spread_bits_2d(x) | (spread_bits_2d(y) << 1);
}
#endif

// Compute Morton code from normalized position [0, 1]^DIM
// Position is normalized to integer range [0, 2^MORTON_BITS_PER_DIM - 1]
inline uint64_t encode(const real pos[DIM],
                       const real domain_min[DIM],
                       const real domain_size[DIM]) {
    // Normalize position to [0, 1] range, then scale to integer
    const uint64_t max_coord = (1ULL << MORTON_BITS_PER_DIM) - 1;

#if DIM == 1
    real normalized = (pos[0] - domain_min[0]) / domain_size[0];
    // Clamp to [0, 1] to handle boundary particles
    normalized = std::max(0.0, std::min(1.0, normalized));
    return static_cast<uint64_t>(normalized * max_coord);

#elif DIM == 2
    uint64_t ix[2];
    for (int d = 0; d < 2; ++d) {
        real normalized = (pos[d] - domain_min[d]) / domain_size[d];
        normalized = std::max(0.0, std::min(1.0, normalized));
        ix[d] = static_cast<uint64_t>(normalized * max_coord);
    }
    return encode_2d(ix[0], ix[1]);

#else // DIM == 3
    uint64_t ix[3];
    for (int d = 0; d < 3; ++d) {
        real normalized = (pos[d] - domain_min[d]) / domain_size[d];
        normalized = std::max(0.0, std::min(1.0, normalized));
        ix[d] = static_cast<uint64_t>(normalized * max_coord);
    }
    return encode_3d(ix[0], ix[1], ix[2]);
#endif
}

// Structure to hold Morton code and original index
struct MortonPair {
    uint64_t code;
    int index;

    bool operator<(const MortonPair& other) const {
        return code < other.code;
    }
};

// Compute Morton codes for all particles and return sorted permutation
// Returns vector where result[i] = original index of particle that should be at position i
inline std::vector<int> compute_sort_permutation(
    const std::vector<SPHParticle>& particles,
    const real domain_min[DIM],
    const real domain_size[DIM])
{
    const int n = static_cast<int>(particles.size());
    if (n == 0) return {};

    // Compute Morton codes
    std::vector<MortonPair> morton_pairs(n);

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        morton_pairs[i].code = encode(particles[i].pos.get_array(), domain_min, domain_size);
        morton_pairs[i].index = i;
    }

    // Sort by Morton code
    std::sort(morton_pairs.begin(), morton_pairs.end());

    // Extract permutation
    std::vector<int> permutation(n);
    for (int i = 0; i < n; ++i) {
        permutation[i] = morton_pairs[i].index;
    }

    return permutation;
}

// Apply permutation to reorder particles in-place
// After this, particles[i] will contain the particle that was at permutation[i]
inline void apply_permutation(std::vector<SPHParticle>& particles,
                              const std::vector<int>& permutation)
{
    const int n = static_cast<int>(particles.size());
    if (n == 0) return;

    // Create reordered copy
    std::vector<SPHParticle> reordered(n);

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        reordered[i] = particles[permutation[i]];
        // Update particle ID to reflect new position
        reordered[i].id = i;
    }

    // Swap back
    particles.swap(reordered);
}

// Compute inverse permutation
// If permutation[i] = j, then inverse[j] = i
// This maps old indices to new indices
inline std::vector<int> compute_inverse_permutation(const std::vector<int>& permutation)
{
    const int n = static_cast<int>(permutation.size());
    std::vector<int> inverse(n);

    for (int i = 0; i < n; ++i) {
        inverse[permutation[i]] = i;
    }

    return inverse;
}

// Sort particles by Morton code (main entry point)
// This reorders particles for cache-friendly access during tree traversal
// Returns the permutation that was applied (for debugging/validation)
inline std::vector<int> sort_particles_by_morton(
    std::vector<SPHParticle>& particles,
    const real domain_min[DIM],
    const real domain_size[DIM])
{
    auto permutation = compute_sort_permutation(particles, domain_min, domain_size);
    apply_permutation(particles, permutation);
    return permutation;
}

} // namespace morton
} // namespace sph
