// Performance benchmark for SPH simulation optimizations
// Measures tree construction and neighbor search timing

#include <gtest/gtest.h>
#include <vector>
#include <random>
#include <chrono>
#include <iostream>
#include <iomanip>

#include "defines.hpp"
#include "particle.hpp"
#include "vector_type.hpp"
#include "bhtree.hpp"
#include "parameters.hpp"

using namespace sph;

class PerformanceBenchmark : public ::testing::Test {
protected:
    std::shared_ptr<SPHParameters> param;
    std::unique_ptr<BHTree> tree;
    std::vector<SPHParticle> particles;

    void SetUp() override {
        param = std::make_shared<SPHParameters>();
        param->physics.gamma = 1.6667;
        param->physics.neighbor_number = 50;
        param->tree.max_level = 20;
        param->tree.leaf_particle_num = 16;
        param->gravity.constant = 0.00430091;
        param->gravity.theta = 0.5;

        tree = std::make_unique<BHTree>();
        tree->initialize(param);
    }

    SPHParticle create_particle(int id, real x, real y, real z, real h = 0.02) {
        SPHParticle p;
        p.id = id;
        p.pos[0] = x;
        p.pos[1] = y;
        p.pos[2] = z;
        p.vel[0] = p.vel[1] = p.vel[2] = 0.0;
        p.mass = 1.0;
        p.sml = h;
        p.dens = 1.0;
        p.pres = 1.0;
        p.ene = 1.0;
        p.next = nullptr;
        return p;
    }

    // Create particles in a sphere with density profile (like Bonnor-Ebert)
    void create_sphere_particles(int n, real radius = 1.0, unsigned seed = 42) {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<real> dist(-radius, radius);

        particles.clear();
        particles.reserve(n);

        int id = 0;
        while (static_cast<int>(particles.size()) < n) {
            real x = dist(rng);
            real y = dist(rng);
            real z = dist(rng);
            real r2 = x*x + y*y + z*z;
            if (r2 < radius * radius) {
                // Variable smoothing length based on local density (denser at center)
                real r = std::sqrt(r2);
                real h = 0.01 + 0.03 * (r / radius);  // h ranges from 0.01 to 0.04
                particles.push_back(create_particle(id++, x, y, z, h));
            }
        }
    }

    // Create particles on a random distribution
    void create_random_particles(int n, real box_size = 1.0, unsigned seed = 42) {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<real> dist(0.0, box_size);
        std::uniform_real_distribution<real> h_dist(0.01, 0.05);

        particles.clear();
        particles.reserve(n);
        for (int i = 0; i < n; ++i) {
            real h = h_dist(rng);
            particles.push_back(create_particle(i, dist(rng), dist(rng), dist(rng), h));
        }
    }
};

// Benchmark tree construction
TEST_F(PerformanceBenchmark, TreeConstruction) {
    const std::vector<int> particle_counts = {10000, 50000, 100000, 200000};

    std::cout << "\n=== Tree Construction Benchmark ===" << std::endl;
    std::cout << std::setw(15) << "Particles" << std::setw(15) << "Time (ms)"
              << std::setw(15) << "Rate (M/s)" << std::endl;
    std::cout << std::string(45, '-') << std::endl;

    for (int n : particle_counts) {
        create_sphere_particles(n);
        tree->resize(particles.size());

        // Warm up
        tree->make(particles, particles.size());
        tree->set_kernel();

        // Benchmark
        const int iterations = 5;
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; ++i) {
            tree->make(particles, particles.size());
            tree->set_kernel();
        }
        auto end = std::chrono::high_resolution_clock::now();

        double ms = std::chrono::duration<double, std::milli>(end - start).count() / iterations;
        double rate = (n / 1e6) / (ms / 1000.0);

        std::cout << std::setw(15) << n << std::setw(15) << std::fixed << std::setprecision(2) << ms
                  << std::setw(15) << rate << std::endl;
    }
}

// Benchmark neighbor search
TEST_F(PerformanceBenchmark, NeighborSearch) {
    const std::vector<int> particle_counts = {10000, 50000, 100000, 200000};

    std::cout << "\n=== Neighbor Search Benchmark ===" << std::endl;
    std::cout << std::setw(15) << "Particles" << std::setw(15) << "Time (ms)"
              << std::setw(15) << "Searches/s" << std::setw(15) << "Avg Neighbors" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    for (int n : particle_counts) {
        create_sphere_particles(n);
        tree->resize(particles.size());
        tree->make(particles, particles.size());
        tree->set_kernel();

        std::vector<int> neighbor_list(neighbor_list_size);

        // Sample particles for neighbor search
        const int num_samples = std::min(1000, n);
        std::vector<int> sample_indices;
        for (int i = 0; i < num_samples; ++i) {
            sample_indices.push_back((i * n) / num_samples);
        }

        // Warm up
        for (int idx : sample_indices) {
            tree->neighbor_search(particles[idx], neighbor_list, particles, true);
        }

        // Benchmark
        const int iterations = 3;
        long long total_neighbors = 0;
        auto start = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < iterations; ++iter) {
            for (int idx : sample_indices) {
                int nn = tree->neighbor_search(particles[idx], neighbor_list, particles, true);
                if (iter == 0) total_neighbors += nn;
            }
        }
        auto end = std::chrono::high_resolution_clock::now();

        double ms = std::chrono::duration<double, std::milli>(end - start).count();
        double searches_per_sec = (num_samples * iterations) / (ms / 1000.0);
        double avg_neighbors = static_cast<double>(total_neighbors) / num_samples;

        std::cout << std::setw(15) << n << std::setw(15) << std::fixed << std::setprecision(2) << ms
                  << std::setw(15) << std::setprecision(0) << searches_per_sec
                  << std::setw(15) << std::setprecision(1) << avg_neighbors << std::endl;
    }
}

// Full simulation step benchmark (tree + all neighbor searches)
TEST_F(PerformanceBenchmark, FullNeighborSearchPass) {
    const std::vector<int> particle_counts = {10000, 50000, 100000};

    std::cout << "\n=== Full Neighbor Search Pass Benchmark ===" << std::endl;
    std::cout << "(All particles, with OpenMP parallelization)" << std::endl;
    std::cout << std::setw(15) << "Particles" << std::setw(15) << "Time (ms)"
              << std::setw(15) << "Rate (kp/s)" << std::endl;
    std::cout << std::string(45, '-') << std::endl;

    for (int n : particle_counts) {
        create_sphere_particles(n);
        tree->resize(particles.size());
        tree->make(particles, particles.size());
        tree->set_kernel();

        // Warm up
        #pragma omp parallel
        {
            std::vector<int> local_neighbor_list(neighbor_list_size);
            #pragma omp for
            for (int i = 0; i < n; ++i) {
                tree->neighbor_search(particles[i], local_neighbor_list, particles, true);
            }
        }

        // Benchmark
        const int iterations = 3;
        auto start = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < iterations; ++iter) {
            #pragma omp parallel
            {
                std::vector<int> local_neighbor_list(neighbor_list_size);
                #pragma omp for
                for (int i = 0; i < n; ++i) {
                    tree->neighbor_search(particles[i], local_neighbor_list, particles, true);
                }
            }
        }
        auto end = std::chrono::high_resolution_clock::now();

        double ms = std::chrono::duration<double, std::milli>(end - start).count() / iterations;
        double rate = (n / 1e3) / (ms / 1000.0);

        std::cout << std::setw(15) << n << std::setw(15) << std::fixed << std::setprecision(2) << ms
                  << std::setw(15) << std::setprecision(0) << rate << std::endl;
    }
}

// Report optimization flags status
TEST_F(PerformanceBenchmark, OptimizationStatus) {
    std::cout << "\n=== Optimization Flags Status ===" << std::endl;

#ifdef USE_MORTON_ORDERING
    std::cout << "✓ USE_MORTON_ORDERING: ENABLED" << std::endl;
#else
    std::cout << "✗ USE_MORTON_ORDERING: DISABLED" << std::endl;
#endif

#ifdef USE_ITERATIVE_TRAVERSAL
    std::cout << "✓ USE_ITERATIVE_TRAVERSAL: ENABLED" << std::endl;
#else
    std::cout << "✗ USE_ITERATIVE_TRAVERSAL: DISABLED" << std::endl;
#endif

#ifdef SPH_USE_THREAD_LOCAL_NEIGHBOR_LIST
    std::cout << "✓ SPH_USE_THREAD_LOCAL_NEIGHBOR_LIST: ENABLED" << std::endl;
#else
    std::cout << "✗ SPH_USE_THREAD_LOCAL_NEIGHBOR_LIST: DISABLED" << std::endl;
#endif

#ifdef SPH_USE_DISTANCE_CACHING
    std::cout << "✓ SPH_USE_DISTANCE_CACHING: ENABLED" << std::endl;
#else
    std::cout << "✗ SPH_USE_DISTANCE_CACHING: DISABLED" << std::endl;
#endif

#ifdef SPH_USE_DYNAMIC_SCHEDULING
    std::cout << "✓ SPH_USE_DYNAMIC_SCHEDULING: ENABLED" << std::endl;
#else
    std::cout << "✗ SPH_USE_DYNAMIC_SCHEDULING: DISABLED" << std::endl;
#endif

    std::cout << "\nNeighbor list size: " << neighbor_list_size << std::endl;
    std::cout << "Dimensions: " << DIM << std::endl;

    SUCCEED();
}
