// Benchmark comparison: Optimized vs Baseline
// This test measures performance WITH optimizations enabled
// To compare, build with/without optimization flags

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

class BenchmarkComparison : public ::testing::Test {
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
                real r = std::sqrt(r2);
                real h = 0.01 + 0.03 * (r / radius);
                particles.push_back(create_particle(id++, x, y, z, h));
            }
        }
    }
};

// Main benchmark: Full tree build + all neighbor searches
// This measures what matters for actual simulation performance
TEST_F(BenchmarkComparison, SimulationStep) {
    const int N = 100000;  // Fixed particle count for consistent comparison
    const int iterations = 10;

    create_sphere_particles(N);
    tree->resize(particles.size());

    // Warm up
    tree->make(particles, particles.size());
    tree->set_kernel();
    #pragma omp parallel
    {
        std::vector<int> local_neighbor_list(neighbor_list_size);
        #pragma omp for
        for (int i = 0; i < N; ++i) {
            tree->neighbor_search(particles[i], local_neighbor_list, particles, true);
        }
    }

    // Benchmark: tree construction + all neighbor searches
    auto start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < iterations; ++iter) {
        tree->make(particles, particles.size());
        tree->set_kernel();
        #pragma omp parallel
        {
            std::vector<int> local_neighbor_list(neighbor_list_size);
            #pragma omp for
            for (int i = 0; i < N; ++i) {
                tree->neighbor_search(particles[i], local_neighbor_list, particles, true);
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();

    double total_ms = std::chrono::duration<double, std::milli>(end - start).count();
    double per_step_ms = total_ms / iterations;

    std::cout << "\n╔══════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║         PERFORMANCE BENCHMARK RESULTS                        ║" << std::endl;
    std::cout << "╠══════════════════════════════════════════════════════════════╣" << std::endl;
    std::cout << "║ Particles: " << std::setw(10) << N << "                                     ║" << std::endl;
    std::cout << "║ Iterations: " << std::setw(9) << iterations << "                                     ║" << std::endl;
    std::cout << "║ Time per step: " << std::fixed << std::setprecision(2) << std::setw(8) << per_step_ms << " ms                              ║" << std::endl;
    std::cout << "╠══════════════════════════════════════════════════════════════╣" << std::endl;
    std::cout << "║ Optimization Flags:                                          ║" << std::endl;
#ifdef USE_MORTON_ORDERING
    std::cout << "║   [✓] Morton Ordering                                        ║" << std::endl;
#else
    std::cout << "║   [✗] Morton Ordering                                        ║" << std::endl;
#endif
#ifdef USE_ITERATIVE_TRAVERSAL
    std::cout << "║   [✓] Iterative Traversal                                    ║" << std::endl;
#else
    std::cout << "║   [✗] Iterative Traversal                                    ║" << std::endl;
#endif
#ifdef SPH_USE_THREAD_LOCAL_NEIGHBOR_LIST
    std::cout << "║   [✓] Thread-Local Neighbor List                             ║" << std::endl;
#else
    std::cout << "║   [✗] Thread-Local Neighbor List                             ║" << std::endl;
#endif
#ifdef SPH_USE_DISTANCE_CACHING
    std::cout << "║   [✓] Distance Caching                                       ║" << std::endl;
#else
    std::cout << "║   [✗] Distance Caching                                       ║" << std::endl;
#endif
#ifdef SPH_USE_DYNAMIC_SCHEDULING
    std::cout << "║   [✓] Dynamic Scheduling                                     ║" << std::endl;
#else
    std::cout << "║   [✗] Dynamic Scheduling                                     ║" << std::endl;
#endif
    std::cout << "╚══════════════════════════════════════════════════════════════╝" << std::endl;

    SUCCEED();
}
