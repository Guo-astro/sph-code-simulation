#include <cassert>
#include <iostream>
#include <array>

#include "parameters.hpp"
#include "bhtree.hpp"
#include "openmp.hpp"
#include "exception.hpp"
#include "periodic.hpp"

#ifdef USE_MORTON_ORDERING
#include "morton.hpp"
#endif

namespace sph
{

void BHTree::initialize(std::shared_ptr<SPHParameters> param)
{
    m_max_level         = param->tree.max_level;
    m_leaf_particle_num = param->tree.leaf_particle_num;
    m_root.clear();
    m_root.level = 1;
    m_is_periodic = param->periodic.is_valid;
    if(m_is_periodic) {
        m_range_max = param->periodic.range_max;
        m_range_min = param->periodic.range_min;
        m_root.center = (m_range_max + m_range_min) * 0.5;
        auto range = m_range_max - m_range_min;
        real l = 0.0;
        for(int i = 0; i < DIM; ++i) {
            if(l < range[i]) {
                l = range[i];
            }
        }
        m_root.edge = l;
    }
    m_periodic = std::make_shared<Periodic>();
    m_periodic->initialize(param);

    if(param->gravity.is_valid) {
        m_g_constant = param->gravity.constant;
        m_theta = param->gravity.theta;
        m_theta2 = m_theta * m_theta;
        m_softening_type = param->gravity.softening_type;
        m_use_fixed_softening = param->gravity.use_fixed_softening;
        m_fixed_softening = param->gravity.fixed_softening;
    }
}

void BHTree::resize(const int particle_num, const int tree_size)
{
    // Allow re-initialization by resetting nodes if already allocated
    if(m_nodes.get() != nullptr) {
        m_nodes.reset();
    }

    m_node_size = particle_num * tree_size;
    m_nodes = std::shared_ptr<BHNode>(new BHNode[m_node_size], std::default_delete<BHNode[]>());

#pragma omp parallel for
    for(int i = 0; i < m_node_size; ++i) {
        m_nodes.get()[i].clear();
    }
}

void BHTree::make(std::vector<SPHParticle> & particles, const int particle_num)
{
    m_root.root_clear();

    if(particle_num == 0) {
        return;  // No particles to process
    }

    if(!m_is_periodic) {
        omp_real r_min[DIM];
        omp_real r_max[DIM];
        for(int i = 0; i < DIM; ++i) {
            r_min[i].get() = std::numeric_limits<real>::max();
            r_max[i].get() = std::numeric_limits<real>::lowest();
        }

#pragma omp parallel for
        for(int i = 0; i < particle_num; ++i) {
            auto & r_i = particles[i].pos;
            for(int j = 0; j < DIM; ++j) {
                if(r_min[j].get() > r_i[j]) {
                    r_min[j].get() = r_i[j];
                }
                if(r_max[j].get() < r_i[j]) {
                    r_max[j].get() = r_i[j];
                }
            }
        }

        vec_t range_min, range_max;
        for(int i = 0; i < DIM; ++i) {
            range_min[i] = r_min[i].min();
            range_max[i] = r_max[i].max();
        }

        m_root.center = (range_max + range_min) * 0.5;
        auto range = range_max - range_min;
        real l = 0.0;
        for(int i = 0; i < DIM; ++i) {
            if(l < range[i]) {
                l = range[i];
            }
        }
        m_root.edge = l;
    }

#pragma omp parallel for
    for(int i = 0; i < particle_num - 1; ++i) {
        particles[i].next = &particles[i + 1];
    }
    particles[particle_num - 1].next = nullptr;
    m_root.first = &particles[0];

    int remaind = m_node_size;
    auto * p = m_nodes.get();
    m_root.create_tree(p, remaind, m_max_level, m_leaf_particle_num);
}

void BHTree::set_kernel()
{
    m_root.set_kernel();
}

int BHTree::neighbor_search(const SPHParticle & p_i, std::vector<int> & neighbor_list, const std::vector<SPHParticle> & particles, const bool is_ij)
{
    int n_neighbor = 0;
    int max_neighbors = neighbor_list.size();
    m_root.neighbor_search(p_i, neighbor_list, n_neighbor, max_neighbors, is_ij, m_periodic.get());

    // CRITICAL DEBUG: Log the values to see what's happening
    if(n_neighbor > max_neighbors || n_neighbor < 0) {
        std::cerr << "[DEBUG] OVERFLOW! particle=" << p_i.id 
                  << " n_neighbor=" << n_neighbor 
                  << " max_neighbors=" << max_neighbors 
                  << " list.size()=" << neighbor_list.size() << std::endl;
        THROW_ERROR("n_neighbor (", n_neighbor, ") exceeds max_neighbors (", max_neighbors, ") for particle ", p_i.id);
    }

    const auto & pos_i = p_i.pos;
    std::sort(neighbor_list.begin(), neighbor_list.begin() + n_neighbor, [&](const int a, const int b) {
        const vec_t r_ia = m_periodic->calc_r_ij(pos_i, particles[a].pos);
        const vec_t r_ib = m_periodic->calc_r_ij(pos_i, particles[b].pos);
        return abs2(r_ia) < abs2(r_ib);
    });
    return n_neighbor;
}

void BHTree::tree_force(SPHParticle & p_i)
{
    p_i.phi = 0.0;
    p_i.grav_acc = vec_t(0.0);  // Initialize gravity acceleration
    m_root.calc_force(p_i, m_theta2, m_g_constant, m_periodic.get(),
                      m_softening_type, m_use_fixed_softening, m_fixed_softening);
    // Note: grav_acc is stored but NOT added to acc here.
    // This allows the gravity-aware Riemann solver to use grav_acc,
    // and then gravity is added to acc AFTER fluid force calculation.
}

// --------------------------------------------------------------------------------------------------------------- //

void BHTree::BHNode::create_tree(BHNode * & nodes, int & remaind, const int max_level, const int leaf_particle_num)
{
    std::fill(childs, childs + NCHILD, nullptr);

    auto * pp = first;
    do {
        auto * pnext = pp->next;
        assign(pp, nodes, remaind);
        pp = pnext;
    } while(pp != nullptr);

    int num_child = 0;
    real mass_before = mass;  // DEBUG
    for(int i = 0; i < NCHILD; ++i) {
        auto * child = childs[i];
        if(child) {
            ++num_child;
            
            // First, recursively build the tree or mark as leaf
            if(child->num > leaf_particle_num && level < max_level) {
                // CRITICAL FIX: Clear child's mass before recursion
                // because assign() already added mass, but recursion will
                // recalculate it from grandchildren
                child->mass = 0.0;
                child->m_center = 0.0;
                child->create_tree(nodes, remaind, max_level, leaf_particle_num);
                // After recursion, child->m_center is already the true center of mass
            } else {
                child->is_leaf = true;
                // For leaf nodes, m_center was accumulated as sum(mass*pos) in assign()
                // We need to normalize it to get the true center of mass
                if(child->mass > 0.0) {
                    child->m_center /= child->mass;
                }
            }
            
            // BUG FIX: Accumulate child mass into parent
            // Without this, root and internal nodes have mass=0, breaking gravity!
            mass += child->mass;
            
            // Now accumulate the child's center of mass (already normalized)
            // multiplied by its mass to get mass-weighted sum
            m_center += child->m_center * child->mass;
        }
    }
    
    // Finally, normalize parent's accumulated mass-weighted center
    if(mass > 0.0) {
        m_center /= mass;
    }
}

void BHTree::BHNode::assign(SPHParticle * particle, BHNode * & nodes, int & remaind)
{
    auto & p_i = *particle;
    const auto & pos = p_i.pos;

    int index = 0;
    int mask = 1;
    for(int i = 0; i < DIM; ++i) {
        if(pos[i] > center[i]) {
            index |= mask;
        }
        mask <<= 1;
    }

    auto * child = childs[index];
    if(!child) {
        if(remaind < 0) {
            THROW_ERROR("There is no free node.");
        }
        childs[index] = nodes;
        child = childs[index];
        ++nodes;
        --remaind;
        child->clear();
        child->level = level + 1;
        child->edge = edge * 0.5;

        int a = 1;
        real b = 2.0;
        for(int i = 0; i < DIM; ++i) {
            child->center[i] = center[i] + ((index & a) * b - 1.0) * edge * 0.25;
            a <<= 1;
            b *= 0.5;
        }
    }

    child->num++;
    child->mass += p_i.mass;
    child->m_center += pos * p_i.mass;
    p_i.next = child->first;
    child->first = particle;
}

real BHTree::BHNode::set_kernel()
{
    real kernel = 0.0;
    if(is_leaf) {
        auto * p = first;
        while(p) {
            const real h = p->sml;
            if(h > kernel) {
                kernel = h;
            }
            p = p->next;
        }
    } else {
        for(int i = 0; i < NCHILD; ++i) {
            auto * child = childs[i];
            if(child) {
                const real h = child->set_kernel();
                if(h > kernel) {
                    kernel = h;
                }
            }
        }
    }

    kernel_size = kernel;
    return kernel_size;
}

void BHTree::BHNode::neighbor_search(const SPHParticle & p_i, std::vector<int> & neighbor_list, int & n_neighbor, int max_neighbors, const bool is_ij, const Periodic * periodic)
{
    const vec_t & r_i = p_i.pos;
    const real h = is_ij ? std::max(p_i.sml, kernel_size) : p_i.sml;
    const real h2 = h * h;
    const real l2 = sqr(edge * 0.5 + h);
    const vec_t d = periodic->calc_r_ij(r_i, center);
    real dx2_max = sqr(d[0]);
    for(int i = 1; i < DIM; ++i) {
        const real dx2 = sqr(d[i]);
        if(dx2 > dx2_max) {
            dx2_max = dx2;
        }
    }

    if(dx2_max <= l2) {
        if(is_leaf) {
            auto * p = first;
            while(p) {
                const vec_t & r_j = p->pos;
                const vec_t r_ij = periodic->calc_r_ij(r_i, r_j);
                const real r2 = abs2(r_ij);
                if(r2 < h2) {
                    if(n_neighbor >= max_neighbors) {
                        THROW_ERROR("Neighbor list overflow: increase neighbor_list_size in defines.hpp");
                    }
                    neighbor_list[n_neighbor] = p->id;
                    ++n_neighbor;
                }
                p = p->next;
            }
        } else {
            for(int i = 0; i < NCHILD; ++i) {
                if(childs[i]) {
                    childs[i]->neighbor_search(p_i, neighbor_list, n_neighbor, max_neighbors, is_ij, periodic);
                }
            }
        }
    }
}

 // Hernquist & Katz (1989)
inline real f(const real r, const real h)
{
    const real e = h * 0.5;
    const real u = r / e;
    
    if(u < 1.0) {
        return (-0.5 * u * u * (1.0 / 3.0 - 3.0 / 20 * u * u + u * u * u / 20) + 1.4) / e;
    } else if(u < 2.0) {
        return -1.0 / (15 * r) + (-u * u * (4.0 / 3.0 - u + 0.3 * u * u - u * u * u / 30) + 1.6) / e;
    } else {
        return 1 / r;
    }
}

inline real g(const real r, const real h)
{
    const real e = h * 0.5;
    const real u = r / e;
    
    if(u < 1.0) {
        return (4.0 / 3.0 - 1.2 * u * u + 0.5 * u * u * u) / (e * e * e);
    } else if(u < 2.0) {
        return (-1.0 / 15 + 8.0 / 3 * u * u * u - 3 * u * u * u * u + 1.2 * u * u * u * u * u - u * u * u * u * u * u / 6.0) / (r * r * r);
    } else {
        return 1 / (r * r * r);
    }
}

// Wendland C4 gravitational potential kernel (3D only)
// Derived from solving ∇²φ̃ = -4πG W for Wendland C4 kernel
// The potential is: φ(r) = -G ∫ M(<r')/r'² dr' where M(<r') is enclosed mass
// For W_C4(q) = (495/32π)/h³ (1-q)⁶(1 + 6q + 35/3 q²), q = r/h, 0 ≤ q ≤ 1
// Numerically integrated to get polynomial fit: φ(q)*h = a₀ + a₁q + a₂q² + ... + a₉q⁹
inline real wendland_phi(const real r, const real h)
{
#if DIM == 3
    const real q = r / h;
    if (q >= 1.0) {
        return 1.0 / r;
    }
    const real q2 = q * q;
    const real q3 = q2 * q;
    const real q4 = q2 * q2;
    const real q5 = q4 * q;
    const real q6 = q3 * q3;
    const real q7 = q6 * q;
    const real q8 = q4 * q4;
    const real q9 = q8 * q;
    
    // Coefficients from numerical integration of Poisson equation (verified)
    // φ(0)*h ≈ 3.44 corresponds to enclosed mass integral
    const real a0 =  3.4374743761;
    const real a1 = -0.0031873250;  // ≈ 0 (boundary condition)
    const real a2 = -10.2154807743;
    const real a3 = -1.1577720555;
    const real a4 =  36.1013669755;
    const real a5 = -26.3399094060;
    const real a6 = -44.1079372114;
    const real a7 =  82.6543766683;
    const real a8 = -50.5921624056;
    const real a9 =  11.2232565249;
    
    return (a0 + a1*q + a2*q2 + a3*q3 + a4*q4 + a5*q5 + a6*q6 + a7*q7 + a8*q8 + a9*q9) / h;
#else
    return 1.0 / (r + 1e-10);
#endif
}

// Wendland C4 gravitational force kernel: g(r) = -dφ/dr / r
// This is derived from the derivative of wendland_phi
// Force: F = -G m₁ m₂ g(r) r̂
// g(r) = (1/h³) × [-(d/dq)(φ*h)/q] where the derivative gives the force direction
inline real wendland_g(const real r, const real h)
{
#if DIM == 3
    const real q = r / h;
    if (q >= 1.0) {
        return 1.0 / (r * r * r);
    }
    if (q < 1e-10) {
        // At q=0, the force is zero (symmetric mass distribution)
        return 0.0;
    }
    const real q2 = q * q;
    const real q3 = q2 * q;
    const real q4 = q2 * q2;
    const real q5 = q4 * q;
    const real q6 = q3 * q3;
    const real q7 = q6 * q;
    
    // Derivative coefficients: bₙ = n × aₙ (from d(φ*h)/dq = b₁ + b₂q + ...)
    // Then g(r) = -(b₁ + b₂q + b₃q² + ... + b₉q⁸) / (h³ × q)
    const real b1 = -0.0031873250;   // 1 × a1 ≈ 0
    const real b2 = -20.4309615486;  // 2 × a2
    const real b3 = -3.4733161665;   // 3 × a3
    const real b4 = 144.4054679020;  // 4 × a4
    const real b5 = -131.6995470300; // 5 × a5
    const real b6 = -264.6476232684; // 6 × a6
    const real b7 = 578.5806366781;  // 7 × a7
    const real b8 = -404.7372992448; // 8 × a8
    const real b9 = 101.0093087241;  // 9 × a9
    
    const real denom = h * h * h;
    // g(r) = -(b₁/q + b₂ + b₃q + b₄q² + ... + b₉q⁷) / h³
    return -(b1/q + b2 + b3*q + b4*q2 + b5*q3 + b6*q4 + b7*q5 + b8*q6 + b9*q7) / denom;
#else
    return 1.0 / (r * r * r + 1e-30);
#endif
}

void BHTree::BHNode::calc_force(SPHParticle & p_i, const real theta2, const real g_constant, const Periodic * periodic,
                                GravitySofteningType softening_type, bool use_fixed_softening, real fixed_softening)
{
    const vec_t & r_i = p_i.pos;
    const real l2 = edge * edge;
    const vec_t d = periodic->calc_r_ij(r_i, m_center);
    const real d2 = abs2(d);

    // Check if we need to open this node:
    // 1. Barnes-Hut criterion: l^2 > theta^2 * d^2 (cell not small enough)
    // 2. Distance safety: d < l (particle might be inside cell)
    const bool theta_open = (l2 > theta2 * d2);
    const bool distance_open = (d2 < l2);
    
    if(theta_open || distance_open) {
        // Must open node - either due to Barnes-Hut or safety
        if(is_leaf) {
            // Leaf node: compute direct pairwise forces
            auto * p = first;
            while(p) {
                if(p->id != p_i.id) {  // Skip self-interaction
                    const vec_t & r_j = p->pos;
                    const vec_t r_ij = periodic->calc_r_ij(r_i, r_j);
                    const real r = std::abs(r_ij);
                    
                    if (softening_type == GravitySofteningType::WENDLAND_C4) {
                        if (use_fixed_softening) {
                            p_i.phi -= g_constant * p->mass * wendland_phi(r, fixed_softening);
                            p_i.grav_acc -= r_ij * (g_constant * p->mass * wendland_g(r, fixed_softening));
                        } else {
                            const real h_ij = 0.5 * (p_i.sml + p->sml);
                            p_i.phi -= g_constant * p->mass * wendland_phi(r, h_ij);
                            p_i.grav_acc -= r_ij * (g_constant * p->mass * wendland_g(r, h_ij));
                        }
                    } else {
                        // Hernquist-Katz (default)
                        if (use_fixed_softening) {
                            const real h_fixed = fixed_softening * 2.0;
                            p_i.phi -= g_constant * p->mass * f(r, h_fixed);
                            p_i.grav_acc -= r_ij * (g_constant * p->mass * g(r, h_fixed));
                        } else {
                            p_i.phi -= g_constant * p->mass * (f(r, p_i.sml) + f(r, p->sml)) * 0.5;
                            p_i.grav_acc -= r_ij * (g_constant * p->mass * (g(r, p_i.sml) + g(r, p->sml)) * 0.5);
                        }
                    }
                }
                p = p->next;
            }
        } else {
            // Internal node: recurse to children
            for(int i = 0; i < NCHILD; ++i) {
                if(childs[i]) {
                    childs[i]->calc_force(p_i, theta2, g_constant, periodic, softening_type, use_fixed_softening, fixed_softening);
                }
            }
        }
    } else {
        // Can use monopole approximation safely
        const real r_inv = 1.0 / std::sqrt(d2);
        p_i.phi -= g_constant * mass * r_inv;
        p_i.grav_acc -= d * (g_constant * mass * pow3(r_inv));
    }
}

// =============================================================================
// Iterative Tree Traversal Implementations
// =============================================================================

#ifdef USE_ITERATIVE_TRAVERSAL

// Maximum tree depth of 64 levels supports 2^64 nodes
constexpr int ITERATIVE_MAX_TREE_DEPTH = 64;

int BHTree::neighbor_search_iterative(const SPHParticle & p_i, std::vector<int> & neighbor_list,
                                       const std::vector<SPHParticle> & particles, const bool is_ij)
{
    int n_neighbor = 0;
    const int max_neighbors = static_cast<int>(neighbor_list.size());
    const vec_t & r_i = p_i.pos;

    // Thread-local stack to avoid recursion
    thread_local std::array<BHNode*, ITERATIVE_MAX_TREE_DEPTH> stack;
    int stack_top = 0;

    // Start with root node
    stack[stack_top++] = &m_root;

    while (stack_top > 0) {
        BHNode* node = stack[--stack_top];

        // Compute kernel size for this search
        const real h = is_ij ? std::max(p_i.sml, node->kernel_size) : p_i.sml;
        const real h2 = h * h;
        const real l2 = sqr(node->edge * 0.5 + h);

        // AABB overlap test: check if search sphere overlaps node bounding box
        const vec_t d = m_periodic->calc_r_ij(r_i, node->center);
        real dx2_max = sqr(d[0]);
        for (int i = 1; i < DIM; ++i) {
            const real dx2 = sqr(d[i]);
            if (dx2 > dx2_max) {
                dx2_max = dx2;
            }
        }

        if (dx2_max > l2) {
            continue;  // No overlap, skip this node
        }

        if (node->is_leaf) {
            // Process particles in leaf node
            auto * p = node->first;
            while (p) {
                // Scalar distance check with periodic boundary handling
                const vec_t & r_j = p->pos;
                const vec_t r_ij = m_periodic->calc_r_ij(r_i, r_j);
                const real r2 = abs2(r_ij);
                if (r2 < h2) {
                    if (n_neighbor >= max_neighbors) {
                        THROW_ERROR("Neighbor list overflow: increase neighbor_list_size in defines.hpp");
                    }
                    neighbor_list[n_neighbor++] = p->id;
                }
                p = p->next;
            }
        } else {
            // Push children onto stack (reverse order for depth-first traversal)
            for (int c = NCHILD - 1; c >= 0; --c) {
                if (node->childs[c]) {
                    if (stack_top >= ITERATIVE_MAX_TREE_DEPTH) {
                        THROW_ERROR("Tree depth exceeded ITERATIVE_MAX_TREE_DEPTH (", ITERATIVE_MAX_TREE_DEPTH, ")");
                    }
                    stack[stack_top++] = node->childs[c];
                }
            }
        }
    }

    // Overflow check
    if (n_neighbor > max_neighbors || n_neighbor < 0) {
        THROW_ERROR("n_neighbor (", n_neighbor, ") exceeds max_neighbors (", max_neighbors, ") for particle ", p_i.id);
    }

    // Sort neighbors by distance (same as recursive version)
    const auto & pos_i = p_i.pos;
    std::sort(neighbor_list.begin(), neighbor_list.begin() + n_neighbor, [&](const int a, const int b) {
        const vec_t r_ia = m_periodic->calc_r_ij(pos_i, particles[a].pos);
        const vec_t r_ib = m_periodic->calc_r_ij(pos_i, particles[b].pos);
        return abs2(r_ia) < abs2(r_ib);
    });

    return n_neighbor;
}

void BHTree::tree_force_iterative(SPHParticle & p_i)
{
    p_i.phi = 0.0;
    p_i.grav_acc = vec_t(0.0);

    const vec_t & r_i = p_i.pos;

    // Thread-local stack
    thread_local std::array<BHNode*, ITERATIVE_MAX_TREE_DEPTH> stack;
    int stack_top = 0;

    stack[stack_top++] = &m_root;

    while (stack_top > 0) {
        BHNode* node = stack[--stack_top];

        const real l2 = node->edge * node->edge;
        const vec_t d = m_periodic->calc_r_ij(r_i, node->m_center);
        const real d2 = abs2(d);

        // Barnes-Hut opening criterion
        const bool theta_open = (l2 > m_theta2 * d2);
        const bool distance_open = (d2 < l2);

        if (theta_open || distance_open) {
            // Must open node
            if (node->is_leaf) {
                // Direct pairwise forces
                auto * p = node->first;
                while (p) {
                    if (p->id != p_i.id) {
                        const vec_t & r_j = p->pos;
                        const vec_t r_ij = m_periodic->calc_r_ij(r_i, r_j);
                        const real r = std::abs(r_ij);

                        if (m_softening_type == GravitySofteningType::WENDLAND_C4) {
                            if (m_use_fixed_softening) {
                                p_i.phi -= m_g_constant * p->mass * wendland_phi(r, m_fixed_softening);
                                p_i.grav_acc -= r_ij * (m_g_constant * p->mass * wendland_g(r, m_fixed_softening));
                            } else {
                                const real h_ij = 0.5 * (p_i.sml + p->sml);
                                p_i.phi -= m_g_constant * p->mass * wendland_phi(r, h_ij);
                                p_i.grav_acc -= r_ij * (m_g_constant * p->mass * wendland_g(r, h_ij));
                            }
                        } else {
                            // Hernquist-Katz
                            if (m_use_fixed_softening) {
                                const real h_fixed = m_fixed_softening * 2.0;
                                p_i.phi -= m_g_constant * p->mass * f(r, h_fixed);
                                p_i.grav_acc -= r_ij * (m_g_constant * p->mass * g(r, h_fixed));
                            } else {
                                p_i.phi -= m_g_constant * p->mass * (f(r, p_i.sml) + f(r, p->sml)) * 0.5;
                                p_i.grav_acc -= r_ij * (m_g_constant * p->mass * (g(r, p_i.sml) + g(r, p->sml)) * 0.5);
                            }
                        }
                    }
                    p = p->next;
                }
            } else {
                // Push children (reverse order for depth-first)
                for (int c = NCHILD - 1; c >= 0; --c) {
                    if (node->childs[c]) {
                        if (stack_top >= ITERATIVE_MAX_TREE_DEPTH) {
                            THROW_ERROR("Tree depth exceeded ITERATIVE_MAX_TREE_DEPTH (", ITERATIVE_MAX_TREE_DEPTH, ")");
                        }
                        stack[stack_top++] = node->childs[c];
                    }
                }
            }
        } else {
            // Use monopole approximation
            const real r_inv = 1.0 / std::sqrt(d2);
            p_i.phi -= m_g_constant * node->mass * r_inv;
            p_i.grav_acc -= d * (m_g_constant * node->mass * pow3(r_inv));
        }
    }
}

#endif // USE_ITERATIVE_TRAVERSAL

// =============================================================================
// Morton Code Particle Reordering
// =============================================================================

#ifdef USE_MORTON_ORDERING

void BHTree::reorder_particles_by_morton(std::vector<SPHParticle> & particles, const int particle_num)
{
    if (particle_num <= 0) return;

    // Compute domain bounds
    real domain_min[DIM];
    real domain_max[DIM];
    real domain_size[DIM];

    for (int d = 0; d < DIM; ++d) {
        if (m_is_periodic) {
            domain_min[d] = m_range_min[d];
            domain_max[d] = m_range_max[d];
        } else {
            domain_min[d] = m_root.center[d] - m_root.edge * 0.5;
            domain_max[d] = m_root.center[d] + m_root.edge * 0.5;
        }
        domain_size[d] = domain_max[d] - domain_min[d];
        if (domain_size[d] <= 0) {
            domain_size[d] = 1.0;  // Avoid division by zero
        }
    }

    // Reorder particles by Morton code
    morton::sort_particles_by_morton(particles, domain_min, domain_size);
}

#endif // USE_MORTON_ORDERING

}