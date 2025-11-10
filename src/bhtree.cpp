#include <cassert>
#include <iostream>

#include "parameters.hpp"
#include "bhtree.hpp"
#include "openmp.hpp"
#include "exception.hpp"
#include "periodic.hpp"

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
    }
}

void BHTree::resize(const int particle_num, const int tree_size)
{
    assert(m_nodes.get() == nullptr);

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
    m_root.calc_force(p_i, m_theta2, m_g_constant, m_periodic.get());
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

void BHTree::BHNode::calc_force(SPHParticle & p_i, const real theta2, const real g_constant, const Periodic * periodic)
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
                    p_i.phi -= g_constant * p->mass * (f(r, p_i.sml) + f(r, p->sml)) * 0.5;
                    // Gravity should ADD attractive force
                    // r_ij points FROM j TO i (outward from j's perspective)  
                    // Attractive force points FROM i TO j (opposite of r_ij)
                    // So force = -r_ij * |force|, and we ADD it
                    p_i.acc -= r_ij * (g_constant * p->mass * (g(r, p_i.sml) + g(r, p->sml)) * 0.5);
                }
                p = p->next;
            }
        } else {
            // Internal node: recurse to children
            for(int i = 0; i < NCHILD; ++i) {
                if(childs[i]) {
                    childs[i]->calc_force(p_i, theta2, g_constant, periodic);
                }
            }
        }
    } else {
        // Can use monopole approximation safely
        const real r_inv = 1.0 / std::sqrt(d2);
        p_i.phi -= g_constant * mass * r_inv;
        p_i.acc -= d * (g_constant * mass * pow3(r_inv));
    }
}

}