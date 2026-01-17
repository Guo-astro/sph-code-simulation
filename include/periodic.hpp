#pragma once

#include <memory>
#include <cassert>

#include "vector_type.hpp"
#include "parameters.hpp"

namespace sph
{

class Periodic {
    bool m_is_valid;
    bool m_per_dim[DIM];  // Per-dimension periodic flags
    vec_t m_max;
    vec_t m_min;
    vec_t m_range;

public:
    Periodic() : m_is_valid(false) {
        for(int i = 0; i < DIM; ++i) m_per_dim[i] = false;
    }

    // Check if periodic boundaries are enabled
    bool is_valid() const { return m_is_valid; }

    // Check if periodic in specific dimension
    bool is_periodic_in_dim(int dim) const { return m_is_valid && m_per_dim[dim]; }

    void initialize(std::shared_ptr<SPHParameters> param)
    {
        if(param->periodic.is_valid) {
            m_is_valid = true;
            for(int i = 0; i < DIM; ++i) {
                m_max[i] = param->periodic.range_max[i];
                m_min[i] = param->periodic.range_min[i];
                m_per_dim[i] = param->periodic.per_dimension[i];
            }
            m_range = m_max - m_min;
        } else {
            m_is_valid = false;
            for(int i = 0; i < DIM; ++i) m_per_dim[i] = false;
        }
    }

    vec_t calc_r_ij(const vec_t & r_i, const vec_t & r_j) const
    {
        if(m_is_valid) {
            vec_t r_ij;
            for(int i = 0; i < DIM; ++i) {
                real diff = r_i[i] - r_j[i];
                if(m_per_dim[i]) {
                    // Apply minimum image convention only in periodic dimensions
                    real diff_plus = diff + m_range[i];
                    real diff_minus = diff - m_range[i];
                    // Pick the one with smallest absolute value
                    if(std::abs(diff) <= std::abs(diff_plus) && std::abs(diff) <= std::abs(diff_minus)) {
                        r_ij[i] = diff;
                    } else if(std::abs(diff_plus) <= std::abs(diff_minus)) {
                        r_ij[i] = diff_plus;
                    } else {
                        r_ij[i] = diff_minus;
                    }
                } else {
                    // Non-periodic dimension: simple difference
                    r_ij[i] = diff;
                }
            }
            return r_ij;
        } else {
            return r_i - r_j;
        }
    }

    void apply(vec_t & r) const
    {
        if(m_is_valid) {
            for(int i = 0; i < DIM; ++i) {
                // Only wrap in periodic dimensions
                if(m_per_dim[i]) {
                    if(r[i] < m_min[i]) {
                        r[i] += m_range[i];
                    } else if(r[i] > m_max[i]) {
                        r[i] -= m_range[i];
                    }
                }
            }
        }
    }
};

}