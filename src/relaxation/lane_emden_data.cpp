#include "relaxation/lane_emden_data.hpp"
#include "exception.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <algorithm>  // For std::lower_bound

namespace sph
{

LaneEmdenData::LaneEmdenData()
    : m_xi_1(0.0)
    , m_dtheta_1(0.0)
    , m_loaded(false)
{
}

void LaneEmdenData::load_from_file(const std::string& filename)
{
    std::ifstream infile(filename);
    
    if(!infile.is_open()) {
        THROW_ERROR("Cannot open Lane-Emden solution file: " + filename);
    }
    
    // Clear existing data
    m_xi_array.clear();
    m_theta_array.clear();
    m_dtheta_array.clear();
    
    std::string line;
    while(std::getline(infile, line)) {
        // Parse header comments
        if(line[0] == '#') {
            if(line.find("xi_1") != std::string::npos) {
                sscanf(line.c_str(), "# xi_1 = %lf", &m_xi_1);
            } else if(line.find("dtheta_1") != std::string::npos) {
                sscanf(line.c_str(), "# dtheta_1 = %lf", &m_dtheta_1);
            }
            continue;
        }
        
        // Parse data: xi theta dtheta/dxi
        real xi, theta, dtheta;
        std::istringstream iss(line);
        if(iss >> xi >> theta >> dtheta) {
            m_xi_array.push_back(xi);
            m_theta_array.push_back(theta);
            m_dtheta_array.push_back(dtheta);
        }
    }
    infile.close();
    
    if(m_xi_array.empty()) {
        THROW_ERROR("Failed to read Lane-Emden solution data from: " + filename);
    }
    
    m_loaded = true;
    
    std::cout << "LaneEmdenData: Loaded " << m_xi_array.size() << " points from " << filename << std::endl;
    std::cout << "LaneEmdenData: ξ₁ = " << m_xi_1 << ", |dθ/dξ|_{ξ₁} = " << std::abs(m_dtheta_1) << std::endl;
}

real LaneEmdenData::get_theta(real xi) const
{
    if(!m_loaded) {
        THROW_ERROR("LaneEmdenData: Data not loaded, call load_from_file() first");
    }
    
    if(xi >= m_xi_1) return 0.0;
    if(xi <= 0.0) return 1.0;
    
    // Binary search for interpolation interval (O(log n) vs O(n))
    auto it = std::lower_bound(m_xi_array.begin(), m_xi_array.end(), xi);
    
    // Handle edge cases
    if(it == m_xi_array.end()) return 0.0;  // Beyond data range
    if(it == m_xi_array.begin()) return m_theta_array[0];  // At or before first point
    
    // Get indices for interpolation
    size_t i1 = std::distance(m_xi_array.begin(), it);
    size_t i0 = i1 - 1;
    
    // Linear interpolation between points
    const real frac = (xi - m_xi_array[i0]) / (m_xi_array[i1] - m_xi_array[i0]);
    return m_theta_array[i0] * (1.0 - frac) + m_theta_array[i1] * frac;
}

real LaneEmdenData::dtheta_dxi(real xi) const
{
    if(!m_loaded) {
        THROW_ERROR("LaneEmdenData: Data not loaded, call load_from_file() first");
    }
    
    if(xi >= m_xi_1) return 0.0;
    if(xi <= 0.0) return 0.0;  // dθ/dξ = 0 at center
    
    // Binary search for interpolation interval (O(log n) vs O(n))
    auto it = std::lower_bound(m_xi_array.begin(), m_xi_array.end(), xi);
    
    // Handle edge cases
    if(it == m_xi_array.end()) return 0.0;  // Beyond data range
    if(it == m_xi_array.begin()) return m_dtheta_array[0];  // At or before first point
    
    // Get indices for interpolation
    size_t i1 = std::distance(m_xi_array.begin(), it);
    size_t i0 = i1 - 1;
    
    // Linear interpolation between points
    const real frac = (xi - m_xi_array[i0]) / (m_xi_array[i1] - m_xi_array[i0]);
    return m_dtheta_array[i0] * (1.0 - frac) + m_dtheta_array[i1] * frac;
}

} // namespace sph
