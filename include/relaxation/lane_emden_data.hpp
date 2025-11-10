#pragma once

#include <vector>
#include <string>
#include "defines.hpp"

namespace sph
{

/**
 * @brief Loads and interpolates Lane-Emden equation solution data
 * 
 * Reads the numerical solution from data/lane_emden/n1.5_3d.dat
 * and provides interpolation methods for θ(ξ) and dθ/dξ
 */
class LaneEmdenData
{
public:
    LaneEmdenData();
    
    /**
     * @brief Load Lane-Emden solution from file
     * @param filename Path to data file (default: data/lane_emden/n1.5_3d.dat)
     */
    void load_from_file(const std::string& filename = "data/lane_emden/n1.5_3d.dat");
    
    /**
     * @brief Get θ(ξ) by linear interpolation
     * @param xi Dimensionless radius
     * @return theta value at xi
     */
    real get_theta(real xi) const;
    
    /**
     * @brief Get dθ/dξ by linear interpolation
     * @param xi Dimensionless radius
     * @return derivative of theta at xi
     */
    real dtheta_dxi(real xi) const;
    
    /**
     * @brief Get first zero of θ (surface radius in dimensionless units)
     */
    real get_xi_1() const { return m_xi_1; }
    
    /**
     * @brief Get |dθ/dξ| at ξ₁ (used for mass calculation)
     */
    real get_dtheta_1() const { return m_dtheta_1; }
    
    /**
     * @brief Get number of data points loaded
     */
    int get_num_points() const { return m_xi_array.size(); }
    
private:
    std::vector<real> m_xi_array;      // ξ values
    std::vector<real> m_theta_array;   // θ(ξ) values
    std::vector<real> m_dtheta_array;  // dθ/dξ values
    
    real m_xi_1;       // First zero of θ
    real m_dtheta_1;   // dθ/dξ at ξ₁
    
    bool m_loaded;
};

} // namespace sph
