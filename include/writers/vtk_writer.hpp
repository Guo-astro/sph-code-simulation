#ifndef SPH_VTK_WRITER_HPP
#define SPH_VTK_WRITER_HPP

#include "defines.hpp"
#include "particle.hpp"
#include "output_metadata.hpp"
#include <string>
#include <fstream>
#include <vector>

namespace sph {

/**
 * @brief Writer for VTK (Visualization Toolkit) format files
 * 
 * Writes SPH particle data in legacy VTK format for visualization
 * in ParaView, VisIt, and other VTK-compatible tools.
 * 
 * Format: VTK Legacy (UNSTRUCTURED_GRID with POINTS)
 * https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 */
class VTKWriter {
public:
    /**
     * @brief Constructor
     * @param binary If true, write binary VTK (faster, smaller), else ASCII
     */
    explicit VTKWriter(bool binary = true);
    
    /**
     * @brief Destructor
     */
    ~VTKWriter();
    
    /**
     * @brief Open VTK file for writing
     * @param filepath Path to output file
     * @param metadata Simulation metadata
     * @return true if successful
     */
    bool open(const std::string& filepath, const OutputMetadata& metadata);
    
    /**
     * @brief Write particle data to VTK file
     * @param particles Vector of particle pointers
     * @return true if successful
     */
    bool write_particles(const std::vector<SPHParticle*>& particles);
    
    /**
     * @brief Close the VTK file
     */
    void close();
    
    /**
     * @brief Check if file is currently open
     * @return true if file is open
     */
    bool is_open() const { return m_file.is_open(); }

private:
    std::ofstream m_file;              ///< Output file stream
    bool m_binary;                     ///< Binary mode flag
    OutputMetadata m_metadata;         ///< Cached metadata
    
    /**
     * @brief Write VTK header
     */
    void write_header();
    
    /**
     * @brief Write particle positions as VTK POINTS
     */
    void write_points(const std::vector<SPHParticle*>& particles);
    
    /**
     * @brief Write scalar fields as VTK POINT_DATA
     */
    void write_point_data(const std::vector<SPHParticle*>& particles);
    
    /**
     * @brief Write vector fields as VTK vectors
     */
    void write_vectors(const std::vector<SPHParticle*>& particles);
    
    /**
     * @brief Write a single float in binary format (big-endian)
     */
    void write_binary_float(float value);
    
    /**
     * @brief Write a single int in binary format (big-endian)
     */
    void write_binary_int(int value);
    
    /**
     * @brief Swap bytes for big-endian format
     */
    template<typename T>
    T swap_endian(T value);
};

} // namespace sph

#endif // SPH_VTK_WRITER_HPP
