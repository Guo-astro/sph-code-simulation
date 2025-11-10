#ifndef SPH_HDF5_WRITER_HPP
#define SPH_HDF5_WRITER_HPP

#include "defines.hpp"
#include "particle.hpp"
#include "output_metadata.hpp"
#include <string>
#include <vector>
#include <hdf5.h>

namespace sph {

/**
 * @brief HDF5 output writer with compression and hierarchical structure
 * 
 * Writes SPH particle data to HDF5 format with:
 * - Hierarchical organization (/metadata, /particles, /energy)
 * - Compression (gzip)
 * - Efficient binary storage
 * - Self-describing format
 * - Python/yt compatible
 * 
 * File structure:
 * /metadata/
 *   - simulation_name, timestamp, step, time, particle_count
 *   - unit_system (JSON as attribute)
 *   - parameters (JSON as attribute)
 *   - checkpoint (JSON as attribute, if applicable)
 * /particles/
 *   - id [N]
 *   - pos [N, 3]
 *   - vel [N, 3]
 *   - acc [N, 3]
 *   - mass, dens, pres, ene, sml, sound [N]
 *   - alpha, balsara, gradh, phi, neighbor [N]
 * /energy/
 *   - kinetic, thermal, potential, total
 */
class HDF5Writer {
public:
    /**
     * @brief Constructor
     * @param compression Compression level (0-9, 0=none, 9=max, default=6)
     */
    explicit HDF5Writer(int compression = 6);
    
    /**
     * @brief Destructor - closes file if open
     */
    ~HDF5Writer();
    
    // Delete copy constructor and assignment (HDF5 handles not copyable)
    HDF5Writer(const HDF5Writer&) = delete;
    HDF5Writer& operator=(const HDF5Writer&) = delete;
    
    /**
     * @brief Open HDF5 file for writing
     * @param filepath Path to output HDF5 file
     * @param metadata Metadata to write
     * @return true if successful, false otherwise
     */
    bool open(const std::string& filepath, const OutputMetadata& metadata);
    
    /**
     * @brief Write particle data to HDF5 file
     * @param particles Vector of SPH particles
     * @return true if successful, false otherwise
     */
    bool write_particles(const std::vector<SPHParticle*>& particles);
    
    /**
     * @brief Close the HDF5 file
     */
    void close();
    
    /**
     * @brief Check if file is currently open
     * @return true if open, false otherwise
     */
    bool is_open() const { return m_file_id >= 0; }
    
    /**
     * @brief Read particles from HDF5 file (for resume)
     * @param filepath Path to HDF5 file
     * @param particles Output vector of particles (will be cleared and filled)
     * @return true if successful, false otherwise
     */
    static bool read_particles(const std::string& filepath, std::vector<SPHParticle*>& particles);
    
    /**
     * @brief Read metadata from HDF5 file (for resume)
     * @param filepath Path to HDF5 file
     * @param metadata Output metadata
     * @return true if successful, false otherwise
     */
    static bool read_metadata(const std::string& filepath, OutputMetadata& metadata);
    
    /**
     * @brief Get file extension for HDF5 format
     * @return ".h5"
     */
    static std::string get_extension() { return ".h5"; }
    
private:
    hid_t m_file_id;           ///< HDF5 file handle
    int m_compression;         ///< Compression level (0-9)
    
    /**
     * @brief Create and write metadata group
     * @param metadata Simulation metadata
     * @return true if successful
     */
    bool write_metadata_group(const OutputMetadata& metadata);
    
    /**
     * @brief Create and write particles group
     * @param particles SPH particles
     * @return true if successful
     */
    bool write_particles_group(const std::vector<SPHParticle*>& particles);
    
    /**
     * @brief Create and write energy group
     * @param metadata Metadata containing energy values
     * @return true if successful
     */
    bool write_energy_group(const OutputMetadata& metadata);
    
    /**
     * @brief Write string attribute to HDF5 group
     * @param loc_id Location ID (file or group)
     * @param name Attribute name
     * @param value String value
     * @return true if successful
     */
    static bool write_string_attribute(hid_t loc_id, const char* name, const std::string& value);
    
    /**
     * @brief Read string attribute from HDF5 group
     * @param loc_id Location ID (file or group)
     * @param name Attribute name
     * @param value Output string value
     * @return true if successful
     */
    static bool read_string_attribute(hid_t loc_id, const char* name, std::string& value);
};

} // namespace sph

#endif // SPH_HDF5_WRITER_HPP
