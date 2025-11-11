#ifndef SPH_OUTPUT_MANAGER_HPP
#define SPH_OUTPUT_MANAGER_HPP

#include "defines.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "output_metadata.hpp"
#include "units.hpp"
#include "writers/csv_writer.hpp"
#include "writers/hdf5_writer.hpp"
#include "writers/vtk_writer.hpp"
#include <string>
#include <memory>
#include <fstream>

namespace sph {

// Forward declarations
class Simulation;
struct SPHParameters;

/**
 * @brief Manages all output operations for SPH simulations
 * 
 * Coordinates multiple output formats (CSV, HDF5, VTK) and handles:
 * - Snapshot output for visualization and resume (SSOT)
 * - Energy logging
 * - Metadata management
 * 
 * Design follows RAII and uses smart pointers for resource management.
 */
class OutputManager {
public:
    /**
     * @brief Constructor
     * @param config Output configuration
     * @param units Unit system for conversions
     * @param output_dir Base output directory
     */
    OutputManager(const OutputConfig& config, 
                  const UnitSystem& units,
                  const std::string& output_dir);
    
    /**
     * @brief Destructor - closes all open files
     */
    ~OutputManager();
    
    // Delete copy constructor and assignment (file handles not copyable)
    OutputManager(const OutputManager&) = delete;
    OutputManager& operator=(const OutputManager&) = delete;
    
    /**
     * @brief Initialize output system (create directories, open energy file)
     * @return true if successful
     */
    bool initialize();
    
    /**
     * @brief Write snapshot output (for visualization and resume)
     * @param sim Simulation state
     * @param params Simulation parameters
     * @param count Snapshot counter
     * @return true if successful
     */
    bool write_snapshot(std::shared_ptr<Simulation> sim,
                       std::shared_ptr<SPHParameters> params,
                       int count);
    
    /**
     * @brief Load snapshot file for resume (SSOT mode)
     * @param filepath Path to snapshot file (.csv, .h5, or .vtk)
     * @param sim Simulation instance to populate with loaded data
     * @param output_meta Output metadata including physics parameters (optional output)
     * @return true if successful, false otherwise
     * 
     * This function loads particle data and metadata from a snapshot file.
     * When resuming, physics parameters from output_meta should override config values (SSOT).
     */
    bool load_for_resume(const std::string& filepath,
                         std::shared_ptr<Simulation> sim,
                         OutputMetadata* output_meta = nullptr);
    
    /**
     * @brief Write energy values to energy.dat file
     * @param time Current simulation time
     * @param kinetic Kinetic energy
     * @param thermal Thermal energy
     * @param potential Gravitational potential energy
     * @return true if successful
     */
    bool write_energy(real time, real kinetic, real thermal, real potential);
    
    /**
     * @brief Get output configuration
     * @return const reference to OutputConfig
     */
    const OutputConfig& get_config() const { return m_config; }
    
    /**
     * @brief Get unit system
     * @return const reference to UnitSystem
     */
    const UnitSystem& get_units() const { return m_units; }
    
private:
    OutputConfig m_config;             ///< Output configuration
    UnitSystem m_units;                ///< Unit system
    std::string m_output_dir;          ///< Base output directory
    
    std::unique_ptr<CSVWriter> m_csv_writer;    ///< CSV writer
    std::unique_ptr<HDF5Writer> m_hdf5_writer;  ///< HDF5 writer
    std::unique_ptr<VTKWriter> m_vtk_writer;    ///< VTK writer
    std::ofstream m_energy_file;                ///< Energy log file
    
    /**
     * @brief Build metadata from simulation state and parameters
     * @param sim Simulation state
     * @param params Simulation parameters
     * @param step Timestep number
     * @return OutputMetadata structure
     */
    OutputMetadata build_metadata(std::shared_ptr<Simulation> sim,
                                  std::shared_ptr<SPHParameters> params,
                                  int step) const;
    
    /**
     * @brief Compute total energies from particle data
     * @param particles Vector of particles
     * @param use_gravity Whether gravity is enabled
     * @param kinetic Output kinetic energy
     * @param thermal Output thermal energy
     * @param potential Output potential energy
     */
    void compute_energies(const std::vector<SPHParticle*>& particles,
                         bool use_gravity,
                         real& kinetic,
                         real& thermal,
                         real& potential) const;
    
    /**
     * @brief Generate filename for snapshot
     * @param prefix Filename prefix ("snapshot")
     * @param count File counter
     * @param extension File extension (".csv", ".h5", or ".vtk")
     * @return Full filepath
     */
    std::string generate_filename(const std::string& prefix,
                                  int count,
                                  const std::string& extension) const;
};

} // namespace sph

#endif // SPH_OUTPUT_MANAGER_HPP
