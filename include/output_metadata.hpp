#ifndef SPH_OUTPUT_METADATA_HPP
#define SPH_OUTPUT_METADATA_HPP

#include "defines.hpp"
#include "units.hpp"
#include "parameters.hpp"
#include <nlohmann/json.hpp>
#include <string>
#include <ctime>

namespace sph {

/**
 * @brief Comprehensive metadata for SPH simulation outputs
 * 
 * Contains all information needed to understand and reproduce a simulation snapshot,
 * including provenance, simulation state, physics parameters, and unit definitions.
 */
struct OutputMetadata {
    // === Provenance ===
    std::string simulation_name;      ///< Descriptive name of simulation
    std::string timestamp;            ///< ISO 8601 timestamp of output generation
    int format_version;               ///< Output format version number
    
    // === Simulation State ===
    int step;                         ///< Timestep number
    real time_code;                   ///< Simulation time in code units
    real time_physical;               ///< Simulation time in physical units
    int particle_count;               ///< Number of particles
    
    // === Physics Parameters ===
    SPHType sph_type;                 ///< SPH formulation (SSPH, DISPH, GSPH)
    KernelType kernel_type;           ///< Kernel function type
    real gamma;                       ///< Adiabatic index
    real gravitational_constant;      ///< G in code units
    int neighbor_number;              ///< Target neighbor count
    bool use_gravity;                 ///< Whether gravity is enabled
    bool use_balsara;                 ///< Balsara switch enabled
    bool use_time_dependent_av;       ///< Time-dependent AV enabled
    
    // === Unit System ===
    UnitSystem units;                 ///< Unit conversion system
    
    // === Relaxation-Specific Metadata (for Lane-Emden resume) ===
    bool is_relaxation = false;         ///< Is this from a relaxation phase?
    int relaxation_step = 0;            ///< Current relaxation iteration
    int relaxation_total_steps = 0;     ///< Total relaxation steps planned
    real accumulated_time = 0.0;        ///< Accumulated physical time during relaxation
    real alpha_scaling = 1.0;           ///< Scaling factor for relaxation
    real rho_center = 1.0;              ///< Central density
    real K = 1.0;                       ///< Polytropic constant
    real R = 1.0;                       ///< Star radius
    real M_total = 1.0;                 ///< Total mass
    
    // === Energy Diagnostics ===
    real kinetic_energy;              ///< Total kinetic energy (code units)
    real thermal_energy;              ///< Total thermal energy (code units)
    real potential_energy;            ///< Total gravitational potential energy (code units)
    real total_energy;                ///< Total energy (code units)
    
    /**
     * @brief Default constructor
     */
    OutputMetadata();
    
    /**
     * @brief Generate ISO 8601 timestamp for current time
     * @return Timestamp string (e.g., "2025-11-11T14:30:00Z")
     */
    static std::string generate_timestamp();
    
    /**
     * @brief Serialize metadata to JSON
     * @return JSON object with all metadata fields
     */
    nlohmann::json to_json() const;
    
    /**
     * @brief Deserialize metadata from JSON
     * @param j JSON object
     * @return OutputMetadata instance
     */
    static OutputMetadata from_json(const nlohmann::json& j);
    
    /**
     * @brief Get SPH type name as string
     * @return String representation ("Standard SPH", "DISPH", "Godunov SPH")
     */
    std::string get_sph_type_name() const;
    
    /**
     * @brief Get kernel type name as string
     * @return String representation ("Cubic Spline", "Wendland C4")
     */
    std::string get_kernel_type_name() const;
};

/**
 * @brief Configuration for output system
 * 
 * Supports two JSON formats:
 * 1. Array of format objects: [{"type": "csv", "precision": 16}, {"type": "hdf5", "compression": 6}]
 * 2. Simple array (legacy): ["csv", "hdf5"] with separate option fields
 */
struct OutputConfig {
    // Output formats - specified as array in JSON
    std::vector<std::string> formats;  ///< Output formats to enable
    bool enable_energy_file;           ///< Write energy.dat file
    
    // CSV options
    int csv_precision;                ///< Floating point precision for CSV
    
    // HDF5 options
    int hdf5_compression;             ///< Compression level (0-9, 0=none, 9=max)
    bool hdf5_single_file_series;     ///< Store time series in single HDF5 file (future)
    
    // VTK options
    bool vtk_binary;                  ///< Use binary VTK format (vs ASCII)
    
    // Unit system for output
    UnitSystem::Type output_unit_type; ///< Which unit system to use for output
    
    /**
     * @brief Check if a format is enabled
     * @param format Format name ("csv", "hdf5", "vtk")
     * @return true if format is enabled
     */
    bool is_format_enabled(const std::string& format) const;
    
    /**
     * @brief Default constructor with sensible defaults
     */
    OutputConfig();
    
    /**
     * @brief Load configuration from JSON
     * @param j JSON configuration object
     * @return OutputConfig instance
     */
    static OutputConfig from_json(const nlohmann::json& j);
    
    /**
     * @brief Serialize configuration to JSON
     * @return JSON object
     */
    nlohmann::json to_json() const;
};

} // namespace sph

#endif // SPH_OUTPUT_METADATA_HPP
