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
 * @brief Checkpoint metadata structure
 * 
 * Stores information about the checkpoint state for resuming simulations
 */
struct CheckpointMetadata {
    int version = 1;                    // Checkpoint format version
    real time = 0.0;                    // Current simulation time
    int step = 0;                       // Current step/loop count
    int particle_num = 0;               // Number of particles
    
    // Relaxation-specific metadata
    bool is_relaxation = false;         // Is this a relaxation checkpoint?
    int relaxation_step = 0;            // Current relaxation iteration
    int relaxation_total_steps = 0;     // Total relaxation steps planned
    real accumulated_time = 0.0;        // Accumulated physical time during relaxation
    
    // Lane-Emden relaxation parameters (needed for resume)
    real alpha_scaling = 1.0;           // Scaling factor for relaxation
    real rho_center = 1.0;              // Central density
    real K = 1.0;                       // Polytropic constant
    real R = 1.0;                       // Star radius
    real M_total = 1.0;                 // Total mass
    
    // Configuration fingerprint (for validation)
    std::string config_hash;            // Hash of configuration file
    std::string preset_name;            // Preset name if using Lane-Emden
};

/**
 * @brief Configuration for checkpoint behavior
 */
struct CheckpointConfig {
    bool enabled = false;                           // Enable checkpointing
    int save_interval = 100;                        // Save every N steps/iterations
    std::string directory = "checkpoints";          // Checkpoint directory
    bool save_on_exit = true;                       // Always save when exiting
    int max_checkpoints = 5;                        // Keep only N most recent checkpoints
    bool compress = false;                          // Compress checkpoint files (future)
    
    // Auto-resume settings
    bool auto_resume = false;                       // Automatically resume from latest checkpoint
    std::string resume_file;                        // Specific checkpoint to resume from
};

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
    
    // === Resume/Checkpoint Data (optional) ===
    bool is_checkpoint;               ///< True if this is a checkpoint file
    CheckpointMetadata checkpoint_data; ///< Checkpoint-specific metadata (if is_checkpoint=true)
    
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
 */
struct OutputConfig {
    // Output formats to enable
    bool enable_csv;                  ///< Write CSV files
    bool enable_hdf5;                 ///< Write HDF5 files
    bool enable_energy_file;          ///< Write energy.dat file
    
    // CSV options
    int csv_precision;                ///< Floating point precision for CSV
    
    // HDF5 options
    int hdf5_compression;             ///< Compression level (0-9, 0=none, 9=max)
    bool hdf5_single_file_series;     ///< Store time series in single HDF5 file (future)
    
    // Unit system for output
    UnitSystem::Type output_unit_type; ///< Which unit system to use for output
    
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
