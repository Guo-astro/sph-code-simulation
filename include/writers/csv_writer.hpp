#ifndef SPH_CSV_WRITER_HPP
#define SPH_CSV_WRITER_HPP

#include "defines.hpp"
#include "particle.hpp"
#include "output_metadata.hpp"
#include <string>
#include <fstream>
#include <vector>
#include <memory>

namespace sph {

/**
 * @brief CSV output writer with self-documenting headers
 * 
 * Writes SPH particle data to CSV format with comprehensive metadata headers.
 * Format is human-readable and can be easily loaded by pandas, Excel, etc.
 * 
 * Header includes:
 * - Simulation metadata (name, timestamp, version)
 * - Unit system information
 * - Physics parameters
 * - Simulation state (step, time, particle count)
 * - Energy diagnostics
 * - Column definitions
 */
class CSVWriter {
public:
    /**
     * @brief Constructor
     * @param precision Number of decimal digits for floating point output (default: 16)
     */
    explicit CSVWriter(int precision = 16);
    
    /**
     * @brief Destructor - closes file if open
     */
    ~CSVWriter();
    
    // Delete copy constructor and assignment (file handle not copyable)
    CSVWriter(const CSVWriter&) = delete;
    CSVWriter& operator=(const CSVWriter&) = delete;
    
    /**
     * @brief Open CSV file for writing
     * @param filepath Path to output CSV file
     * @param metadata Metadata to write in header
     * @return true if successful, false otherwise
     */
    bool open(const std::string& filepath, const OutputMetadata& metadata);
    
    /**
     * @brief Write particle data to CSV file
     * @param particles Vector of SPH particles
     * @return true if successful, false otherwise
     */
    bool write_particles(const std::vector<SPHParticle*>& particles);
    
    /**
     * @brief Close the CSV file
     */
    void close();
    
    /**
     * @brief Check if file is currently open
     * @return true if open, false otherwise
     */
    bool is_open() const { return m_file.is_open(); }
    
    /**
     * @brief Get file extension for CSV format
     * @return ".csv"
     */
    static std::string get_extension() { return ".csv"; }
    
    /**
     * @brief Read particles from CSV checkpoint file
     * @param filepath Path to CSV checkpoint file
     * @param particles Vector to store loaded particles (as pointers, caller owns memory)
     * @return true if successful, false otherwise
     */
    static bool read_particles(const std::string& filepath, std::vector<SPHParticle*>& particles);
    
    /**
     * @brief Read metadata from CSV checkpoint file
     * @param filepath Path to CSV checkpoint file
     * @param metadata Metadata structure to populate
     * @return true if successful, false otherwise
     */
    static bool read_metadata(const std::string& filepath, OutputMetadata& metadata);
    
private:
    std::ofstream m_file;      ///< Output file stream
    int m_precision;           ///< Floating point precision
    
    /**
     * @brief Write comprehensive header with metadata
     * @param metadata Simulation metadata
     */
    void write_header(const OutputMetadata& metadata);
    
    /**
     * @brief Write column names
     */
    void write_column_names();
};

} // namespace sph

#endif // SPH_CSV_WRITER_HPP
