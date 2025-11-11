#include "writers/vtk_writer.hpp"
#include "logger.hpp"
#include <iomanip>
#include <cstring>

namespace sph {

VTKWriter::VTKWriter(bool binary)
    : m_binary(binary)
{
}

VTKWriter::~VTKWriter() {
    close();
}

bool VTKWriter::open(const std::string& filepath, const OutputMetadata& metadata) {
    if (m_file.is_open()) {
        close();
    }
    
    m_metadata = metadata;
    
    // Open file in binary mode if requested
    std::ios_base::openmode mode = std::ios::out;
    if (m_binary) {
        mode |= std::ios::binary;
    }
    
    m_file.open(filepath, mode);
    if (!m_file.is_open()) {
        WRITE_LOG << "Failed to open VTK file: " << filepath;
        return false;
    }
    
    write_header();
    return true;
}

void VTKWriter::write_header() {
    // VTK Legacy format header
    m_file << "# vtk DataFile Version 3.0\n";
    m_file << "SPH Simulation - " << m_metadata.simulation_name << "\n";
    
    if (m_binary) {
        m_file << "BINARY\n";
    } else {
        m_file << "ASCII\n";
    }
    
    m_file << "DATASET UNSTRUCTURED_GRID\n";
}

bool VTKWriter::write_particles(const std::vector<SPHParticle*>& particles) {
    if (!m_file.is_open()) {
        WRITE_LOG << "VTK file is not open";
        return false;
    }
    
    if (particles.empty()) {
        WRITE_LOG << "No particles to write";
        return false;
    }
    
    const int n = static_cast<int>(particles.size());
    
    // Write points (positions)
    write_points(particles);
    
    // Write cells (each particle is a vertex)
    m_file << "CELLS " << n << " " << (n * 2) << "\n";
    for (int i = 0; i < n; ++i) {
        if (m_binary) {
            write_binary_int(1);  // Number of points in cell
            write_binary_int(i);  // Point ID
        } else {
            m_file << "1 " << i << "\n";
        }
    }
    
    // Write cell types (1 = VTK_VERTEX)
    m_file << "CELL_TYPES " << n << "\n";
    for (int i = 0; i < n; ++i) {
        if (m_binary) {
            write_binary_int(1);
        } else {
            m_file << "1\n";
        }
    }
    
    // Write point data (scalar and vector fields)
    write_point_data(particles);
    
    return true;
}

void VTKWriter::write_points(const std::vector<SPHParticle*>& particles) {
    const int n = static_cast<int>(particles.size());
    
    m_file << "POINTS " << n << " float\n";
    
    for (const auto* p : particles) {
        if (m_binary) {
            write_binary_float(static_cast<float>(p->pos[0]));
#if DIM >= 2
            write_binary_float(static_cast<float>(p->pos[1]));
#else
            write_binary_float(0.0f);
#endif
#if DIM >= 3
            write_binary_float(static_cast<float>(p->pos[2]));
#else
            write_binary_float(0.0f);
#endif
        } else {
            m_file << p->pos[0];
#if DIM >= 2
            m_file << " " << p->pos[1];
#else
            m_file << " 0.0";
#endif
#if DIM >= 3
            m_file << " " << p->pos[2];
#else
            m_file << " 0.0";
#endif
            m_file << "\n";
        }
    }
    
    if (!m_binary) {
        m_file << "\n";
    }
}

void VTKWriter::write_point_data(const std::vector<SPHParticle*>& particles) {
    const int n = static_cast<int>(particles.size());
    
    m_file << "POINT_DATA " << n << "\n";
    
    // Density
    m_file << "SCALARS density float 1\n";
    m_file << "LOOKUP_TABLE default\n";
    for (const auto* p : particles) {
        if (m_binary) {
            write_binary_float(static_cast<float>(p->dens));
        } else {
            m_file << p->dens << "\n";
        }
    }
    if (!m_binary) m_file << "\n";
    
    // Pressure
    m_file << "SCALARS pressure float 1\n";
    m_file << "LOOKUP_TABLE default\n";
    for (const auto* p : particles) {
        if (m_binary) {
            write_binary_float(static_cast<float>(p->pres));
        } else {
            m_file << p->pres << "\n";
        }
    }
    if (!m_binary) m_file << "\n";
    
    // Mass
    m_file << "SCALARS mass float 1\n";
    m_file << "LOOKUP_TABLE default\n";
    for (const auto* p : particles) {
        if (m_binary) {
            write_binary_float(static_cast<float>(p->mass));
        } else {
            m_file << p->mass << "\n";
        }
    }
    if (!m_binary) m_file << "\n";
    
    // Internal energy
    m_file << "SCALARS internal_energy float 1\n";
    m_file << "LOOKUP_TABLE default\n";
    for (const auto* p : particles) {
        if (m_binary) {
            write_binary_float(static_cast<float>(p->ene));
        } else {
            m_file << p->ene << "\n";
        }
    }
    if (!m_binary) m_file << "\n";
    
    // Smoothing length
    m_file << "SCALARS smoothing_length float 1\n";
    m_file << "LOOKUP_TABLE default\n";
    for (const auto* p : particles) {
        if (m_binary) {
            write_binary_float(static_cast<float>(p->sml));
        } else {
            m_file << p->sml << "\n";
        }
    }
    if (!m_binary) m_file << "\n";
    
    // Sound speed
    m_file << "SCALARS sound_speed float 1\n";
    m_file << "LOOKUP_TABLE default\n";
    for (const auto* p : particles) {
        if (m_binary) {
            write_binary_float(static_cast<float>(p->sound));
        } else {
            m_file << p->sound << "\n";
        }
    }
    if (!m_binary) m_file << "\n";
    
    // Particle ID
    m_file << "SCALARS particle_id int 1\n";
    m_file << "LOOKUP_TABLE default\n";
    for (const auto* p : particles) {
        if (m_binary) {
            write_binary_int(p->id);
        } else {
            m_file << p->id << "\n";
        }
    }
    if (!m_binary) m_file << "\n";
    
    // Write vector fields
    write_vectors(particles);
}

void VTKWriter::write_vectors(const std::vector<SPHParticle*>& particles) {
    const int n = static_cast<int>(particles.size());
    
    // Velocity
    m_file << "VECTORS velocity float\n";
    for (const auto* p : particles) {
        if (m_binary) {
            write_binary_float(static_cast<float>(p->vel[0]));
#if DIM >= 2
            write_binary_float(static_cast<float>(p->vel[1]));
#else
            write_binary_float(0.0f);
#endif
#if DIM >= 3
            write_binary_float(static_cast<float>(p->vel[2]));
#else
            write_binary_float(0.0f);
#endif
        } else {
            m_file << p->vel[0];
#if DIM >= 2
            m_file << " " << p->vel[1];
#else
            m_file << " 0.0";
#endif
#if DIM >= 3
            m_file << " " << p->vel[2];
#else
            m_file << " 0.0";
#endif
            m_file << "\n";
        }
    }
    if (!m_binary) m_file << "\n";
    
    // Acceleration
    m_file << "VECTORS acceleration float\n";
    for (const auto* p : particles) {
        if (m_binary) {
            write_binary_float(static_cast<float>(p->acc[0]));
#if DIM >= 2
            write_binary_float(static_cast<float>(p->acc[1]));
#else
            write_binary_float(0.0f);
#endif
#if DIM >= 3
            write_binary_float(static_cast<float>(p->acc[2]));
#else
            write_binary_float(0.0f);
#endif
        } else {
            m_file << p->acc[0];
#if DIM >= 2
            m_file << " " << p->acc[1];
#else
            m_file << " 0.0";
#endif
#if DIM >= 3
            m_file << " " << p->acc[2];
#else
            m_file << " 0.0";
#endif
            m_file << "\n";
        }
    }
    if (!m_binary) m_file << "\n";
}

void VTKWriter::close() {
    if (m_file.is_open()) {
        m_file.close();
    }
}

void VTKWriter::write_binary_float(float value) {
    // VTK binary format uses big-endian
    float swapped = swap_endian(value);
    m_file.write(reinterpret_cast<const char*>(&swapped), sizeof(float));
}

void VTKWriter::write_binary_int(int value) {
    // VTK binary format uses big-endian
    int swapped = swap_endian(value);
    m_file.write(reinterpret_cast<const char*>(&swapped), sizeof(int));
}

template<typename T>
T VTKWriter::swap_endian(T value) {
    union {
        T value;
        unsigned char bytes[sizeof(T)];
    } source, dest;
    
    source.value = value;
    
    for (size_t i = 0; i < sizeof(T); ++i) {
        dest.bytes[i] = source.bytes[sizeof(T) - 1 - i];
    }
    
    return dest.value;
}

} // namespace sph
