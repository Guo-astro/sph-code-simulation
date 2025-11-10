# Output & Unit System Implementation Plan

## Overview
This document provides step-by-step actionable tasks for implementing the new output and unit system. Each task is designed to be completable independently with clear acceptance criteria.

---

## Phase 1: Foundation & Dependencies (Day 1)

### Task 1.1: Setup Dependencies
**Estimated Time**: 1 hour

**Steps**:
1. Update `CMakeLists.txt` to find HDF5:
   ```cmake
   find_package(HDF5 REQUIRED COMPONENTS C)
   include_directories(${HDF5_INCLUDE_DIRS})
   ```

2. Add nlohmann/json (replace Boost PropertyTree gradually):
   ```cmake
   include(FetchContent)
   FetchContent_Declare(json
       URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
   FetchContent_MakeAvailable(json)
   ```

3. Update target linking:
   ```cmake
   target_link_libraries(sph 
       ${HDF5_LIBRARIES}
       nlohmann_json::nlohmann_json)
   ```

4. Test build:
   ```bash
   cd build && cmake .. && make
   ```

**Acceptance Criteria**:
- ✅ CMake finds HDF5 successfully
- ✅ Project compiles with new dependencies
- ✅ Can include `<nlohmann/json.hpp>` in code

**Files Modified**:
- `CMakeLists.txt`
- `src/CMakeLists.txt`

---

### Task 1.2: Create Unit System Header
**Estimated Time**: 2 hours

**Steps**:
1. Create `include/units.hpp`:
   ```cpp
   #pragma once
   #include "defines.hpp"
   #include <string>
   #include <nlohmann/json.hpp>
   
   namespace sph {
   
   class UnitSystem {
   public:
       enum class Type {
           CODE,      // Dimensionless code units
           GALACTIC,  // kpc, M_sun, km/s
           SI,        // m, kg, s
           CGS        // cm, g, s
       };
       
       // Constructors
       UnitSystem();  // Default: CODE units
       UnitSystem(Type type, real code_length = 1.0, 
                  real code_mass = 1.0, real code_velocity = 1.0);
       
       // Factory methods
       static UnitSystem create_galactic(real length_kpc = 1.0,
                                        real mass_msun = 1.0,
                                        real vel_kms = 1.0);
       static UnitSystem create_si();
       static UnitSystem create_cgs();
       static UnitSystem create_code();
       
       // Conversion: code → physical
       real to_physical_length(real code_val) const;
       real to_physical_mass(real code_val) const;
       real to_physical_time(real code_val) const;
       real to_physical_velocity(real code_val) const;
       real to_physical_energy(real code_val) const;
       real to_physical_density(real code_val) const;
       real to_physical_pressure(real code_val) const;
       
       // Conversion: physical → code
       real from_physical_length(real phys_val) const;
       real from_physical_mass(real phys_val) const;
       real from_physical_time(real phys_val) const;
       real from_physical_velocity(real phys_val) const;
       real from_physical_energy(real phys_val) const;
       real from_physical_density(real phys_val) const;
       real from_physical_pressure(real phys_val) const;
       
       // Accessors
       Type get_type() const { return m_type; }
       std::string get_type_name() const;
       std::string get_length_unit_name() const;
       std::string get_mass_unit_name() const;
       std::string get_time_unit_name() const;
       std::string get_velocity_unit_name() const;
       
       real get_length_to_cgs() const { return m_length_to_cgs; }
       real get_mass_to_cgs() const { return m_mass_to_cgs; }
       real get_time_to_cgs() const { return m_time_to_cgs; }
       
       // Serialization
       nlohmann::json to_json() const;
       static UnitSystem from_json(const nlohmann::json& j);
       
   private:
       Type m_type;
       
       // Conversion factors to CGS (base units for conversion)
       real m_length_to_cgs;    // code length → cm
       real m_mass_to_cgs;      // code mass → g
       real m_time_to_cgs;      // code time → s
       
       // Derived conversion factors (computed from above)
       real m_velocity_to_cgs;  // code velocity → cm/s
       real m_energy_to_cgs;    // code energy → erg
       real m_density_to_cgs;   // code density → g/cm³
       real m_pressure_to_cgs;  // code pressure → dyne/cm²
       
       void compute_derived_factors();
   };
   
   // Physical constants in CGS
   namespace constants {
       constexpr real G_cgs = 6.67430e-8;        // cm³ g⁻¹ s⁻²
       constexpr real kpc_to_cm = 3.085677581e21; // cm per kpc
       constexpr real msun_to_g = 1.98841e33;     // g per solar mass
       constexpr real km_to_cm = 1.0e5;           // cm per km
       constexpr real pc_to_cm = 3.085677581e18;  // cm per parsec
   }
   
   } // namespace sph
   ```

2. Create corresponding implementation file `src/units.cpp`

3. Implement all conversion methods using conversion factors

4. Add unit tests in `test/test_units.cpp`

**Acceptance Criteria**:
- ✅ UnitSystem class compiles
- ✅ All conversion methods implemented
- ✅ Unit tests pass for Galactic, SI, CGS conversions
- ✅ JSON serialization/deserialization works

**Files Created**:
- `include/units.hpp`
- `src/units.cpp`
- `test/test_units.cpp`

---

### Task 1.3: Implement Unit System Logic
**Estimated Time**: 3 hours

**Implementation Details**:

```cpp
// src/units.cpp
#include "units.hpp"
#include <cmath>
#include <stdexcept>

namespace sph {

UnitSystem::UnitSystem() : m_type(Type::CODE) {
    m_length_to_cgs = 1.0;
    m_mass_to_cgs = 1.0;
    m_time_to_cgs = 1.0;
    compute_derived_factors();
}

UnitSystem::UnitSystem(Type type, real code_length, real code_mass, real code_velocity)
    : m_type(type) {
    
    switch (type) {
    case Type::CODE:
        m_length_to_cgs = 1.0;
        m_mass_to_cgs = 1.0;
        m_time_to_cgs = 1.0;
        break;
        
    case Type::GALACTIC:
        // code_length is in kpc, code_mass in M_sun, code_velocity in km/s
        m_length_to_cgs = code_length * constants::kpc_to_cm;
        m_mass_to_cgs = code_mass * constants::msun_to_g;
        m_time_to_cgs = m_length_to_cgs / (code_velocity * constants::km_to_cm);
        break;
        
    case Type::SI:
        m_length_to_cgs = code_length * 100.0;  // m to cm
        m_mass_to_cgs = code_mass * 1000.0;      // kg to g
        m_time_to_cgs = code_time;               // s to s
        break;
        
    case Type::CGS:
        m_length_to_cgs = code_length;  // cm to cm
        m_mass_to_cgs = code_mass;      // g to g
        m_time_to_cgs = code_time;      // s to s
        break;
    }
    
    compute_derived_factors();
}

void UnitSystem::compute_derived_factors() {
    m_velocity_to_cgs = m_length_to_cgs / m_time_to_cgs;
    m_energy_to_cgs = m_mass_to_cgs * m_velocity_to_cgs * m_velocity_to_cgs;
    m_density_to_cgs = m_mass_to_cgs / (m_length_to_cgs * m_length_to_cgs * m_length_to_cgs);
    m_pressure_to_cgs = m_energy_to_cgs / (m_length_to_cgs * m_length_to_cgs * m_length_to_cgs);
}

real UnitSystem::to_physical_length(real code_val) const {
    return code_val * m_length_to_cgs;
}

// ... implement all other conversion methods similarly

nlohmann::json UnitSystem::to_json() const {
    nlohmann::json j;
    j["type"] = static_cast<int>(m_type);
    j["type_name"] = get_type_name();
    j["length_unit"] = get_length_unit_name();
    j["mass_unit"] = get_mass_unit_name();
    j["time_unit"] = get_time_unit_name();
    j["velocity_unit"] = get_velocity_unit_name();
    j["length_to_cgs"] = m_length_to_cgs;
    j["mass_to_cgs"] = m_mass_to_cgs;
    j["time_to_cgs"] = m_time_to_cgs;
    j["velocity_to_cgs"] = m_velocity_to_cgs;
    j["energy_to_cgs"] = m_energy_to_cgs;
    j["density_to_cgs"] = m_density_to_cgs;
    j["pressure_to_cgs"] = m_pressure_to_cgs;
    return j;
}

} // namespace sph
```

**Test Cases**:
```cpp
TEST(UnitSystem, GalacticLengthConversion) {
    auto units = UnitSystem::create_galactic(1.0, 1.0, 1.0);
    EXPECT_DOUBLE_EQ(units.to_physical_length(1.0), 3.085677581e21);
}

TEST(UnitSystem, VelocityConsistency) {
    auto units = UnitSystem::create_galactic(1.0, 1.0, 1.0);
    real vel_code = 1.0;  // 1 km/s in code units
    real vel_cgs = units.to_physical_velocity(vel_code);
    EXPECT_DOUBLE_EQ(vel_cgs, 1.0e5);  // 1 km/s = 10^5 cm/s
}
```

**Acceptance Criteria**:
- ✅ All unit conversions mathematically correct
- ✅ Galactic units: 1 kpc = 3.086e21 cm verified
- ✅ SI/CGS conversions accurate
- ✅ Round-trip conversions work (code→phys→code)

---

## Phase 2: CSV Output Format (Day 2)

### Task 2.1: Create Output Metadata Structure
**Estimated Time**: 1 hour

**Steps**:
1. Create `include/output_metadata.hpp`:
   ```cpp
   #pragma once
   #include "defines.hpp"
   #include "parameters.hpp"
   #include "units.hpp"
   #include "checkpoint.hpp"
   #include <string>
   #include <nlohmann/json.hpp>
   
   namespace sph {
   
   struct OutputMetadata {
       // Provenance
       std::string simulation_name;
       std::string timestamp;
       int format_version = 1;
       
       // State
       int step;
       real time_code;
       real time_physical;
       int particle_count;
       
       // Physics (will be populated from SPHParameters)
       real gamma;
       real G;
       int neighbor_number;
       std::string sph_type;
       std::string kernel_type;
       
       // Units
       UnitSystem units;
       
       // Resume data (optional)
       bool is_checkpoint = false;
       CheckpointMetadata checkpoint_data;
       
       // Energy
       real kinetic_energy = 0.0;
       real thermal_energy = 0.0;
       real potential_energy = 0.0;
       real total_energy = 0.0;
       
       // Methods
       nlohmann::json to_json() const;
       static OutputMetadata from_json(const nlohmann::json& j);
       
       std::string generate_timestamp() const;
   };
   
   } // namespace sph
   ```

2. Implement in `src/output_metadata.cpp`

**Acceptance Criteria**:
- ✅ OutputMetadata structure defined
- ✅ JSON serialization implemented
- ✅ Timestamp generation works (ISO 8601 format)

**Files Created**:
- `include/output_metadata.hpp`
- `src/output_metadata.cpp`

---

### Task 2.2: Create CSV Writer Class
**Estimated Time**: 3 hours

**Steps**:
1. Create `include/writers/csv_writer.hpp`:
   ```cpp
   #pragma once
   #include "output_metadata.hpp"
   #include "particle.hpp"
   #include <fstream>
   #include <vector>
   #include <iomanip>
   
   namespace sph {
   
   class CSVWriter {
   public:
       struct Config {
           int precision = 16;
           bool include_header = true;
           char delimiter = ',';
       };
       
       CSVWriter(const Config& config = Config());
       ~CSVWriter();
       
       bool open(const std::string& filepath, const OutputMetadata& metadata);
       bool write_particles(const std::vector<SPHParticle>& particles);
       bool close();
       
       std::string get_extension() const { return ".csv"; }
       
   private:
       Config m_config;
       std::ofstream m_file;
       bool m_is_open = false;
       
       void write_header(const OutputMetadata& metadata);
       void write_column_names();
   };
   
   } // namespace sph
   ```

2. Implement in `src/writers/csv_writer.cpp`:
   - Header with metadata comments
   - Column names
   - Particle data with high precision

3. Add directory `src/writers/` to `src/CMakeLists.txt`

**Header Format Example**:
```
# SPH Simulation Output - CSV Format v1.0
# Timestamp: 2025-11-10T21:30:00Z
# Simulation: Lane-Emden n=1.5 polytrope
#
# === Unit System ===
# Type: Galactic
# Length: kpc (3.086e+21 cm)
# Mass: M_sun (1.989e+33 g)
# Time: Myr (3.156e+13 s)
# Velocity: km/s (1.000e+05 cm/s)
#
# === Simulation State ===
# Step: 95000
# Time (code): 112.045
# Time (physical): 1.234 Myr
# Particle Count: 5400
#
# === Physics Parameters ===
# Gamma: 1.667
# G: 1.0
# Neighbor Number: 50
# SPH Type: Standard SPH
# Kernel: Wendland
#
# === Energy ===
# Kinetic: 123.45 [code units]
# Thermal: 234.56 [code units]
# Potential: -345.67 [code units]
# Total: 12.34 [code units]
#
# === Columns ===
id,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,acc_x,acc_y,acc_z,mass,dens,pres,ene,sml,sound,alpha,balsara,gradh,phi,neighbor
```

**Acceptance Criteria**:
- ✅ CSV file created with complete header
- ✅ All 20 particle fields written
- ✅ High precision output (16 digits)
- ✅ Valid CSV format (can load in pandas)

**Files Created**:
- `include/writers/csv_writer.hpp`
- `src/writers/csv_writer.cpp`

---

### Task 2.3: Test CSV Writer
**Estimated Time**: 1 hour

**Steps**:
1. Create test file `test/test_csv_writer.cpp`
2. Write test particles
3. Read back and verify
4. Test with Python pandas

**Test Code**:
```cpp
TEST(CSVWriter, WriteAndVerify) {
    // Create test particles
    std::vector<SPHParticle> particles(100);
    for (int i = 0; i < 100; ++i) {
        particles[i].id = i;
        particles[i].pos[0] = i * 0.1;
        particles[i].mass = 1.0;
        // ... fill other fields
    }
    
    // Create metadata
    OutputMetadata meta;
    meta.simulation_name = "Test";
    meta.step = 0;
    meta.particle_count = 100;
    meta.units = UnitSystem::create_galactic();
    
    // Write
    CSVWriter writer;
    ASSERT_TRUE(writer.open("test_output.csv", meta));
    ASSERT_TRUE(writer.write_particles(particles));
    ASSERT_TRUE(writer.close());
    
    // Verify file exists and has correct number of lines
    std::ifstream file("test_output.csv");
    ASSERT_TRUE(file.is_open());
    std::string line;
    int data_lines = 0;
    while (std::getline(file, line)) {
        if (line[0] != '#' && !line.empty()) data_lines++;
    }
    EXPECT_EQ(data_lines, 101);  // header + 100 particles
}
```

**Python Verification**:
```python
import pandas as pd
df = pd.read_csv('test_output.csv', comment='#')
assert len(df) == 100
assert 'pos_x' in df.columns
assert df['id'].iloc[0] == 0
```

**Acceptance Criteria**:
- ✅ C++ test passes
- ✅ Python can read CSV with pandas
- ✅ All fields present and correct

---

## Phase 3: HDF5 Output Format (Days 3-4)

### Task 3.1: Create HDF5 Writer Class
**Estimated Time**: 4 hours

**Steps**:
1. Create `include/writers/hdf5_writer.hpp`:
   ```cpp
   #pragma once
   #include "output_metadata.hpp"
   #include "particle.hpp"
   #include <hdf5.h>
   #include <vector>
   #include <string>
   
   namespace sph {
   
   class HDF5Writer {
   public:
       struct Config {
           int compression_level = 6;  // 0-9, gzip
           bool single_file_series = false;
       };
       
       HDF5Writer(const Config& config = Config());
       ~HDF5Writer();
       
       bool open(const std::string& filepath, const OutputMetadata& metadata);
       bool write_particles(const std::vector<SPHParticle>& particles);
       bool write_metadata(const OutputMetadata& metadata);
       bool close();
       
       std::string get_extension() const { return ".h5"; }
       
       // For reading/resume
       static bool read_particles(const std::string& filepath,
                                 std::vector<SPHParticle>& particles);
       static bool read_metadata(const std::string& filepath,
                                OutputMetadata& metadata);
       
   private:
       Config m_config;
       hid_t m_file_id = -1;
       bool m_is_open = false;
       
       bool write_metadata_group(const OutputMetadata& metadata);
       bool write_particles_group(const std::vector<SPHParticle>& particles);
       bool write_energy_group(const OutputMetadata& metadata);
       
       hid_t create_compressed_dataset(hid_t group_id, const char* name,
                                       hsize_t size, hid_t datatype);
   };
   
   } // namespace sph
   ```

2. Implement in `src/writers/hdf5_writer.cpp`

**HDF5 Structure**:
```
/
├── metadata/
│   ├── simulation_name (string)
│   ├── timestamp (string)
│   ├── format_version (int)
│   ├── step (int)
│   ├── time_code (double)
│   ├── time_physical (double)
│   ├── particle_count (int)
│   ├── unit_system/  (JSON as string attribute)
│   ├── parameters/   (JSON as string attribute)
│   └── checkpoint/   (JSON as string attribute, if checkpoint)
│
├── particles/
│   ├── id [N] (int32, compressed)
│   ├── pos [N, 3] (float64, compressed)
│   ├── vel [N, 3] (float64, compressed)
│   ├── acc [N, 3] (float64, compressed)
│   ├── mass [N] (float64, compressed)
│   ├── dens [N] (float64, compressed)
│   ├── pres [N] (float64, compressed)
│   ├── ene [N] (float64, compressed)
│   ├── sml [N] (float64, compressed)
│   ├── sound [N] (float64, compressed)
│   ├── alpha [N] (float64, compressed)
│   ├── balsara [N] (float64, compressed)
│   ├── gradh [N] (float64, compressed)
│   ├── phi [N] (float64, compressed)
│   └── neighbor [N] (int32, compressed)
│
└── energy/
    ├── kinetic (double)
    ├── thermal (double)
    ├── potential (double)
    └── total (double)
```

**Implementation Example**:
```cpp
bool HDF5Writer::write_particles_group(const std::vector<SPHParticle>& particles) {
    hsize_t N = particles.size();
    hsize_t dims_scalar[1] = {N};
    hsize_t dims_vector[2] = {N, 3};
    
    // Create group
    hid_t group = H5Gcreate(m_file_id, "/particles", 
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // Write id (int32, compressed)
    std::vector<int> ids(N);
    for (size_t i = 0; i < N; ++i) ids[i] = particles[i].id;
    
    hid_t dataset = create_compressed_dataset(group, "id", N, H5T_NATIVE_INT);
    H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids.data());
    H5Dclose(dataset);
    
    // Write pos (3D vector, float64, compressed)
    std::vector<double> pos(N * 3);
    for (size_t i = 0; i < N; ++i) {
        for (int d = 0; d < 3; ++d) {
            pos[i*3 + d] = particles[i].pos[d];
        }
    }
    
    hid_t dataspace = H5Screate_simple(2, dims_vector, NULL);
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_deflate(dcpl, m_config.compression_level);
    H5Pset_chunk(dcpl, 2, dims_vector);
    
    dataset = H5Dcreate(group, "pos", H5T_NATIVE_DOUBLE, dataspace,
                       H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos.data());
    
    H5Dclose(dataset);
    H5Pclose(dcpl);
    H5Sclose(dataspace);
    
    // ... repeat for all other fields
    
    H5Gclose(group);
    return true;
}
```

**Acceptance Criteria**:
- ✅ HDF5 file created with correct structure
- ✅ All particle data written with compression
- ✅ Metadata stored as attributes
- ✅ Can read back with h5py in Python

**Files Created**:
- `include/writers/hdf5_writer.hpp`
- `src/writers/hdf5_writer.cpp`

---

### Task 3.2: Implement HDF5 Reading for Resume
**Estimated Time**: 3 hours

**Steps**:
1. Implement `HDF5Writer::read_particles()`:
   - Open file
   - Read all datasets in `/particles/`
   - Reconstruct SPHParticle vector

2. Implement `HDF5Writer::read_metadata()`:
   - Read `/metadata/` attributes
   - Deserialize JSON strings
   - Reconstruct OutputMetadata

3. Handle checkpoint data if present

**Implementation**:
```cpp
bool HDF5Writer::read_particles(const std::string& filepath,
                                std::vector<SPHParticle>& particles) {
    hid_t file = H5Fopen(filepath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return false;
    
    // Get particle count
    hid_t dataset = H5Dopen(file, "/particles/id", H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t dims[1];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    hsize_t N = dims[0];
    
    particles.resize(N);
    
    // Read id
    std::vector<int> ids(N);
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids.data());
    for (size_t i = 0; i < N; ++i) particles[i].id = ids[i];
    H5Dclose(dataset);
    H5Sclose(dataspace);
    
    // Read pos
    dataset = H5Dopen(file, "/particles/pos", H5P_DEFAULT);
    std::vector<double> pos(N * 3);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos.data());
    for (size_t i = 0; i < N; ++i) {
        for (int d = 0; d < 3; ++d) {
            particles[i].pos[d] = pos[i*3 + d];
        }
    }
    H5Dclose(dataset);
    
    // ... read all other fields similarly
    
    H5Fclose(file);
    return true;
}
```

**Acceptance Criteria**:
- ✅ Can read HDF5 files written by writer
- ✅ Particle data matches original
- ✅ Metadata correctly reconstructed
- ✅ Can resume simulation from HDF5 checkpoint

---

### Task 3.3: Test HDF5 Writer & Reader
**Estimated Time**: 2 hours

**C++ Tests**:
```cpp
TEST(HDF5Writer, WriteReadRoundTrip) {
    // Create test data
    std::vector<SPHParticle> particles_orig(1000);
    // ... initialize particles
    
    OutputMetadata meta_orig;
    // ... initialize metadata
    
    // Write
    HDF5Writer writer;
    writer.open("test.h5", meta_orig);
    writer.write_particles(particles_orig);
    writer.close();
    
    // Read back
    std::vector<SPHParticle> particles_read;
    OutputMetadata meta_read;
    HDF5Writer::read_particles("test.h5", particles_read);
    HDF5Writer::read_metadata("test.h5", meta_read);
    
    // Verify
    ASSERT_EQ(particles_read.size(), particles_orig.size());
    for (size_t i = 0; i < particles_orig.size(); ++i) {
        EXPECT_DOUBLE_EQ(particles_read[i].pos[0], particles_orig[i].pos[0]);
        EXPECT_DOUBLE_EQ(particles_read[i].mass, particles_orig[i].mass);
        // ... check all fields
    }
}
```

**Python Tests**:
```python
import h5py
import numpy as np

with h5py.File('test.h5', 'r') as f:
    # Check structure
    assert 'metadata' in f
    assert 'particles' in f
    assert 'energy' in f
    
    # Check particle data
    pos = f['particles/pos'][:]
    assert pos.shape == (1000, 3)
    
    # Check metadata
    step = f['metadata/step'][()]
    assert isinstance(step, int)
    
    # Check compression
    dataset = f['particles/mass']
    assert dataset.compression == 'gzip'
```

**Acceptance Criteria**:
- ✅ Round-trip test passes
- ✅ Python can read all fields
- ✅ Compression verified
- ✅ File size ~4x smaller than uncompressed

---

## Phase 4: Output Manager Integration (Day 5)

### Task 4.1: Create OutputManager Class
**Estimated Time**: 3 hours

**Steps**:
1. Create `include/output_manager.hpp`:
   ```cpp
   #pragma once
   #include "output_metadata.hpp"
   #include "simulation.hpp"
   #include "units.hpp"
   #include "writers/csv_writer.hpp"
   #include "writers/hdf5_writer.hpp"
   #include <memory>
   #include <fstream>
   
   namespace sph {
   
   class OutputManager {
   public:
       struct Config {
           bool enable_csv = true;
           bool enable_hdf5 = true;
           bool enable_energy_file = true;
           
           int csv_precision = 16;
           int hdf5_compression = 6;
           
           UnitSystem::Type output_units = UnitSystem::Type::CODE;
       };
       
       OutputManager(const Config& config,
                    const UnitSystem& units,
                    const std::string& output_dir);
       ~OutputManager();
       
       // Snapshot output (visualization)
       void write_snapshot(std::shared_ptr<Simulation> sim, 
                          std::shared_ptr<SPHParameters> params,
                          int count);
       
       // Checkpoint output (resume-capable)
       void write_checkpoint(std::shared_ptr<Simulation> sim,
                            std::shared_ptr<SPHParameters> params,
                            const CheckpointMetadata& meta);
       
       // Load for resume
       bool load_for_resume(const std::string& filepath,
                           std::shared_ptr<Simulation> sim,
                           OutputMetadata& meta);
       
       // Energy tracking
       void write_energy(std::shared_ptr<Simulation> sim);
       
       // Accessors
       const Config& get_config() const { return m_config; }
       const UnitSystem& get_units() const { return m_units; }
       
   private:
       Config m_config;
       UnitSystem m_units;
       std::string m_output_dir;
       
       std::unique_ptr<CSVWriter> m_csv_writer;
       std::unique_ptr<HDF5Writer> m_hdf5_writer;
       std::ofstream m_energy_file;
       
       OutputMetadata build_metadata(std::shared_ptr<Simulation> sim,
                                     std::shared_ptr<SPHParameters> params);
       
       void compute_energies(std::shared_ptr<Simulation> sim,
                           OutputMetadata& meta);
   };
   
   } // namespace sph
   ```

2. Implement in `src/output_manager.cpp`

**Key Methods**:
```cpp
void OutputManager::write_snapshot(std::shared_ptr<Simulation> sim,
                                   std::shared_ptr<SPHParameters> params,
                                   int count) {
    // Build metadata
    OutputMetadata meta = build_metadata(sim, params);
    meta.step = count;
    meta.is_checkpoint = false;
    
    compute_energies(sim, meta);
    
    const auto& particles = sim->get_particles();
    
    // Write CSV if enabled
    if (m_config.enable_csv) {
        std::string csv_path = m_output_dir + "/snapshot_" 
                             + std::to_string(count) + ".csv";
        m_csv_writer->open(csv_path, meta);
        m_csv_writer->write_particles(particles);
        m_csv_writer->close();
    }
    
    // Write HDF5 if enabled
    if (m_config.enable_hdf5) {
        std::string h5_path = m_output_dir + "/snapshot_" 
                            + std::to_string(count) + ".h5";
        m_hdf5_writer->open(h5_path, meta);
        m_hdf5_writer->write_particles(particles);
        m_hdf5_writer->close();
    }
}

void OutputManager::write_checkpoint(std::shared_ptr<Simulation> sim,
                                     std::shared_ptr<SPHParameters> params,
                                     const CheckpointMetadata& chk_meta) {
    OutputMetadata meta = build_metadata(sim, params);
    meta.is_checkpoint = true;
    meta.checkpoint_data = chk_meta;
    
    compute_energies(sim, meta);
    
    // Always use HDF5 for checkpoints (complete + compressed)
    std::string h5_path = m_output_dir + "/checkpoint_" 
                        + std::to_string(chk_meta.relaxation_step) + ".h5";
    
    m_hdf5_writer->open(h5_path, meta);
    m_hdf5_writer->write_particles(sim->get_particles());
    m_hdf5_writer->close();
}
```

**Acceptance Criteria**:
- ✅ OutputManager integrates all writers
- ✅ Can write snapshots in multiple formats
- ✅ Checkpoints use HDF5 format
- ✅ Metadata correctly populated

**Files Created**:
- `include/output_manager.hpp`
- `src/output_manager.cpp`

---

### Task 4.2: Update Solver to Use OutputManager
**Estimated Time**: 2 hours

**Steps**:
1. Modify `include/solver.hpp`:
   ```cpp
   class Solver {
       // ... existing members ...
       
       std::shared_ptr<OutputManager> m_output_manager;
       OutputManager::Config m_output_config;
       UnitSystem m_unit_system;
       
       // Remove old m_output (Output class)
   };
   ```

2. Modify `src/solver.cpp`:
   - Initialize OutputManager in constructor
   - Replace `m_output->output_particle()` with `m_output_manager->write_snapshot()`
   - Replace checkpoint saves with `m_output_manager->write_checkpoint()`
   - Update resume logic to use `m_output_manager->load_for_resume()`

**Configuration Reading**:
```cpp
void Solver::read_parameterfile(const char* filename) {
    // ... existing parameter reading ...
    
    // Read unit system configuration
    std::string unit_system_str = input.get<std::string>("units.system", "code");
    if (unit_system_str == "galactic") {
        real length_kpc = input.get<real>("units.code_length", 1.0);
        real mass_msun = input.get<real>("units.code_mass", 1.0);
        real vel_kms = input.get<real>("units.code_velocity", 1.0);
        m_unit_system = UnitSystem::create_galactic(length_kpc, mass_msun, vel_kms);
    } else if (unit_system_str == "si") {
        m_unit_system = UnitSystem::create_si();
    } else if (unit_system_str == "cgs") {
        m_unit_system = UnitSystem::create_cgs();
    } else {
        m_unit_system = UnitSystem::create_code();
    }
    
    // Read output configuration
    m_output_config.enable_csv = input.get<bool>("output.csv.enabled", true);
    m_output_config.enable_hdf5 = input.get<bool>("output.hdf5.enabled", true);
    m_output_config.csv_precision = input.get<int>("output.csv.precision", 16);
    m_output_config.hdf5_compression = input.get<int>("output.hdf5.compression", 6);
}

void Solver::initialize() {
    // ... existing initialization ...
    
    // Create OutputManager instead of Output
    m_output_manager = std::make_shared<OutputManager>(
        m_output_config, m_unit_system, m_output_dir
    );
}
```

**Update Run Loop**:
```cpp
void Solver::run() {
    // ...
    
    if (t > t_out) {
        m_output_manager->write_snapshot(m_sim, m_param, loop);
        t_out += m_param->time.output;
    }
    
    if (t > t_ene) {
        m_output_manager->write_energy(m_sim);
        t_ene += m_param->time.energy;
    }
    
    // Checkpoint
    if (m_checkpoint_config.enabled && should_save(loop)) {
        CheckpointMetadata chk_meta;
        // ... fill checkpoint metadata ...
        m_output_manager->write_checkpoint(m_sim, m_param, chk_meta);
    }
}
```

**Acceptance Criteria**:
- ✅ Solver compiles with OutputManager
- ✅ Old Output class removed
- ✅ Configuration read from JSON
- ✅ Snapshots written in new format
- ✅ Checkpoints use HDF5

**Files Modified**:
- `include/solver.hpp`
- `src/solver.cpp`

---

### Task 4.3: Remove Old Checkpoint System
**Estimated Time**: 30 minutes

**Steps**:
1. Delete `include/checkpoint.hpp` completely
2. Delete `src/checkpoint.cpp` completely
3. Remove references from CMakeLists.txt
4. No backward compatibility code needed

**Files Deleted**:
- `include/checkpoint.hpp`
- `src/checkpoint.cpp`

**Files Modified**:
- `CMakeLists.txt` (remove checkpoint.cpp)
- `src/CMakeLists.txt` (if needed)

**Acceptance Criteria**:
- ✅ No checkpoint.hpp/cpp files exist
- ✅ No `.chk` files generated
- ✅ Clean compilation
- ✅ No backward compatibility code

---

## Phase 5: Documentation & Cleanup (Day 6)

### Task 5.1: Create User Documentation
**Estimated Time**: 2 hours

**Steps**:
1. Create `docs/USER_GUIDE_OUTPUT.md`:
   - How to configure output formats
   - How to set unit systems  
   - Example configurations
   - Python reading examples

2. Update `README.md` with new output section

**Acceptance Criteria**:
- ✅ User guide complete with examples
- ✅ README updated

**Files Created**:
- `docs/USER_GUIDE_OUTPUT.md`

**Files Modified**:
- `README.md`

---

### Task 5.2: Update Example Configurations
**Estimated Time**: 1 hour

**Steps**:
1. Update all JSON configs to use new output format:
   - `sample/lane_emden/lane_emden.json`
   - `sample/evrard/evrard.json`
   - Any other example configs

2. Add output and units sections:
   ```json
   {
     "units": {
       "system": "code"
     },
     "output": {
       "formats": ["csv", "hdf5"],
       "csv": {
         "enabled": true,
         "precision": 16
       },
       "hdf5": {
         "enabled": true,
         "compression": 6
       }
     }
   }
   ```

3. Remove old `checkpoint` sections from configs

**Acceptance Criteria**:
- ✅ All example configs updated
- ✅ No references to old checkpoint format
- ✅ Configs use new output format

**Files Modified**:
- All `sample/*/config.json` files
- All `sample/*/*.json` files


---

### Task 5.3: Remove Legacy Files and Folders
**Estimated Time**: 2 hours

**Steps**:

1. **Delete old checkpoint system completely**:
   ```bash
   # Remove checkpoint header and source
   git rm include/checkpoint.hpp
   git rm src/checkpoint.cpp
   ```

2. **Delete old output system**:
   ```bash
   # Remove old output implementation
   git rm include/output.hpp
   git rm src/output.cpp
   ```

3. **Clean up legacy documentation**:
   ```bash
   # Remove old checkpoint docs
   git rm CHECKPOINT_IMPLEMENTATION.md
   git rm docs/old_relaxation_config.json.bak
   ```

4. **Remove legacy data files** (if any exist):
   ```bash
   # Clean up any old .dat and .chk files
   find . -name "*.chk" -type f -delete
   find . -name "*.dat" -type f -delete
   find . -name "*checkpoint*.json" -type f -delete
   ```

5. **Update CMakeLists.txt** to remove references:
   ```cmake
   # Remove these lines if they exist:
   # ${CMAKE_SOURCE_DIR}/src/checkpoint.cpp
   # ${CMAKE_SOURCE_DIR}/src/output.cpp
   ```

6. **Remove old conversion/migration tools** (since no backward compatibility needed):
   ```bash
   rm -rf tools/convert_old_outputs.py  # Created in Task 5.1 but not needed
   rm -rf docs/MIGRATION_GUIDE.md       # Not needed without migration
   ```

7. **Clean sample/lane_emden/** of old outputs:
   ```bash
   cd sample/lane_emden
   find . -name "*.dat" -delete
   find . -name "*.chk" -delete
   find . -name "checkpoint_*.json" -delete
   ```

8. **Update all JSON configs** to remove checkpoint section (use output instead):
   - Edit `sample/lane_emden/lane_emden.json`
   - Remove `checkpoint` section completely
   - Ensure only `output` section exists

9. **Remove backward compatibility code** from new implementations:
   - In `HDF5Writer::read()`, remove any code for reading old .chk format
   - In `OutputManager`, remove checkpoint migration logic

10. **Update .gitignore**:
    ```bash
    # Add to .gitignore to prevent accidentally committing old formats
    echo "*.chk" >> .gitignore
    echo "*.dat" >> .gitignore
    echo "checkpoint_*.json" >> .gitignore
    ```

**Files to Delete**:
- `include/checkpoint.hpp` - Old binary checkpoint system
- `src/checkpoint.cpp` - Old checkpoint implementation
- `include/output.hpp` - Old ASCII output system
- `src/output.cpp` - Old output implementation
- `CHECKPOINT_IMPLEMENTATION.md` - Old documentation
- `docs/old_relaxation_config.json.bak` - Legacy backup
- `tools/convert_old_outputs.py` - Not needed (no migration)
- `docs/MIGRATION_GUIDE.md` - Not needed (no backward compatibility)
- All `*.chk`, `*.dat`, `checkpoint_*.json` files in data directories

**Directories to Clean**:
- `lane_emden/results/polytrope_n1.5_3d/*.dat` - Old outputs
- `lane_emden/results/polytrope_n1.5_3d/*.chk` - Old checkpoints
- Any other result directories with legacy formats

**Acceptance Criteria**:
- ✅ No checkpoint.hpp/cpp files exist
- ✅ No output.hpp/cpp files exist
- ✅ No .chk or .dat files in repository
- ✅ CMakeLists.txt compiles without legacy files
- ✅ All sample configs updated to new format
- ✅ git status shows clean removal
- ✅ No backward compatibility code remains
- ✅ .gitignore prevents legacy formats

**Files Modified**:
- `CMakeLists.txt` (remove old source references)
- `src/CMakeLists.txt` (remove old source references)
- `sample/lane_emden/lane_emden.json` (remove checkpoint section)
- `.gitignore` (add legacy format patterns)

**Files Deleted**:
- `include/checkpoint.hpp`
- `src/checkpoint.cpp`
- `include/output.hpp`
- `src/output.cpp`
- `CHECKPOINT_IMPLEMENTATION.md`
- `docs/old_relaxation_config.json.bak`
- `tools/convert_old_outputs.py`
- `docs/MIGRATION_GUIDE.md`
- All `*.chk`, `*.dat` files in results

**Verification Commands**:
```bash
# Verify no legacy files remain
find . -name "*.chk" -o -name "checkpoint.hpp" -o -name "output.hpp"
# Should return nothing

# Verify clean compile
cd build && cmake .. && make clean && make -j8
# Should compile successfully without errors

# Verify git status
git status
# Should show cleanly deleted files
```

---

## Phase 6: Testing & Validation (Day 7)

### Task 6.1: Integration Tests
**Estimated Time**: 3 hours

**Steps**:
1. Run Lane-Emden simulation with new output
2. Verify outputs match expected format
3. Test resume from checkpoint
4. Compare numerical results with old format

**Test Script**:
```bash
#!/bin/bash
# test/integration_test_output.sh

echo "=== Integration Test: New Output System ==="

# Test 1: Run with new output
echo "Test 1: Running simulation with CSV+HDF5 output..."
./build/sph lane_emden
if [ $? -ne 0 ]; then
    echo "FAILED: Simulation crashed"
    exit 1
fi

# Test 2: Verify CSV format
echo "Test 2: Verifying CSV output..."
python3 -c "
import pandas as pd
df = pd.read_csv('lane_emden/results/polytrope_n1.5_3d/snapshot_00000.csv', comment='#')
assert len(df) == 5400, 'Wrong particle count'
assert 'pos_x' in df.columns, 'Missing pos_x column'
print('CSV format OK')
"

# Test 3: Verify HDF5 format
echo "Test 3: Verifying HDF5 output..."
python3 -c "
import h5py
with h5py.File('lane_emden/results/polytrope_n1.5_3d/snapshot_00000.h5', 'r') as f:
    assert 'metadata' in f, 'Missing metadata group'
    assert 'particles' in f, 'Missing particles group'
    assert f['particles/pos'].shape == (5400, 3), 'Wrong pos shape'
    print('HDF5 format OK')
"

# Test 4: Resume from checkpoint
echo "Test 4: Testing resume from checkpoint..."
# Modify config to resume
python3 -c "
import json
with open('sample/lane_emden/lane_emden.json', 'r') as f:
    config = json.load(f)
config['checkpoint']['resumeFile'] = 'lane_emden/results/polytrope_n1.5_3d/checkpoint_005000.h5'
with open('sample/lane_emden/lane_emden_resume_test.json', 'w') as f:
    json.dump(config, f, indent=2)
"

cp sample/lane_emden/lane_emden_resume_test.json sample/lane_emden/lane_emden.json
./build/sph lane_emden

if [ $? -ne 0 ]; then
    echo "FAILED: Resume failed"
    exit 1
fi

echo "=== All tests passed ==="
```

**Acceptance Criteria**:
- ✅ All integration tests pass
- ✅ Resume works correctly
- ✅ Numerical results consistent

---

### Task 6.2: Performance Benchmarks
**Estimated Time**: 2 hours

**Steps**:
1. Benchmark write times for different formats
2. Benchmark read times
3. Compare file sizes
4. Document results

**Benchmark Script**:
```python
# test/benchmark_output.py
import time
import h5py
import pandas as pd
import subprocess
import json

def benchmark_write():
    # Run simulation with timing
    configs = [
        ('csv_only', {'csv': True, 'hdf5': False}),
        ('hdf5_only', {'csv': False, 'hdf5': True}),
        ('both', {'csv': True, 'hdf5': True}),
    ]
    
    results = {}
    for name, formats in configs:
        # Update config
        # Run simulation
        # Measure time
        results[name] = time_taken
    
    return results

def benchmark_read():
    # Time reading CSV
    # Time reading HDF5
    # Compare
    pass

if __name__ == '__main__':
    write_results = benchmark_write()
    read_results = benchmark_read()
    
    print("Write Performance:")
    print(json.dumps(write_results, indent=2))
    print("\nRead Performance:")
    print(json.dumps(read_results, indent=2))
```

**Acceptance Criteria**:
- ✅ Performance metrics documented
- ✅ HDF5 faster than CSV for read
- ✅ HDF5 smaller than CSV (with compression)

---

### Task 6.3: Update Example Configurations
**Estimated Time**: 1 hour

**Steps**:
1. Update all sample JSON files with new output/unit sections
2. Create example configurations for different use cases
3. Test each example

**Example Configs**:
```json
// sample/lane_emden/lane_emden_galactic.json
{
  "units": {
    "system": "galactic",
    "code_length": 10.0,
    "code_mass": 1e6,
    "code_velocity": 1.0
  },
  "output": {
    "formats": ["hdf5"],
    "hdf5": {
      "enabled": true,
      "compression": 9
    }
  }
}
```

**Acceptance Criteria**:
- ✅ All samples updated
- ✅ Each sample tested
- ✅ Documentation matches examples

**Files Modified**:
- All `sample/**/*.json` files

---

## Summary Checklist

### Dependencies
- [ ] HDF5 library installed
- [ ] nlohmann/json integrated
- [ ] CMake configuration updated
- [ ] Project compiles

### Core Implementation
- [ ] UnitSystem class complete
- [ ] CSVWriter complete
- [ ] HDF5Writer complete
- [ ] OutputManager complete
- [ ] Solver integration complete
- [ ] Old checkpoint system removed

### Testing
- [ ] Unit tests pass (units, CSV, HDF5)
- [ ] Integration tests pass
- [ ] Resume from checkpoint works
- [ ] Performance benchmarks done

### Documentation
- [ ] Design document complete ✅
- [ ] User guide written
- [ ] README updated
- [ ] Code comments added

### Cleanup
- [ ] Legacy files removed (checkpoint.hpp/cpp, output.hpp/cpp)
- [ ] Old .dat and .chk files deleted
- [ ] Example configs updated to new format
- [ ] No backward compatibility code remains

---

## Post-Implementation Tasks

1. **Code Review**:
   - Review all new code for best practices
   - Check memory leaks (valgrind)
   - Verify thread safety (OpenMP)

2. **User Testing**:
   - Test with Lane-Emden relaxation
   - Verify Python analysis tools work
   - Test resume capability

3. **Performance Optimization**:
   - Profile I/O operations
   - Optimize hot paths
   - Consider parallel HDF5 for future MPI support

4. **Future Enhancements**:
   - VTK writer for ParaView
   - Parallel I/O support
   - In-situ visualization hooks

---

## Estimated Timeline

| Phase | Days | Tasks |
|-------|------|-------|
| Phase 1: Foundation | 1 | Dependencies, UnitSystem |
| Phase 2: CSV Output | 1 | CSVWriter, tests |
| Phase 3: HDF5 Output | 2 | HDF5Writer, reader, tests |
| Phase 4: Integration | 1 | OutputManager, Solver update |
| Phase 5: Documentation & Cleanup | 1 | Docs, config updates, legacy removal |
| Phase 6: Testing | 1 | Integration tests, benchmarks |
| **Total** | **7 days** | |

---

## Risk Mitigation

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| HDF5 build issues | Medium | High | Provide Docker container with dependencies |
| Performance regression | Low | Medium | Benchmark early, optimize as needed |
| User adoption | Medium | Low | Excellent documentation, examples |
| Breaking changes | High | Low | **No backward compatibility needed**, clean break |
| Bugs in new code | Medium | Medium | Comprehensive testing, code review |

---

## Success Criteria

✅ **Functionality**:
- All output formats working
- Resume capability maintained
- Unit conversions accurate

✅ **Performance**:
- Write time < 2x old system
- HDF5 file size < 50% of old .dat
- Read time faster than old system

✅ **Usability**:
- Easy configuration
- Clear documentation
- Python examples work

✅ **Quality**:
- All tests passing
- No memory leaks
- Thread-safe

---

*End of Implementation Plan*
