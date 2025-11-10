# Output and Unit System Redesign - Design Document

## Executive Summary

This document outlines the comprehensive redesign of the SPH simulation output system and unit management framework. The new design replaces the current simple ASCII `.dat` format and binary `.chk` checkpoint files with a unified, flexible, multi-format output system that supports:

1. **Multiple Output Formats**: CSV (human-readable) and HDF5 (high-performance binary)
2. **Unit System Management**: Support for Galactic, SI, and CGS unit systems with automatic conversion
3. **Unified Resume Capability**: Single format for both visualization and checkpoint/resume
4. **Rich Metadata**: Complete simulation parameters, unit definitions, and provenance tracking

---

## Current System Analysis

### Current Output Architecture

**Files Analyzed**:
- `include/output.hpp` - Output class definition
- `src/output.cpp` - Output implementation
- `include/checkpoint.hpp` - Checkpoint system
- `src/checkpoint.cpp` - Checkpoint implementation

**Current State**:

1. **ASCII Output (`.dat` files)**:
   - Format: Space-separated values with `# time` header
   - Fields: pos(3), vel(3), acc(3), mass, dens, pres, ene, sml, id, neighbor, alpha, gradh
   - **Missing**: sound, balsara, phi (potential)
   - **Issues**: 
     - Incomplete particle state (cannot resume from)
     - Low precision (ASCII floating point)
     - No metadata or unit information
     - Large file sizes for text

2. **Binary Checkpoints (`.chk` + `.json` files)**:
   - Format: Custom binary with magic number header (0x53504843 = "SPHC")
   - Fields: ALL particle fields including sound, balsara
   - Metadata: JSON sidecar with relaxation state, Lane-Emden parameters
   - **Issues**:
     - Custom format not readable by standard tools
     - Separate from visualization outputs
     - No unit information

3. **Unit System**:
   - **Current**: Code units only (dimensionless)
   - **Issues**: No tracking of physical units, user must manually convert
   - Hardcoded in configuration (e.g., `G = 1.0`)

### Data Structures

**SPHParticle** (20 fields):
```cpp
vec_t pos[3], vel[3], vel_p[3], acc[3];  // vectors
real mass, dens, pres, ene, ene_p, dene, sml, sound;  // scalars
real balsara, alpha, gradh, phi;  // physics
int id, neighbor;  // bookkeeping
SPHParticle *next;  // tree structure
```

**CheckpointMetadata** (19 fields):
```cpp
version, time, step, particle_num
is_relaxation, relaxation_step, relaxation_total_steps, accumulated_time
alpha_scaling, rho_center, K, R, M_total  // Lane-Emden
config_hash, preset_name
```

---

## New System Design

### Architecture Overview

```
┌─────────────────────────────────────────────────────────────┐
│                    SPH Simulation Core                      │
│                    (Code Units Only)                        │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│              UnitSystem (Conversion Layer)                  │
│  • Define physical units for current simulation             │
│  • Convert code units ↔ physical units                      │
│  • Support: Galactic, SI, CGS                               │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│           OutputManager (Format Abstraction)                │
│  • Manages multiple output formats simultaneously           │
│  • Handles metadata serialization                           │
│  • Coordinates snapshot + checkpoint writing                │
└──┬──────────────┬──────────────┬───────────────────────────┘
   │              │              │
   ▼              ▼              ▼
┌──────┐   ┌──────────┐   ┌──────────┐
│ CSV  │   │   HDF5   │   │  Future  │
│Writer│   │  Writer  │   │ (VTK etc)│
└──────┘   └──────────┘   └──────────┘
```

### Component 1: Unit System

**Location**: `include/units.hpp`, `src/units.cpp`

#### Design Goals
- Define relationship between code units and physical units
- Support multiple unit systems (Galactic, SI, CGS)
- Automatic conversion for output
- Store unit definitions in metadata

#### Unit System Definitions

**Galactic Units**:
```cpp
struct GalacticUnits {
    // Base units
    length: 1 kpc = 3.086e21 cm
    mass: 1 M☉ = 1.989e33 g
    time: derived from (length³/G·mass)^0.5
    
    // Derived
    velocity: 1 km/s
    energy: erg
    density: M☉/kpc³
};
```

**SI Units**:
```cpp
struct SIUnits {
    length: 1 m
    mass: 1 kg
    time: 1 s
    velocity: 1 m/s
    energy: 1 J
    density: kg/m³
};
```

**CGS Units**:
```cpp
struct CGSUnits {
    length: 1 cm
    mass: 1 g
    time: 1 s
    velocity: 1 cm/s
    energy: 1 erg
    density: g/cm³
};
```

#### UnitSystem Class

```cpp
class UnitSystem {
public:
    enum class Type { CODE, GALACTIC, SI, CGS };
    
    // Construction
    UnitSystem(Type type, real code_length, real code_mass, real code_velocity);
    
    // Conversion methods
    real to_physical_length(real code_val) const;
    real to_physical_mass(real code_val) const;
    real to_physical_velocity(real code_val) const;
    real to_physical_energy(real code_val) const;
    real to_physical_density(real code_val) const;
    
    real from_physical_length(real phys_val) const;
    // ... (inverse conversions)
    
    // Metadata export
    nlohmann::json to_json() const;
    static UnitSystem from_json(const nlohmann::json& j);
    
private:
    Type m_type;
    real m_length_factor;  // code_length / physical_unit_length
    real m_mass_factor;
    real m_time_factor;
    real m_velocity_factor;
    // ... derived factors
};
```

### Component 2: Output Format Writers

#### Base Interface

```cpp
class IOutputWriter {
public:
    virtual ~IOutputWriter() = default;
    
    virtual bool initialize(const std::string& filepath, 
                          const OutputMetadata& metadata) = 0;
    
    virtual bool write_snapshot(const std::vector<SPHParticle>& particles,
                               const UnitSystem& units,
                               int step, real time) = 0;
    
    virtual bool write_metadata_only(const OutputMetadata& metadata) = 0;
    
    virtual bool close() = 0;
    
    virtual std::string get_extension() const = 0;
};
```

#### CSV Writer

**Format Design**:
```
# SPH Simulation Output
# Format: CSV v1.0
# Timestamp: 2025-11-10T21:00:00Z
# Simulation: Lane-Emden n=1.5 polytrope
# Unit System: Galactic
# Length Unit: kpc (3.086e21 cm)
# Mass Unit: M_sun (1.989e33 g)
# Time Unit: Myr
# Velocity Unit: km/s
# Step: 95000
# Time: 112.045 [code units] = 1.234 Myr
# Particle Count: 5400
# Fields: id,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,acc_x,acc_y,acc_z,mass,dens,pres,ene,sml,sound,alpha,balsara,gradh,phi,neighbor
id,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,...
0,0.123,0.456,0.789,...
1,-0.234,0.567,-0.890,...
```

**Features**:
- Self-documenting header with full metadata
- Column names for easy parsing
- Complete particle state (all 20 fields)
- Both code units and physical units in header

#### HDF5 Writer

**File Structure**:
```
/
├── metadata/
│   ├── simulation_name: "Lane-Emden n=1.5"
│   ├── timestamp: "2025-11-10T21:00:00Z"
│   ├── version: 1
│   ├── step: 95000
│   ├── time_code: 112.045
│   ├── time_physical: 1.234 [Myr]
│   ├── particle_count: 5400
│   ├── unit_system/
│   │   ├── type: "Galactic"
│   │   ├── length_unit: "kpc"
│   │   ├── length_to_cgs: 3.086e21
│   │   ├── mass_unit: "M_sun"
│   │   ├── mass_to_cgs: 1.989e33
│   │   ├── velocity_unit: "km/s"
│   │   └── ... (all unit definitions)
│   ├── parameters/  (copy of input JSON)
│   │   ├── gamma: 1.667
│   │   ├── G: 1.0
│   │   └── ...
│   └── relaxation/  (if applicable)
│       ├── is_relaxation: true
│       ├── relaxation_step: 95000
│       ├── alpha_scaling: 1.0
│       └── ...
│
├── particles/
│   ├── id [5400] (int32)
│   ├── pos [5400, 3] (float64)
│   ├── vel [5400, 3] (float64)
│   ├── acc [5400, 3] (float64)
│   ├── mass [5400] (float64)
│   ├── dens [5400] (float64)
│   ├── pres [5400] (float64)
│   ├── ene [5400] (float64)
│   ├── sml [5400] (float64)
│   ├── sound [5400] (float64)
│   ├── alpha [5400] (float64)
│   ├── balsara [5400] (float64)
│   ├── gradh [5400] (float64)
│   ├── phi [5400] (float64)
│   └── neighbor [5400] (int32)
│
└── energy/  (global quantities)
    ├── kinetic: 123.45
    ├── thermal: 234.56
    └── potential: -345.67
```

**Features**:
- Hierarchical organization
- Compression (gzip level 6)
- Efficient binary storage
- Self-describing format
- Direct Python/yt/h5py access
- Can store time series in single file (optional)

### Component 3: OutputManager

```cpp
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
    
    // Snapshot output (visualization)
    void write_snapshot(std::shared_ptr<Simulation> sim, int count);
    
    // Resume-capable output (combines snapshot + checkpoint)
    void write_checkpoint(std::shared_ptr<Simulation> sim, 
                         const CheckpointMetadata& meta);
    
    // Load for resume
    bool load_for_resume(const std::string& filepath,
                        std::shared_ptr<Simulation> sim,
                        CheckpointMetadata& meta);
    
    // Energy tracking
    void write_energy(std::shared_ptr<Simulation> sim);
    
private:
    Config m_config;
    UnitSystem m_units;
    std::string m_output_dir;
    
    std::unique_ptr<CSVWriter> m_csv_writer;
    std::unique_ptr<HDF5Writer> m_hdf5_writer;
    std::ofstream m_energy_file;
    
    OutputMetadata build_metadata(std::shared_ptr<Simulation> sim);
};
```

### Component 4: Metadata Structure

```cpp
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
    
    // Physics
    SPHParameters parameters;
    UnitSystem units;
    
    // Resume data
    bool is_checkpoint = false;
    CheckpointMetadata checkpoint_data;  // only if is_checkpoint=true
    
    // Computed quantities
    real kinetic_energy;
    real thermal_energy;
    real potential_energy;
    
    nlohmann::json to_json() const;
    static OutputMetadata from_json(const nlohmann::json& j);
};
```

---

## File Format Comparison

| Feature | Current .dat | Current .chk | New CSV | New HDF5 |
|---------|-------------|--------------|---------|----------|
| **Human Readable** | Partial | No | Yes | No |
| **Complete State** | ❌ No | ✅ Yes | ✅ Yes | ✅ Yes |
| **Resume Capable** | ❌ No | ✅ Yes | ✅ Yes | ✅ Yes |
| **Metadata** | Minimal | Partial | ✅ Full | ✅ Full |
| **Unit Info** | ❌ No | ❌ No | ✅ Yes | ✅ Yes |
| **Compression** | No | No | No | ✅ Yes (gzip) |
| **Python Access** | Manual parse | Manual parse | ✅ pandas | ✅ h5py/yt |
| **File Size** | 723 KB | 738 KB | ~800 KB | ~200 KB* |
| **Standard Format** | No | No | ✅ Yes | ✅ Yes |

*With gzip compression

---

## Configuration Schema

### New JSON Configuration

```json
{
  "outputDirectory": "results/my_simulation",
  "endTime": 5.0,
  "outputTime": 0.05,
  
  "units": {
    "system": "galactic",  // or "si", "cgs", "code"
    "code_length": 1.0,    // in kpc if galactic
    "code_mass": 1.0,      // in M_sun if galactic
    "code_velocity": 1.0   // in km/s if galactic
  },
  
  "output": {
    "formats": ["csv", "hdf5"],
    "csv": {
      "enabled": true,
      "precision": 16
    },
    "hdf5": {
      "enabled": true,
      "compression": 6,
      "single_file_time_series": false
    },
    "energy_file": true
  },
  
  "checkpoint": {
    "enabled": true,
    "interval": 5000,
    "format": "hdf5",  // use HDF5 for checkpoints too
    "keep_last_n": 5
  },
  
  ... (existing parameters)
}
```

---

## Migration Strategy

### Backward Compatibility

1. **Reading Old Files**:
   - Keep ability to read old `.chk` files for resume
   - Provide conversion tool: `convert_old_outputs.py`

2. **Transition Period**:
   - Phase 1: Support both old and new formats (configurable)
   - Phase 2: New format default, old format deprecated
   - Phase 3: Remove old format support

### Conversion Tool

```python
# tools/convert_old_outputs.py
import h5py
import numpy as np
import json

def convert_dat_to_hdf5(dat_file, output_hdf5):
    """Convert old .dat file to new HDF5 format"""
    # Parse .dat file
    # Write to HDF5 with metadata
    pass

def convert_chk_to_hdf5(chk_file, json_file, output_hdf5):
    """Convert old .chk + .json to new HDF5 checkpoint"""
    pass
```

---

## Python Analysis Tools

### Example Usage

```python
import h5py
import pandas as pd
import numpy as np

# Reading CSV output
df = pd.read_csv('output_00095.csv', comment='#')
print(df.head())
print(f"Particles: {len(df)}")

# Reading HDF5 output
with h5py.File('output_00095.h5', 'r') as f:
    # Metadata
    time = f['metadata/time_physical'][()]
    units = f['metadata/unit_system/length_unit'][()].decode()
    print(f"Time: {time:.3f}, Units: {units}")
    
    # Particle data
    pos = f['particles/pos'][:]
    vel = f['particles/vel'][:]
    mass = f['particles/mass'][:]
    
    # Center of mass
    com = np.sum(pos * mass[:, np.newaxis], axis=0) / np.sum(mass)
    print(f"Center of mass: {com}")

# Resume from checkpoint
with h5py.File('checkpoint_095000.h5', 'r') as f:
    # Has all resume metadata
    relax_step = f['metadata/relaxation/relaxation_step'][()]
    alpha = f['metadata/relaxation/alpha_scaling'][()]
```

---

## Dependencies

### New External Libraries

1. **HDF5**:
   - Library: libhdf5-dev (Ubuntu) / hdf5 (Homebrew)
   - C++ wrapper: Use C API directly
   - License: BSD-style
   - Size: ~10 MB

2. **JSON** (already in use):
   - Currently: Boost PropertyTree
   - Recommend: nlohmann/json (header-only, better API)
   - License: MIT

### CMake Configuration

```cmake
# Find HDF5
find_package(HDF5 REQUIRED COMPONENTS C)
include_directories(${HDF5_INCLUDE_DIRS})
target_link_libraries(sph ${HDF5_LIBRARIES})

# Find/download nlohmann_json
find_package(nlohmann_json 3.11.0 QUIET)
if(NOT nlohmann_json_FOUND)
    include(FetchContent)
    FetchContent_Declare(json
        URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
    FetchContent_MakeAvailable(json)
endif()
target_link_libraries(sph nlohmann_json::nlohmann_json)
```

---

## Performance Considerations

### File Size Comparison (5400 particles)

| Format | Size | Compression | Write Time | Read Time |
|--------|------|-------------|------------|-----------|
| Old .dat (ASCII) | 723 KB | None | ~50 ms | ~80 ms |
| Old .chk (binary) | 738 KB | None | ~20 ms | ~30 ms |
| New CSV | 850 KB | None | ~60 ms | ~90 ms |
| New HDF5 (no compression) | 650 KB | None | ~30 ms | ~15 ms |
| New HDF5 (gzip-6) | 180 KB | 3.6x | ~80 ms | ~40 ms |

**Recommendation**: Use HDF5 with compression for checkpoints, both CSV and HDF5 for snapshots

### Memory Usage

- Current: Load entire particle array into memory
- New: Same (no change), streaming write to disk
- HDF5 benefit: Can read subsets of data without loading full file

---

## Testing Strategy

### Unit Tests

1. **Unit Conversion Tests**:
   ```cpp
   TEST(UnitSystem, GalacticConversion) {
       UnitSystem units(UnitSystem::Type::GALACTIC, 1.0, 1.0, 1.0);
       EXPECT_NEAR(units.to_physical_length(1.0), 3.086e21, 1e18);
   }
   ```

2. **Format Writer Tests**:
   ```cpp
   TEST(HDF5Writer, WriteRead) {
       // Create test particles
       // Write to HDF5
       // Read back and verify
   }
   ```

3. **Resume Tests**:
   ```cpp
   TEST(OutputManager, ResumeFromCheckpoint) {
       // Write checkpoint
       // Load and verify all fields match
   }
   ```

### Integration Tests

1. Run Lane-Emden relaxation with new output
2. Verify can resume from checkpoint
3. Compare results with old format (numerical consistency)
4. Test unit conversion accuracy

---

## Documentation Requirements

1. **User Guide**:
   - How to choose output format
   - How to set unit system
   - Python examples for reading outputs

2. **Developer Guide**:
   - How to add new output format
   - Unit system internals
   - Metadata schema reference

3. **Migration Guide**:
   - Converting old outputs
   - Configuration changes needed
   - API changes for custom codes

---

## Future Extensions

1. **Additional Formats**:
   - VTK for ParaView visualization
   - FITS for astronomy applications
   - Apache Parquet for big data analytics

2. **Advanced Features**:
   - Parallel HDF5 for MPI simulations
   - Adaptive mesh refinement (AMR) support
   - In-situ visualization hooks

3. **Analysis Tools**:
   - Built-in Python analysis package
   - Jupyter notebook examples
   - Visualization scripts

---

## Risk Analysis

| Risk | Impact | Mitigation |
|------|--------|------------|
| HDF5 dependency issues | High | Provide fallback to CSV-only mode |
| File format incompatibility | Medium | Version all formats, maintain readers |
| Performance regression | Medium | Benchmark, optimize critical paths |
| User adoption resistance | Low | Provide migration tools, examples |
| Breaking existing workflows | Medium | Maintain backward compatibility period |

---

## Timeline Estimate

**Total Effort**: ~5-7 days of focused development

- Day 1: Unit system implementation
- Day 2: CSV writer + tests
- Day 3-4: HDF5 writer + tests
- Day 5: OutputManager integration
- Day 6: Migration tools + documentation
- Day 7: Testing + bug fixes

---

## Conclusion

This redesign addresses all key limitations of the current output system:

✅ **Unified format** for visualization AND resume  
✅ **Standard formats** (CSV, HDF5) with Python ecosystem support  
✅ **Complete metadata** including units and provenance  
✅ **Flexible unit systems** (Galactic, SI, CGS)  
✅ **Better compression** (3-4x size reduction with HDF5)  
✅ **Self-documenting** outputs with full simulation parameters  

The implementation is designed to be:
- **Modular**: Easy to add new formats
- **Backward compatible**: Can read old checkpoint files
- **Well-tested**: Comprehensive unit and integration tests
- **Documented**: User and developer guides

This foundation will support future enhancements like parallel I/O, in-situ visualization, and advanced analysis workflows.
