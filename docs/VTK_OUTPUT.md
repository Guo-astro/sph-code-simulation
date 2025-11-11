# VTK Output Format

## Overview

The SPH code now supports VTK (Visualization Toolkit) Legacy format output for visualization in ParaView, VisIt, and other scientific visualization tools.

## Configuration

VTK output is configured through the `output.formats` array in your JSON configuration:

```json
{
  "output": {
    "formats": [
      {"type": "csv", "precision": 16},
      {"type": "hdf5", "compression": 6},
      {"type": "vtk", "binary": true}
    ],
    "enableEnergyFile": true
  }
}
```

### VTK Format Options

- `type`: Must be `"vtk"`
- `binary`: Boolean flag
  - `true`: Binary VTK format (smaller files, faster I/O)
  - `false`: ASCII VTK format (human-readable, easier debugging)

## Output Files

VTK files are written alongside CSV and HDF5 outputs:
- **Snapshots**: `snapshot_XXXX.vtk`
- **Checkpoints**: `checkpoint_XXXX.vtk`

## VTK File Structure

The VTK writer outputs UNSTRUCTURED_GRID format with:

### Point Data (Positions)
- 3D coordinates for each particle
- Automatically zero-padded for 1D and 2D simulations

### Scalar Fields
- **density**: Mass density (ρ)
- **pressure**: Hydrodynamic pressure (P)
- **mass**: Particle mass
- **internal_energy**: Specific internal energy (u)
- **smoothing_length**: SPH smoothing length (h)
- **sound_speed**: Local sound speed (cs)
- **particle_id**: Unique particle identifier

### Vector Fields
- **velocity**: 3D velocity vector (v)
- **acceleration**: 3D acceleration vector (a)

## Visualization Workflow

### ParaView
1. Open ParaView
2. File → Open → Select `.vtk` files
3. Click "Apply" in Properties panel
4. Choose coloring variable from toolbar (e.g., "density", "pressure")
5. Use Filters menu for advanced visualization (e.g., Glyph, Slice, Contour)

### Python with PyVista
```python
import pyvista as pv

# Load VTK file
mesh = pv.read('snapshot_0000.vtk')

# Quick plot
mesh.plot(scalars='density', cmap='viridis')

# Extract data
density = mesh['density']
positions = mesh.points
```

## Dimension Handling

The VTK writer automatically handles different spatial dimensions:
- **1D**: Positions output as (x, 0, 0)
- **2D**: Positions output as (x, y, 0)
- **3D**: Positions output as (x, y, z)

Vector fields (velocity, acceleration) follow the same padding scheme.

## File Size Comparison

Typical file sizes for 675 particles:
- **CSV**: ~296 KB (human-readable text)
- **HDF5**: ~119 KB (compressed binary, compression=6)
- **VTK Binary**: ~51 KB (uncompressed binary)

Binary VTK is recommended for production runs. ASCII VTK is useful for debugging.

## Technical Details

### Endianness
Binary VTK files use big-endian byte order per VTK legacy format specification. The writer automatically handles byte-swapping on little-endian systems (x86, ARM).

### Precision
All floating-point data is written as 32-bit floats to match VTK standard practices and reduce file size. Particle IDs are written as 32-bit integers.

### Cell Structure
Particles are represented as individual VERTEX cells in VTK_VERTEX format, allowing direct point-based visualization.

## Implementation

### Source Files
- **Header**: `include/writers/vtk_writer.hpp`
- **Implementation**: `src/writers/vtk_writer.cpp`
- **Integration**: `include/output_manager.hpp`, `src/output_manager.cpp`

### Class Interface
```cpp
class VTKWriter {
public:
    VTKWriter(bool binary = true);
    void open(const std::string& filename, const OutputMetadata& metadata);
    void write_particles(const std::vector<SPHParticle*>& particles);
    void close();
};
```

## Example Configurations

### Production Run (all formats)
```json
"output": {
  "formats": [
    {"type": "csv", "precision": 16},
    {"type": "hdf5", "compression": 6},
    {"type": "vtk", "binary": true}
  ],
  "enableEnergyFile": true
}
```

### Visualization Only
```json
"output": {
  "formats": [
    {"type": "vtk", "binary": true}
  ],
  "enableEnergyFile": false
}
```

### Debugging
```json
"output": {
  "formats": [
    {"type": "csv", "precision": 16},
    {"type": "vtk", "binary": false}
  ],
  "enableEnergyFile": true
}
```

## References

- VTK File Format: https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
- ParaView: https://www.paraview.org/
- VisIt: https://visit-dav.github.io/visit-website/
- PyVista: https://docs.pyvista.org/
