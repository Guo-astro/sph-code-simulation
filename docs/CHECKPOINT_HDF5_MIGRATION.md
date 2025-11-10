# Checkpoint System Migration to HDF5

## Summary
The checkpoint system has been migrated from the legacy binary `.chk` format to the new unified HDF5-based output system. This ensures consistency across all output formats (snapshots, checkpoints, energy logs) and provides better compression and metadata support.

## Changes Made

### 1. Removed Legacy Checkpoint Class Usage
- **Removed** `m_checkpoint` member from `Solver` class
- **Kept** `m_checkpoint_config` for configuration
- All checkpoint operations now go through `OutputManager`

### 2. Updated OutputManager

#### Constructor Changes
```cpp
// Old
OutputManager(const OutputConfig& config, const UnitSystem& units, const std::string& output_dir);

// New
OutputManager(const OutputConfig& config, const UnitSystem& units, 
              const std::string& output_dir, const std::string& checkpoint_dir = "");
```

The checkpoint directory parameter defaults to `output_dir/checkpoints` if empty.

#### New Methods
- `generate_checkpoint_filename()` - Generates checkpoint filenames in checkpoint directory
- Checkpoint directory creation in `initialize()`

### 3. Checkpoint File Format

#### Old Format (REMOVED)
- Binary `.chk` files
- Custom serialization
- Limited metadata
- No compression
- Example: `checkpoint_095000.chk`

#### New Format (CURRENT)
- HDF5 `.h5` files (compressed binary)
- CSV `.csv` files (human-readable fallback)
- Full metadata in HDF5 attributes
- HDF5 compression (6-9x reduction)
- Example: `checkpoint_095000.h5`

### 4. File Sizes
For 5400 particles:
- **CSV**: 2.3 MB (human-readable, metadata in header)
- **HDF5**: 535 KB (compressed binary, ~77% smaller)
- **Old .chk**: ~500 KB (binary, no metadata, DEPRECATED)

### 5. HDF5 Checkpoint Structure
```
checkpoint_NNNNNN.h5
├── /metadata (group)
│   ├── format_version (attribute)
│   ├── particle_count (attribute)
│   ├── time_code (attribute)
│   ├── time_physical (attribute)
│   ├── step (attribute)
│   ├── gamma (attribute)
│   ├── kernel_type (attribute)
│   ├── metadata_json (attribute - full metadata)
│   └── checkpoint_json (attribute - checkpoint-specific)
├── /particles (group)
│   ├── x, y, z (datasets)
│   ├── vx, vy, vz (datasets)
│   ├── mass, rho, pres, ene (datasets)
│   ├── h (smoothing length)
│   └── ... (all particle properties)
└── /energy (group)
    ├── kinetic (dataset)
    ├── thermal (dataset)
    ├── potential (dataset)
    └── total (dataset)
```

## Configuration Changes

### Old Configuration (INVALID)
```json
{
  "checkpoint": {
    "resumeFile": "path/to/checkpoint_095000.chk"  // ❌ Will not work
  }
}
```

### New Configuration (CORRECT)
```json
{
  "checkpoint": {
    "enabled": true,
    "saveInterval": 5000,
    "directory": "lane_emden/results/polytrope_n1.5_3d/checkpoints",
    "saveOnExit": true,
    "maxCheckpoints": 5,
    "autoResume": true,
    "resumeFile": "lane_emden/results/polytrope_n1.5_3d/checkpoints/checkpoint_095000.h5"  // ✅ HDF5 format
  }
}
```

## Code Changes

### Solver.cpp
```cpp
// Old
m_checkpoint->save(m_sim, chk_meta);

// New
m_output_manager->write_checkpoint(m_sim, m_param, chk_meta, loop);
```

```cpp
// Old
m_checkpoint->load(resume_file, m_sim, resume_metadata);

// New
m_output_manager->load_for_resume(resume_file, m_sim, resume_metadata);
```

```cpp
// Old
m_output_manager = std::make_shared<OutputManager>(output_config, m_units, m_output_dir);

// New
std::string checkpoint_dir = m_checkpoint_config.directory;
m_output_manager = std::make_shared<OutputManager>(output_config, m_units, m_output_dir, checkpoint_dir);
```

## Benefits

1. **Consistency**: All outputs use the same unified system
2. **Compression**: 77% file size reduction vs CSV, similar to old .chk but with metadata
3. **Metadata**: Full simulation state and parameters stored in HDF5 attributes
4. **Portability**: HDF5 is cross-platform and well-supported
5. **Resumability**: Complete state preserved for exact resumption
6. **Debugging**: CSV checkpoints available as human-readable fallback

## Testing

### Test Resume from Checkpoint
```bash
# Create checkpoint
./build/sph lane_emden

# Resume from checkpoint (edit config to add resumeFile)
./build/sph sample/lane_emden/lane_emden_resume_test.json
```

### Verify Checkpoint Contents
```bash
# View HDF5 structure
h5dump -H lane_emden/results/polytrope_n1.5_3d/checkpoints/checkpoint_0000.h5

# View metadata
h5dump -A lane_emden/results/polytrope_n1.5_3d/checkpoints/checkpoint_0000.h5
```

## Migration Notes

1. **Old .chk files cannot be read** - Must regenerate checkpoints
2. **Configuration must specify .h5 extension** for resumeFile
3. **Checkpoint directory** is created separately from output directory
4. **Both CSV and HDF5** checkpoints are created if both are enabled
5. **HDF5 is primary** - Resume uses HDF5 format for efficiency

## Backward Compatibility

⚠️ **BREAKING CHANGE**: Old `.chk` checkpoint files are no longer supported. Simulations must be restarted or checkpoints regenerated.

To migrate:
1. Update configuration to use `.h5` extension for `resumeFile`
2. Remove old `.chk` checkpoint references
3. Let simulation create new HDF5 checkpoints
4. Verify resume functionality with new format

## Implementation Status

✅ **Complete**
- OutputManager accepts checkpoint directory
- Checkpoint files written to separate directory
- HDF5 checkpoint write/read working
- Resume from HDF5 checkpoint validated
- CSV checkpoint fallback functional
- Metadata preservation confirmed
- File compression verified
