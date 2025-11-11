# Output Configuration Examples

## Design Philosophy

**Each format is a self-contained object with its own configuration.**

The output configuration uses an array of format objects where each object specifies the format type and its specific options. This keeps related settings together and makes the configuration cleaner and more maintainable.

## Example Configurations

### CSV Only (Human-Readable Text)
```json
{
  "output": {
    "formats": [
      {
        "type": "csv",
        "precision": 16
      }
    ],
    "enableEnergyFile": true
  }
}
```

### HDF5 Only (Compressed Binary)
```json
{
  "output": {
    "formats": [
      {
        "type": "hdf5",
        "compression": 6
      }
    ],
    "enableEnergyFile": true
  }
}
```

### CSV + HDF5 (Both Formats)
```json
{
  "output": {
    "formats": [
      {
        "type": "csv",
        "precision": 16
      },
      {
        "type": "hdf5",
        "compression": 6
      }
    ],
    "enableEnergyFile": true
  }
}
```

### VTK Only (For Visualization)
```json
{
  "output": {
    "formats": [
      {
        "type": "vtk",
        "binary": true
      }
    ],
    "enableEnergyFile": true
  }
}
```

### All Formats
```json
{
  "output": {
    "formats": [
      {
        "type": "csv",
        "precision": 16
      },
      {
        "type": "hdf5",
        "compression": 6
      },
      {
        "type": "vtk",
        "binary": true
      }
    ],
    "enableEnergyFile": true
  }
}
```

## Configuration Options

### Format Objects

Each format object contains:
- **`type`** (string, required): Format name - `"csv"`, `"hdf5"`, or `"vtk"`
- Format-specific options (see below)

### CSV Format Options
- **`precision`** (integer, default: 16): Decimal precision for floating-point values
  - Example: `{"type": "csv", "precision": 16}`

### HDF5 Format Options
- **`compression`** (integer, default: 6): Compression level (0-9)
  - 0 = no compression
  - 6 = balanced (recommended)
  - 9 = maximum compression (slower)
  - Example: `{"type": "hdf5", "compression": 6}`

### VTK Format Options
- **`binary`** (boolean, default: true): Use binary format vs ASCII
  - `true` = smaller files, faster I/O (recommended)
  - `false` = human-readable ASCII
  - Example: `{"type": "vtk", "binary": true}`

### General Options
- **`enableEnergyFile`** (boolean, default: true): Write energy.dat file with global diagnostics

## Legacy Formats (Deprecated but Still Supported)

### Legacy Format 1: Simple String Array
```json
{
  "output": {
    "formats": ["csv", "hdf5"],
    "csvPrecision": 16,
    "hdf5Compression": 6,
    "enableEnergyFile": true
  }
}
```

### Legacy Format 2: Boolean Flags
```json
{
  "output": {
    "enableCSV": true,
    "enableHDF5": false,
    "enableVTK": false,
    "csvPrecision": 16,
    "enableEnergyFile": true
  }
}
```

**Recommendation**: Migrate to the new object-based format for cleaner, self-contained configs.

## Best Practices

1. **Development/Debugging**: Use CSV only for easy inspection
   ```json
   "formats": [{"type": "csv", "precision": 16}]
   ```

2. **Production Runs**: Use HDF5 for smaller files and faster I/O
   ```json
   "formats": [{"type": "hdf5", "compression": 6}]
   ```

3. **Visualization Workflows**: Add VTK when you need ParaView/VisIt
   ```json
   "formats": [
     {"type": "hdf5", "compression": 6},
     {"type": "vtk", "binary": true}
   ]
   ```

4. **Redundancy/Safety**: Use both CSV and HDF5 if you want a human-readable backup
   ```json
   "formats": [
     {"type": "csv", "precision": 16},
     {"type": "hdf5", "compression": 6}
   ]
   ```

## Advantages of Object-Based Format

✅ **Self-contained**: Each format's options are clearly grouped with the format  
✅ **No clutter**: Only include formats you actually use  
✅ **Type-safe**: Each format object is validated independently  
✅ **Extensible**: Easy to add new format-specific options  
✅ **Clear**: Immediately obvious which option belongs to which format  

## File Size Comparison (typical 100k particles)

| Format | Size | Speed | Human Readable |
|--------|------|-------|----------------|
| CSV    | ~50 MB | Slow | ✅ Yes |
| HDF5 (compressed) | ~5 MB | Fast | ❌ No |
| VTK Binary | ~30 MB | Medium | ❌ No |
| VTK ASCII | ~55 MB | Slow | ✅ Yes |
