# Output & Unit System Redesign - Quick Reference

## 📋 What's Changing

### Before (Current System)
- ❌ ASCII `.dat` files (incomplete particle state)
- ❌ Binary `.chk` files (custom format, non-standard)
- ❌ No unit system support (code units only)
- ❌ Separate output and checkpoint files
- ❌ No metadata

### After (New System)
- ✅ **CSV format** (human-readable, complete, self-documenting)
- ✅ **HDF5 format** (compressed, fast, Python-ready)
- ✅ **Unit system** (Galactic/SI/CGS with automatic conversion)
- ✅ **Unified format** (same file for visualization AND resume)
- ✅ **Rich metadata** (units, parameters, provenance)

---

## 🎯 Key Features

### 1. Multiple Output Formats

**CSV** - Human-readable with headers:
```
# SPH Simulation Output - CSV Format v1.0
# Unit System: Galactic (kpc, M_sun, km/s)
# Step: 95000, Time: 112.045 [code] = 1.234 Myr
id,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,...
0,0.123,0.456,0.789,...
```

**HDF5** - Compressed binary (4x smaller):
```
/metadata/ → simulation info, units, parameters
/particles/ → pos[N,3], vel[N,3], mass[N], ... (compressed)
/energy/ → kinetic, thermal, potential
```

### 2. Flexible Unit Systems

**Configuration**:
```json
{
  "units": {
    "system": "galactic",  // or "si", "cgs", "code"
    "code_length": 1.0,    // in kpc
    "code_mass": 1.0,      // in M_sun
    "code_velocity": 1.0   // in km/s
  }
}
```

**Supported Systems**:
- **Galactic**: kpc, M☉, km/s (astrophysics)
- **SI**: m, kg, s (engineering)
- **CGS**: cm, g, s (physics)
- **Code**: dimensionless (internal)

### 3. Resume from Any Output

**Old way** (separate formats):
```
output_00095.dat  ← visualization only
checkpoint_095000.chk  ← resume only
```

**New way** (unified):
```
snapshot_00095.h5  ← both visualization AND resume
```

Simply set in config:
```json
{
  "checkpoint": {
    "resumeFile": "results/snapshot_00095.h5"
  }
}
```

---

## 📊 File Comparison

| Feature | Old .dat | Old .chk | New CSV | New HDF5 |
|---------|----------|----------|---------|----------|
| Size (5400 particles) | 723 KB | 738 KB | 850 KB | **180 KB** |
| Human readable | Partial | No | ✅ Yes | No |
| Complete state | ❌ No | ✅ Yes | ✅ Yes | ✅ Yes |
| Resume capable | ❌ No | ✅ Yes | ✅ Yes | ✅ Yes |
| Metadata | ❌ No | Partial | ✅ Full | ✅ Full |
| Unit info | ❌ No | ❌ No | ✅ Yes | ✅ Yes |
| Python access | Manual | Manual | pandas | **h5py** |
| Standard format | ❌ No | ❌ No | ✅ Yes | ✅ Yes |

---

## 🐍 Python Usage

### Reading CSV
```python
import pandas as pd

df = pd.read_csv('snapshot_00095.csv', comment='#')
print(df.head())
print(f"Particles: {len(df)}")

# Plot
import matplotlib.pyplot as plt
plt.scatter(df['pos_x'], df['pos_y'], s=1)
plt.show()
```

### Reading HDF5
```python
import h5py
import numpy as np

with h5py.File('snapshot_00095.h5', 'r') as f:
    # Metadata
    time = f['metadata'].attrs['time_physical']
    units = f['metadata/unit_system'].attrs['type_name']
    print(f"Time: {time:.3f} Myr, Units: {units}")
    
    # Particle data (fast!)
    pos = f['particles/pos'][:]  # [N, 3] array
    vel = f['particles/vel'][:]
    mass = f['particles/mass'][:]
    
    # Analysis
    com = np.sum(pos * mass[:, np.newaxis], axis=0) / np.sum(mass)
    print(f"Center of mass: {com}")
```

### yt (astrophysics analysis)
```python
import yt

ds = yt.load('snapshot_00095.h5')
yt.ProjectionPlot(ds, 'z', 'density').save()
```

---

## 🔧 Configuration Examples

### Minimal (defaults)
```json
{
  "outputDirectory": "results/my_simulation",
  "endTime": 5.0
}
```
→ Outputs in code units, both CSV and HDF5

### Galactic Units
```json
{
  "units": {
    "system": "galactic",
    "code_length": 10.0,   // 10 kpc
    "code_mass": 1e6,      // 10^6 M_sun
    "code_velocity": 1.0   // 1 km/s
  },
  "output": {
    "formats": ["hdf5"],   // HDF5 only for speed
    "hdf5": {
      "compression": 9     // maximum compression
    }
  }
}
```

### High-Resolution CSV
```json
{
  "output": {
    "formats": ["csv"],
    "csv": {
      "precision": 20  // extra precision for debugging
    }
  }
}
```

---

## 🚀 Migration Guide

### Important: Clean Break - No Backward Compatibility

**This is a complete redesign.** Old `.dat` and `.chk` files will NOT be supported. You'll start fresh with the new format.

### Step 1: Backup Your Old Data (Optional)
```bash
# If you want to keep old results for reference
mkdir old_format_backup
mv lane_emden/results/*.dat lane_emden/results/*.chk old_format_backup/
```

### Step 2: Update Configuration
Add to your JSON config:
```json
{
  "units": {
    "system": "code"  // or "galactic", "si", "cgs"
  },
  "output": {
    "formats": ["csv", "hdf5"]
  }
}
```

**Remove** old `checkpoint` sections - they're no longer used.

### Step 3: Recompile
```bash
cd build
cmake ..
make -j8
```

### Step 4: Run Fresh Simulation
```bash
./build/sph lane_emden
```

Output files:
```
results/
├── snapshot_00000.csv
├── snapshot_00000.h5
├── snapshot_00001.csv
├── snapshot_00001.h5
└── energy.dat
```

### Step 5: Resume from HDF5
To resume, just point to any HDF5 snapshot:
```json
{
  "output": {
    "resumeFile": "results/snapshot_00095.h5"
  }
}
```

**No conversion tools needed** - start fresh with better formats!

---

## ⏱️ Implementation Timeline

| Phase | Days | Tasks |
|-------|------|-------|
| **Phase 1** | 1 | Foundation: HDF5/JSON deps, UnitSystem class |
| **Phase 2** | 1 | CSV Output: CSVWriter with metadata headers |
| **Phase 3** | 2 | HDF5 Output: HDF5Writer with compression |
| **Phase 4** | 1 | Integration: OutputManager, remove legacy |
| **Phase 5** | 1 | Documentation: User guide, cleanup |
| **Phase 6** | 1 | Testing: Integration tests, benchmarks |

**Total**: ~7 days (1 developer) or ~3 days (2 developers in parallel)

**Key Milestone**: Phase 4 completes legacy removal - clean break from old formats!

---

## 📚 Documentation

### For Users
1. **USER_GUIDE_OUTPUT.md** - How to use new output system
2. **MIGRATION_GUIDE.md** - Migrating from old format
3. **README.md** - Updated with examples

### For Developers
1. **OUTPUT_UNIT_SYSTEM_DESIGN.md** - Full design document
2. **IMPLEMENTATION_PLAN.md** - Step-by-step tasks
3. **API_REFERENCE.md** - Code documentation

---

## ✅ Benefits

### Performance
- **4x smaller files** with HDF5 compression
- **Faster I/O** with binary format
- **Parallel-ready** for future MPI support

### Usability
- **Self-documenting** outputs with full metadata
- **Python-friendly** (pandas, h5py, yt)
- **Standard formats** (CSV, HDF5)

### Scientific
- **Unit conversion** built-in
- **Complete provenance** tracking
- **Reproducibility** with embedded parameters

### Maintenance
- **Single resume format** (no separate .chk files)
- **Extensible** architecture (easy to add formats)
- **Well-tested** with comprehensive test suite

---

## 🆘 Getting Help

### Common Issues

**Q: HDF5 not found during build?**
```bash
# macOS
brew install hdf5

# Ubuntu
sudo apt-get install libhdf5-dev
```

**Q: Can't read old .chk files?**
**A:** Not supported. The new system is a clean break from legacy formats. Start fresh with HDF5 - it's better in every way!

**Q: I have important data in old .dat files!**
**A:** You can manually extract key data if needed, but the old format was incomplete (missing sound speed, balsara factor, gravitational potential). For production work, use the new HDF5 format going forward.

**Q: How to choose between CSV and HDF5?**
- **CSV**: Debugging, small runs, need human readability
- **HDF5**: Production, large runs, Python analysis, resuming

**Q: Which unit system should I use?**
- **Galactic**: Galaxy simulations, stellar systems
- **SI**: Engineering applications
- **CGS**: Physics papers
- **Code**: Code development, unit testing

---

## 🎓 Learning Resources

### HDF5
- [HDF5 Tutorial](https://portal.hdfgroup.org/display/HDF5/Learning+HDF5)
- [h5py Documentation](https://docs.h5py.org/)

### Unit Systems
- [Galactic Units](https://en.wikipedia.org/wiki/Astronomical_unit#Natural_units)
- [CGS Units](https://en.wikipedia.org/wiki/Centimetre–gram–second_system_of_units)

### SPH Output Analysis
- [yt Project](https://yt-project.org/)
- [pynbody](https://pynbody.github.io/pynbody/)

---

*For detailed information, see:*
- `docs/OUTPUT_UNIT_SYSTEM_DESIGN.md` - Complete design
- `docs/IMPLEMENTATION_PLAN.md` - Task breakdown
