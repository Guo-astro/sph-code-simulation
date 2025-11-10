# Sample Test Cases

## Available Samples

### 3D Gravity Tests
- **evrard** - Standard Evrard collapse (equilibrium after shock)
- **evrard_bound** - Bound Evrard (multiple collapse/rebound cycles)
- **lane_emden** - Lane-Emden hydrostatic sphere

### 2D Hydrodynamic Tests
- **khi** - Kelvin-Helmholtz instability
- **gresho_chan_vortex** - Gresho-Chan vortex
- **hydrostatic** - Hydrostatic equilibrium
- **pairing_instability** - Pairing instability

### 1D Tests
- **shock_tube** - Sod shock tube

## Running Samples

### Using Makefile.automation:
```bash
make -f Makefile.automation evrard
make -f Makefile.automation evrard_bound
make -f Makefile.automation lane_emden
```

### Manual:
```bash
./build/sph evrard
```

Results are saved to `sample/<name>/results/`
