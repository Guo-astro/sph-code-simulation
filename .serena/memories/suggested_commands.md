# Suggested Commands

## Build Commands

### Using CMake (Recommended)
```bash
# Initial build
mkdir -p build
cd build
cmake ..
make -j8

# Rebuild
cd build
make -j8

# Clean build
cd build
make clean
make -j8

# Full rebuild
rm -rf build
mkdir build
cd build
cmake ..
make -j8
```

### Using old Makefile (Deprecated)
The old Makefile now redirects to CMake. Use CMake instead.

## Running Simulations

### Sample Simulations
```bash
# Run with sample name (auto-loads from sample/<name>/<name>.json)
./build/sph <sample_name> [num_threads]

# Examples:
./build/sph shock_tube       # 1D shock tube
./build/sph khi             # Kelvin-Helmholtz instability
./build/sph evrard          # Evrard collapse
./build/sph lane_emden      # Lane-Emden polytrope relaxation

# Specify threads (optional, defaults to max available)
./build/sph shock_tube 4
```

### Custom Configuration
```bash
# Run with custom JSON config file
./build/sph path/to/config.json [num_threads]
```

## Lane-Emden Specific Commands

### Using Makefile.lane_emden presets
```bash
# Show help
make lane_emden_help

# List available presets
make lane_emden_list

# Run with preset
make lane_emden_run_<preset>

# Resume from checkpoint
make lane_emden_resume_<preset>

# Animate results
make lane_emden_animate
```

## Testing
```bash
# Currently no formal test suite
# Test by running sample cases and verifying output
./build/sph shock_tube
./build/sph khi
```

## Linting/Formatting
```bash
# No automated linting/formatting configured yet
# Follow coding_rule.instructions.md manually
```

## Git Workflow
```bash
git status
git add <files>
git commit -m "message"
git push
```

## macOS Specific (Darwin)

### Find files
```bash
find . -name "*.cpp"
find . -type f -name "*.hpp"
```

### List directory
```bash
ls -lh
ls -lha  # include hidden
```

### Grep search
```bash
grep -r "pattern" src/
grep -n "pattern" file.cpp  # with line numbers
```

### Process management
```bash
top          # view processes
ps aux       # list all processes
```

## Analysis/Visualization

### Python scripts (Lane-Emden)
```bash
# Located in lane_emden/scripts/
python3 lane_emden/scripts/animate_relaxation.py
python3 lane_emden/scripts/analyze_profile.py
```

### View logs
```bash
tail -f relaxation.log       # follow log in real-time
cat relaxation_debug.log     # view full log
```
