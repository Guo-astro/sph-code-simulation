# Code Structure

## Directory Layout

```
sphcode/
├── include/              # Header files
│   ├── disph/           # Density Independent SPH headers
│   ├── gsph/            # Godunov SPH headers
│   ├── kernel/          # Kernel function headers
│   ├── relaxation/      # Lane-Emden relaxation headers
│   ├── bhtree.hpp       # Barnes-Hut tree
│   ├── checkpoint.hpp   # Checkpoint system (to be replaced)
│   ├── defines.hpp      # Global definitions and types
│   ├── output.hpp       # Output system (to be replaced)
│   ├── solver.hpp       # Main solver class
│   ├── simulation.hpp   # Simulation state
│   ├── particle.hpp     # SPH particle structure
│   └── parameters.hpp   # Simulation parameters
│
├── src/                 # Implementation files
│   ├── disph/          # DISPH implementations
│   ├── gsph/           # GSPH implementations
│   ├── relaxation/     # Lane-Emden relaxation code
│   ├── sample/         # Sample initial conditions
│   ├── solver.cpp      # Solver implementation
│   ├── output.cpp      # Output implementation
│   └── checkpoint.cpp  # Checkpoint implementation
│
├── sample/             # Example configurations
│   ├── lane_emden/     # Lane-Emden polytrope examples
│   ├── shock_tube/     # 1D shock tube
│   ├── khi/            # Kelvin-Helmholtz instability
│   ├── evrard/         # Evrard collapse
│   └── ...
│
├── lane_emden/         # Lane-Emden specific files
│   ├── config/         # Configuration templates
│   ├── results/        # Output results
│   ├── scripts/        # Python analysis scripts
│   └── docs/           # Lane-Emden documentation
│
├── docs/               # General documentation
├── test/               # Test files
├── build/              # CMake build directory
└── CMakeLists.txt      # CMake configuration
```

## Key Classes

### Core Classes
- `Solver`: Main simulation driver
- `Simulation`: Holds particle data and simulation state
- `SPHParticle`: Individual particle structure (20 fields)
- `SPHParameters`: Simulation parameters from JSON config
- `Output`: Current output system (to be replaced)
- `Checkpoint`: Current checkpoint system (to be replaced)

### SPH Methods
- `StandardSPH`: Density-energy formulation
- `DensityIndependentSPH`: Pressure-energy formulation
- `GodunovSPH`: Riemann solver based

### Utilities
- `BHTree`: Barnes-Hut tree for gravity and neighbor search
- `ExhaustiveSearch`: Brute-force neighbor search
- `Module`: Physics module base class
- `Logger`: Logging utilities

## File Naming Convention
- Headers: `snake_case.hpp`
- Source: `snake_case.cpp`
- Classes: PascalCase
- Functions/methods: camelCase (mostly, some snake_case)
- Variables: snake_case or m_ prefix for members
