# Relativistic Riemann Solver - Package Summary

## Overview

A complete Python package for solving relativistic Riemann problems, converted from Fortran 77 and modernized with `uv` package management.

## Package Structure

```
relativistic-riemann/
├── pyproject.toml              # Package configuration (uv-based)
├── README.md                   # Comprehensive documentation
├── LICENSE                     # MIT License
├── examples.py                 # Usage examples
├── .venv/                      # Virtual environment (created by uv)
├── src/
│   └── relativistic_riemann/
│       ├── __init__.py        # Package initialization
│       ├── solver.py          # Core Riemann solver
│       ├── test_cases.py      # Predefined test cases
│       ├── cli.py             # Command-line interface
│       └── animate.py         # Visualization tools
└── tests/
    └── test_solver.py         # Unit tests (13 tests, all passing)
```

## Installation

The package is installed using `uv`:

```bash
cd relativistic-riemann
uv venv                        # Create virtual environment
uv pip install -e .           # Install in development mode
```

## Features

### 1. Core Solver (`solver.py`)
- Exact Riemann solver for special relativistic hydrodynamics
- Based on Martí and Müller (J. Fluid Mech., 1994)
- Handles shocks and rarefactions
- Brent's method for pressure iteration
- Newton-Raphson for rarefaction fans

### 2. Predefined Test Cases (`test_cases.py`)
- SR Sod shock tube (γ=1.4)
- Blast wave (γ=5/3, strong pressure ratio)
- Relativistic shock (γ=4/3, v=0.9c)
- Two shock collision (symmetric)

### 3. Command-Line Tools
- `relativistic-riemann`: Interactive solver
- `animate-riemann`: Animation generator

### 4. Visualization (`animate.py`)
- Animated GIFs showing time evolution
- Static multi-time comparison plots
- Four-panel layout (density, pressure, velocity, energy)

### 5. Comprehensive Testing
- 13 unit tests covering:
  - Solver initialization
  - Physical constraints (v < c)
  - Conservation properties
  - Test case definitions
  - All tests passing ✓

## Usage Examples

### As a Library

```python
from relativistic_riemann import RelativisiticRiemannSolver, test_case_sr_sod

# Load test case
test = test_case_sr_sod()

# Create and configure solver
solver = RelativisiticRiemannSolver(gamma=test['gamma'])
solver.set_initial_states(
    test['pl'], test['rhol'], test['vell'],
    test['pr'], test['rhor'], test['velr']
)

# Solve
x, p, rho, vel, u = solver.solve(t=0.4)
```

### Command Line

```bash
# Interactive solving
relativistic-riemann

# Create animation
animate-riemann --test sod --fps 30 --duration 5

# Static comparison
animate-riemann --test relativistic --static
```

## Test Results

All 13 tests pass successfully:

- ✓ Solver initialization
- ✓ Initial state configuration
- ✓ SR Sod shock tube solution
- ✓ Relativistic shock (high velocity)
- ✓ Conservation at t→0
- ✓ Test case definitions
- ✓ Speed of light constraint (v < c)
- ✓ Positive density and pressure

## Dependencies

Managed by `uv` via `pyproject.toml`:
- numpy >= 1.24.0
- matplotlib >= 3.7.0

Development:
- pytest >= 7.0.0
- black, ruff (code formatting/linting)

## Key Improvements from Original Fortran

1. **Modern Python**: Object-oriented design with clear class structure
2. **Package Management**: Professional setup with `uv` and `pyproject.toml`
3. **Documentation**: Comprehensive docstrings and README
4. **Testing**: Full test suite with physical validation
5. **Visualization**: Built-in animation and plotting tools
6. **Command-Line Interface**: Easy-to-use CLI tools
7. **Reusability**: Can be used as library or standalone

## Generated Files

The package has been tested and generates:
- `solution.dat`: Numerical solution files
- `*.gif`: Animated visualizations
- `*.png`: Static comparison plots
- Example outputs in `examples.py`

## Verification

✓ Package builds successfully with `uv`
✓ All 13 unit tests pass
✓ CLI commands work correctly
✓ Animations generate successfully
✓ Example scripts run without errors

## Location

`/Users/guo/Downloads/sphcode/relativistic-riemann/`

The package is ready for:
- Local use and development
- Publishing to PyPI (if desired)
- Integration into other projects
- Further extension and customization
