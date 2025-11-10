# Code Style and Conventions

## Naming Conventions

### Types and Classes
- **Classes/Structs**: PascalCase
  - Examples: `Solver`, `SPHParticle`, `BHTree`, `CheckpointMetadata`
  
### Functions and Methods
- **Public methods**: camelCase (mostly)
  - Examples: `readParameterfile()`, `makeInitialCondition()`, `integrate()`
- **Some use snake_case**: `make_shock_tube()`, `output_particle()`
- **Prefer camelCase for new code**

### Variables
- **Local variables**: snake_case
  - Examples: `particle_count`, `time_step`, `density`
- **Member variables**: `m_` prefix + snake_case
  - Examples: `m_param`, `m_output`, `m_sim`, `m_timestep`
- **Global constants**: SCREAMING_SNAKE_CASE or kCamelCase

### Files
- **Headers**: snake_case.hpp
  - Examples: `solver.hpp`, `bhtree.hpp`, `fluid_force.hpp`
- **Source**: snake_case.cpp
  - Examples: `solver.cpp`, `checkpoint.cpp`, `output.cpp`

## C++ Standards

### Language Version
- **C++14** standard
- Modern C++ features encouraged:
  - `std::unique_ptr`, `std::shared_ptr`
  - `auto` keyword
  - Range-based for loops
  - `constexpr` for constants
  - `enum class` over plain enums

### Memory Management
- **RAII principle**: Resources owned by objects
- **Smart pointers**: Use `std::shared_ptr`, `std::unique_ptr`
- **Avoid**: Raw `new`/`delete`, manual memory management
- **Current usage**: Some raw pointers exist (legacy)

### Type Aliases
- `real`: typedef for float/double (defined in defines.hpp)
- `vec_t`: typedef for vector type (3D array of real)

## Code Organization

### Header Files
- Include guards: `#ifndef SPH_FILENAME_HPP` / `#define SPH_FILENAME_HPP`
- Namespace: `namespace sph { ... }`
- Forward declarations before includes when possible
- Minimal includes in headers

### Implementation Files
- Include corresponding header first
- Group includes: STL, external libs, project headers
- Keep functions short and focused
- One class per file (generally)

## Comments and Documentation

### Style
- **Explain WHY, not WHAT**: Comment reasoning, not obvious code
- **Important invariants**: Document assumptions and constraints
- **Algorithms**: Reference papers/equations where applicable
- **TODO markers**: Use `// TODO:` for future work

### Format
```cpp
// Single-line comments for brief explanations

/*
 * Multi-line comments for detailed explanations,
 * algorithm descriptions, or complex logic
 */
```

## Macros

### Current Usage
- `DIM`: Dimensionality (1D, 2D, 3D) - set in defines.hpp
- Include guards
- Debug flags: `_DEBUG`

### New Code Rules (per coding_rule.instructions.md)
- **Avoid macros**: Use constexpr, inline functions, templates instead
- **If needed**: SCREAMING_SNAKE_CASE with SPH_ prefix
- **Document**: Why macro is necessary vs language feature
- **Parenthesize**: All parameters and expansion
- **No side effects**: Keep macros pure

## Floating Point

### Best Practices
- **No direct equality**: Never use `x == y` for floats
- **Tolerance comparisons**: Use epsilon-based comparisons
- **Prefer double**: Unless memory/performance critical
- **Guard NaN/Inf**: Use `std::isfinite()` for validation

## Constants

### Definition
```cpp
// Prefer constexpr over #define
constexpr real PI = 3.14159265358979323846;
constexpr int DEFAULT_NEIGHBORS = 32;

// Or enum class
enum class Kernel { CubicSpline, Wendland };
```

### No Magic Numbers
- Define named constants for all numerical values
- Document physical meaning/units

## Error Handling

### Current Approach
- Some use of exceptions (see exception.hpp)
- Logging via Logger class
- Return codes in some places

### Preferred
- Exceptions for errors
- Logging for diagnostics
- `std::optional` for optional results
- Clear error messages

## Testing

### Philosophy (per coding_rule.instructions.md)
- **Test-Driven Development**: Write tests first
- **Behavior-Driven Style**: Given/When/Then structure
- **Coverage**: Boundary conditions, error cases, invariants

### Current State
- Limited formal test suite
- Manual testing via sample cases
- Need to add comprehensive unit tests

## Thread Safety

### OpenMP Usage
- Parallel regions marked with `#pragma omp parallel`
- Critical sections protected
- Race conditions avoided
- Reduction operations for accumulation

### Best Practices
- Minimize shared mutable state
- Use thread-local storage when possible
- Document thread-safety assumptions

## Repository Hygiene

### Keep Clean
- No `.bak`, `.orig`, `.tmp` files committed
- No editor swap files (`.swp`, etc.)
- No unnecessary logs
- Update .gitignore for artifacts

### Before Commit
- Remove temporary files
- Verify no backup files added
- Clean build artifacts not in .gitignore
