# Task Completion Checklist

When a coding task is complete, ensure the following:

## 1. Code Quality

### Build Successfully
```bash
cd build
cmake ..
make -j8
```
- ✅ No compilation errors
- ✅ No warnings (treat warnings as errors: `-Wall -Wextra`)
- ✅ Clean build with optimization flags

### Code Standards
- ✅ Follows naming conventions (PascalCase types, camelCase methods, snake_case variables)
- ✅ No hard-coded strings or magic numbers (use constexpr constants)
- ✅ No raw new/delete (use smart pointers)
- ✅ Proper RAII for resource management
- ✅ Thread-safe for OpenMP (if applicable)

### Documentation
- ✅ Comments explain WHY, not WHAT
- ✅ Important invariants documented
- ✅ Algorithm references included (papers/equations)
- ✅ Root-cause noted for bug fixes
- ✅ Change log entry (if maintaining changelog)

## 2. Testing

### Unit Tests
- ✅ Tests written FIRST (TDD approach)
- ✅ Given/When/Then structure
- ✅ Cover boundary conditions
- ✅ Cover error cases
- ✅ Verify invariants
- ✅ All tests pass

### Integration Tests
```bash
# Run sample cases to verify
./build/sph shock_tube
./build/sph khi
./build/sph lane_emden
```
- ✅ Sample simulations run successfully
- ✅ Output files generated correctly
- ✅ Numerical results reasonable

### Regression Tests
- ✅ Existing functionality not broken
- ✅ Performance not degraded
- ✅ Backward compatibility (if required)

## 3. File Structure

### Correct Placement (Verify with Serena MCP)
- ✅ Headers in correct include/ subdirectory
- ✅ Source in correct src/ subdirectory
- ✅ Tests in test/ directory
- ✅ Configs in appropriate sample/ or config/ directory
- ✅ Filenames match primary type/purpose

### CMakeLists.txt Updated
```cmake
# If new files added:
target_sources(sph PRIVATE
  src/new_file.cpp
)
```
- ✅ New source files added to build
- ✅ New dependencies linked
- ✅ Include directories updated

## 4. Repository Hygiene

### Clean Workspace
```bash
# Check for artifacts
find . -name "*.bak" -o -name "*.orig" -o -name "*.tmp"
git status
```
- ✅ No `.bak`, `.orig`, `.tmp` files
- ✅ No editor swap files (`.swp`, `~` files)
- ✅ No unnecessary log files committed
- ✅ No commented-out code blocks
- ✅ .gitignore updated for new artifacts

### Legacy Code Removal
- ✅ Old code DELETED (not just commented out)
- ✅ No backward compatibility layers introduced
- ✅ Migration done in same change (if practical)
- ✅ No orphaned files left behind

## 5. Design Principles

### Architecture (per coding_rule.instructions.md)
- ✅ Simplest correct design addressing root cause
- ✅ No new compatibility layers
- ✅ Modular with single responsibility
- ✅ Reused existing utilities (checked with Serena MCP)
- ✅ No code repetition (extracted to functions/templates)

### Modern C++ Practices
- ✅ RAII for all resources
- ✅ Smart pointers over raw pointers
- ✅ `constexpr` over `#define`
- ✅ `enum class` over plain enum
- ✅ `std::optional`, `std::variant` where appropriate
- ✅ No implicit conversions
- ✅ Functions small and pure (minimal side effects)

### Macros (if used)
- ✅ Only when NO safer language feature exists
- ✅ SCREAMING_SNAKE_CASE with SPH_ prefix
- ✅ Fully parenthesized
- ✅ Documented with purpose and example
- ✅ No side effects
- ✅ Removal task created (if temporary)

### Floating Point
- ✅ No direct equality comparisons
- ✅ Tolerance-based comparisons
- ✅ `std::isfinite()` guards for NaN/Inf
- ✅ `double` preferred (unless justified)

## 6. Documentation

### User-Facing
- ✅ README.md updated (if public API changed)
- ✅ User guide updated (if applicable)
- ✅ Example configs updated
- ✅ Migration guide (if breaking change)

### Developer-Facing
- ✅ Design documents updated
- ✅ API reference updated
- ✅ Implementation notes added
- ✅ Serena MCP memories updated (if structure changed)

## 7. Version Control

### Git Commit
```bash
git add <files>
git status  # verify only intended files staged
git commit -m "Concise message"
```
- ✅ Only relevant files staged
- ✅ Commit message clear and concise
- ✅ Atomic commits (one logical change per commit)
- ✅ No merge conflicts

### Before Push
- ✅ Local build passes
- ✅ All tests pass locally
- ✅ Code reviewed (self-review minimum)
- ✅ Ready for CI/CD (if configured)

## 8. Special Considerations

### For Output System Redesign
- ✅ Old checkpoint.hpp/cpp DELETED
- ✅ Old output.hpp/cpp DELETED
- ✅ No .chk or .dat files in repo
- ✅ All configs use new output format
- ✅ HDF5 dependency verified
- ✅ Unit system tested (conversions correct)
- ✅ CSV format self-documenting
- ✅ HDF5 compression working
- ✅ Resume from checkpoint tested

### For Lane-Emden Code
- ✅ Relaxation step counter accurate
- ✅ Alpha scaling correct
- ✅ Configuration hash matches
- ✅ Preset system working
- ✅ Python scripts compatible

## Final Verification

```bash
# Clean rebuild
rm -rf build
mkdir build
cd build
cmake ..
make -j8

# Run tests
./build/sph shock_tube
./build/sph lane_emden

# Check workspace
cd ..
find . -name "*.bak" -o -name "*.orig"
git status
```

When ALL checkboxes ✅, the task is COMPLETE.
