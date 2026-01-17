// Morton code particle reordering for cache-friendly tree traversal
//
// This file instantiates the Morton code functions from the header.
// Most functions are inline, but this ensures proper compilation and
// provides a place for any future non-inline implementations.

#include "morton.hpp"

namespace sph {
namespace morton {

// Explicit instantiation check - this ensures the header compiles correctly
// The actual implementations are all inline in morton.hpp

// Future: Add any heavy/non-inline implementations here if needed
// For example, radix sort for very large particle counts could go here

} // namespace morton
} // namespace sph
