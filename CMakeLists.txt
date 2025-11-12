cmake_minimum_required(VERSION 3.13)
project(sphcode CXX)

# Dimension parameter (1D, 2D, or 3D)
# Can be set via: cmake -DSPH_DIM=1 ..
set(SPH_DIM "2" CACHE STRING "Spatial dimension (1, 2, or 3)")
set_property(CACHE SPH_DIM PROPERTY STRINGS "1" "2" "3")
message(STATUS "Building for ${SPH_DIM}D simulations (DIM=${SPH_DIM})")

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)

# Find HDF5 for new output system
find_package(HDF5 REQUIRED COMPONENTS C CXX)
message(STATUS "Found HDF5: ${HDF5_VERSION}")

# Find or fetch nlohmann/json for JSON handling
find_package(nlohmann_json 3.11.0 QUIET)
if(NOT nlohmann_json_FOUND)
    message(STATUS "nlohmann/json not found, fetching from GitHub...")
    include(FetchContent)
    FetchContent_Declare(json
        URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(json)
endif()

add_executable(sph)
target_compile_features(sph PUBLIC cxx_std_14)

# Inject DIM as compile definition
target_compile_definitions(sph PRIVATE DIM=${SPH_DIM})

# Debug build with sanitizers (uncomment for debugging)
option(ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
if(ENABLE_SANITIZERS)
    target_compile_options(sph PUBLIC -fsanitize=address -fsanitize=undefined -g -O1 -fno-omit-frame-pointer)
    target_link_options(sph PUBLIC -fsanitize=address -fsanitize=undefined)
    message(STATUS "Building with AddressSanitizer and UndefinedBehaviorSanitizer enabled")
endif()

target_compile_options(sph
  PUBLIC
    -Wall
    -Wno-sign-compare
    -Wno-maybe-uninitialized
    $<$<NOT:$<BOOL:${ENABLE_SANITIZERS}>>:-funroll-loops>
    # Removed -ffast-math as it can cause crashes with NaN/Inf in shock simulations
    # $<$<NOT:$<BOOL:${ENABLE_SANITIZERS}>>:-ffast-math>
  )
target_include_directories(sph PUBLIC include ${HDF5_INCLUDE_DIRS})
target_link_libraries(sph
  PUBLIC
    OpenMP::OpenMP_CXX
    Boost::boost
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    nlohmann_json::nlohmann_json
  )

add_subdirectory(include)
add_subdirectory(src)
