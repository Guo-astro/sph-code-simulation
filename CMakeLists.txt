cmake_minimum_required(VERSION 3.13)
project(sphcode CXX)

# =============================================================================
# Build type: default to Release for optimized performance
# =============================================================================
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Debug, Release, RelWithDebInfo, MinSizeRel)" FORCE)
    message(STATUS "No build type specified, defaulting to Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Dimension parameter (1D, 2D, or 3D)
# Can be set via: cmake -DSPH_DIM=1 ..
set(SPH_DIM "2" CACHE STRING "Spatial dimension (1, 2, or 3)")
set_property(CACHE SPH_DIM PROPERTY STRINGS "1" "2" "3")
message(STATUS "Building for ${SPH_DIM}D simulations (DIM=${SPH_DIM})")

# =============================================================================
# Performance options
# =============================================================================
option(ENABLE_NATIVE_ARCH "Enable -march=native for best performance on this CPU (disable for portable builds)" ON)
option(ENABLE_LTO "Enable Link-Time Optimization (slower build, faster runtime)" OFF)

# OpenMP setup - handle macOS AppleClang specially
if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # AppleClang needs Homebrew's libomp
    execute_process(
        COMMAND brew --prefix libomp
        OUTPUT_VARIABLE LIBOMP_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(LIBOMP_PREFIX)
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
        message(STATUS "Using Homebrew OpenMP at ${LIBOMP_PREFIX}")
    endif()
endif()
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

# =============================================================================
# Compiler flags for performance
# =============================================================================
# Release-mode optimizations (only when not using sanitizers)
if(NOT ENABLE_SANITIZERS)
    # Use -O3 for aggressive optimization in Release builds
    target_compile_options(sph PUBLIC
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:RelWithDebInfo>:-O2>
    )
    
    # Native architecture tuning (best performance on build machine)
    if(ENABLE_NATIVE_ARCH AND NOT CMAKE_CROSSCOMPILING)
        include(CheckCXXCompilerFlag)
        check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
        if(COMPILER_SUPPORTS_MARCH_NATIVE)
            target_compile_options(sph PUBLIC -march=native)
            message(STATUS "Enabled -march=native for CPU-specific optimizations")
        endif()
    endif()
    
    # Link-Time Optimization (optional, slower build but faster runtime)
    if(ENABLE_LTO)
        include(CheckIPOSupported)
        check_ipo_supported(RESULT LTO_SUPPORTED OUTPUT LTO_ERROR)
        if(LTO_SUPPORTED)
            set_property(TARGET sph PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
            message(STATUS "Enabled Link-Time Optimization (LTO)")
        else()
            message(WARNING "LTO requested but not supported: ${LTO_ERROR}")
        endif()
    endif()
endif()

# =============================================================================
# SIMD detection for tree optimizations
# =============================================================================
include(CheckCXXCompilerFlag)

# Check for AVX2 support
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
if(COMPILER_SUPPORTS_AVX2)
    target_compile_definitions(sph PRIVATE SPH_HAS_AVX2)
    target_compile_options(sph PRIVATE -mavx2 -mfma)
    message(STATUS "Enabled AVX2 SIMD optimizations")
endif()

# Check for AVX-512 support
check_cxx_compiler_flag("-mavx512f" COMPILER_SUPPORTS_AVX512)
if(COMPILER_SUPPORTS_AVX512)
    target_compile_definitions(sph PRIVATE SPH_HAS_AVX512)
    target_compile_options(sph PRIVATE -mavx512f)
    message(STATUS "Enabled AVX-512 SIMD optimizations")
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

# =============================================================================
# Google Test Setup
# =============================================================================
option(BUILD_TESTS "Build unit tests" ON)
if(BUILD_TESTS)
    message(STATUS "Building tests enabled")
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
    
    enable_testing()
    include(GoogleTest)
    
    # Kernel gravity tests (unit tests for 1D/2D/3D kernel gravity functions)
    # We need to link the gravity_force source since the test uses static member functions
    add_executable(test_kernel_gravity 
        tests/test_kernel_gravity.cpp
        src/gravity_force.cpp
    )
    target_compile_features(test_kernel_gravity PUBLIC cxx_std_14)
    target_compile_definitions(test_kernel_gravity PRIVATE DIM=${SPH_DIM})
    target_include_directories(test_kernel_gravity PUBLIC include ${HDF5_INCLUDE_DIRS})
    target_link_libraries(test_kernel_gravity
        GTest::gtest_main
        GTest::gmock
        OpenMP::OpenMP_CXX
        Boost::boost
    )
    target_compile_options(test_kernel_gravity PRIVATE -Wall -Wno-sign-compare)
    gtest_discover_tests(test_kernel_gravity)

    # Neighbor search validation tests (tree vs exhaustive search comparison)
    add_executable(test_neighbor_search
        tests/test_neighbor_search.cpp
        src/bhtree.cpp
        src/exhaustive_search.cpp
    )
    target_compile_features(test_neighbor_search PUBLIC cxx_std_14)
    target_compile_definitions(test_neighbor_search PRIVATE DIM=${SPH_DIM})
    target_include_directories(test_neighbor_search PUBLIC include ${HDF5_INCLUDE_DIRS})
    target_link_libraries(test_neighbor_search
        GTest::gtest_main
        GTest::gmock
        OpenMP::OpenMP_CXX
        Boost::boost
    )
    target_compile_options(test_neighbor_search PRIVATE -Wall -Wno-sign-compare)
    gtest_discover_tests(test_neighbor_search)

    # Shock tube ghost boundary tests (disabled - test file not yet created)
    # add_executable(test_shock_tube_ghost
    #     tests/test_shock_tube_ghost.cpp
    #     src/shock_tube_ghost.cpp
    # )
    # target_compile_features(test_shock_tube_ghost PUBLIC cxx_std_14)
    # target_compile_definitions(test_shock_tube_ghost PRIVATE DIM=${SPH_DIM})
    # target_include_directories(test_shock_tube_ghost PUBLIC include ${HDF5_INCLUDE_DIRS})
    # target_link_libraries(test_shock_tube_ghost
    #     GTest::gtest_main
    #     GTest::gmock
    #     OpenMP::OpenMP_CXX
    #     Boost::boost
    # )
    # target_compile_options(test_shock_tube_ghost PRIVATE -Wall -Wno-sign-compare)
    # gtest_discover_tests(test_shock_tube_ghost)

    # Boundary particle evolution tests (disabled - test file not yet created)
    # add_executable(test_boundary_particle_evolution
    #     tests/test_boundary_particle_evolution.cpp
    #     src/shock_tube_ghost.cpp
    # )
    # target_compile_features(test_boundary_particle_evolution PUBLIC cxx_std_14)
    # target_compile_definitions(test_boundary_particle_evolution PRIVATE DIM=${SPH_DIM})
    # target_include_directories(test_boundary_particle_evolution PUBLIC include ${HDF5_INCLUDE_DIRS})
    # target_link_libraries(test_boundary_particle_evolution
    #     GTest::gtest_main
    #     GTest::gmock
    #     OpenMP::OpenMP_CXX
    #     Boost::boost
    # )
    # target_compile_options(test_boundary_particle_evolution PRIVATE -Wall -Wno-sign-compare)
    # gtest_discover_tests(test_boundary_particle_evolution)

    # Morton code unit tests for tree optimization
    add_executable(test_morton
        tests/test_morton.cpp
        src/morton.cpp
    )
    target_compile_features(test_morton PUBLIC cxx_std_14)
    target_compile_definitions(test_morton PRIVATE DIM=${SPH_DIM})
    target_include_directories(test_morton PUBLIC include ${HDF5_INCLUDE_DIRS})
    target_link_libraries(test_morton
        GTest::gtest_main
        GTest::gmock
        OpenMP::OpenMP_CXX
        Boost::boost
    )
    target_compile_options(test_morton PRIVATE -Wall -Wno-sign-compare)
    gtest_discover_tests(test_morton)

    # 3D Shock tube regression tests
    # Tests numerical reproducibility against stored reference data
    # L2 error threshold: 1e-10 to catch any numerical changes
    add_executable(test_shock_tube_3d_regression
        tests/regression/test_shock_tube_3d_regression.cpp
    )
    target_compile_features(test_shock_tube_3d_regression PUBLIC cxx_std_17)
    target_compile_definitions(test_shock_tube_3d_regression PRIVATE DIM=${SPH_DIM})
    target_include_directories(test_shock_tube_3d_regression PUBLIC include ${HDF5_INCLUDE_DIRS})
    target_link_libraries(test_shock_tube_3d_regression
        GTest::gtest_main
        GTest::gmock
        nlohmann_json::nlohmann_json
        OpenMP::OpenMP_CXX
        Boost::boost
    )
    target_compile_options(test_shock_tube_3d_regression PRIVATE -Wall -Wno-sign-compare)
    gtest_discover_tests(test_shock_tube_3d_regression)
endif()
