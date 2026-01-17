# CMake generated Testfile for 
# Source directory: /Users/kunishigehana/sph-code-simulation
# Build directory: /Users/kunishigehana/sph-code-simulation/build_test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
include("/Users/kunishigehana/sph-code-simulation/build_test/test_kernel_gravity[1]_include.cmake")
include("/Users/kunishigehana/sph-code-simulation/build_test/test_neighbor_search[1]_include.cmake")
include("/Users/kunishigehana/sph-code-simulation/build_test/test_morton[1]_include.cmake")
subdirs("_deps/json-build")
subdirs("include")
subdirs("src")
subdirs("_deps/googletest-build")
