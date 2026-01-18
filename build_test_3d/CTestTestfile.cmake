# CMake generated Testfile for 
# Source directory: /Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation
# Build directory: /Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_test_3d
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_test_3d/test_kernel_gravity[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_test_3d/test_neighbor_search[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_test_3d/test_morton[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_test_3d/test_hdf5_resume[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_test_3d/test_morton_tree_integration[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_test_3d/test_shock_tube_3d_regression[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_test_3d/test_softening_lookup_table[1]_include.cmake")
subdirs("_deps/json-build")
subdirs("include")
subdirs("src")
subdirs("_deps/googletest-build")
