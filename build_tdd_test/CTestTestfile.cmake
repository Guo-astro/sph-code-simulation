# CMake generated Testfile for 
# Source directory: /Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation
# Build directory: /Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_tdd_test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_tdd_test/test_kernel_gravity[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_tdd_test/test_neighbor_search[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_tdd_test/test_morton[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_tdd_test/test_hdf5_resume[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_tdd_test/test_morton_tree_integration[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_tdd_test/test_shock_tube_3d_regression[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_tdd_test/test_softening_lookup_table[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_tdd_test/test_hernquist_katz_lookup_table[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_tdd_test/test_analytic_gravity_3d[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_tdd_test/test_gravity_soa[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_tdd_test/test_gravity_lookup_integration[1]_include.cmake")
include("/Users/yansongguo/personal/research/shamrock-wrapper/sph-code-simulation/build_tdd_test/test_gravity_first_principles[1]_include.cmake")
subdirs("_deps/json-build")
subdirs("include")
subdirs("src")
subdirs("_deps/googletest-build")
