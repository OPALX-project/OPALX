file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

# Tests for Boost/GSL replacement utilities
add_opalx_test(TestBiMap USE_GTEST_MAIN)
add_opalx_test(TestCubicSpline USE_GTEST_MAIN)
add_opalx_test(TestLinearSpline USE_GTEST_MAIN)
add_opalx_test(TestHistogram USE_GTEST_MAIN)
add_opalx_test(TestOpalException USE_GTEST_MAIN)

# Attribute parsing coverage (OPAL parser interface)
add_opalx_test(TestPredefinedString USE_GTEST_MAIN)
