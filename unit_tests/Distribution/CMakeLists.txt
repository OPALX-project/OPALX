file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestGaussian USE_GTEST_MAIN)
add_opalx_test(TestMultiVariateGaussian USE_GTEST_MAIN)
add_opalx_test(TestFlatTop USE_GTEST_MAIN)
add_opalx_test(TestOpalFlatTop USE_GTEST_MAIN)
add_opalx_test(TestFromFile USE_GTEST_MAIN)
add_opalx_test(TestEmittedFromFile USE_GTEST_MAIN)
