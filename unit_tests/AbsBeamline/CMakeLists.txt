file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestConstantEFieldCavity USE_GTEST_MAIN)
add_opalx_test(TestMultipoleTBase USE_GTEST_MAIN)
add_opalx_test(TestMultipoleTStraight USE_GTEST_MAIN)
add_opalx_test(TestMultipoleTCurvedConstRadius USE_GTEST_MAIN)
