file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestFM2DMagnetoStatic USE_GTEST_MAIN)
add_opalx_test(TestConstEzField USE_GTEST_MAIN)
add_opalx_test(TestAstra1DDynamic USE_GTEST_MAIN)
add_opalx_test(TestAstra1DMagnetoStatic USE_GTEST_MAIN)
add_opalx_test(TestFM2DDynamic USE_GTEST_MAIN)