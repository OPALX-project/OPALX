file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestPolynomialTimeDependence USE_GTEST_MAIN)
add_opalx_test(TestSinusoidalTimeDependence USE_GTEST_MAIN)
add_opalx_test(TestSplineTimeDependence USE_GTEST_MAIN)
add_opalx_test(TestQuaternion USE_GTEST_MAIN)
