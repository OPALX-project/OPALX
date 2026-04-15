file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestOpalMultipoleT USE_GTEST_MAIN)
add_opalx_test(TestOpalPolynomialTimeDependence USE_GTEST_MAIN)
add_opalx_test(TestOpalSplineTimeDependence USE_GTEST_MAIN)
add_opalx_test(TestOpalSinusoidalTimeDependence USE_GTEST_MAIN)
add_opalx_test(TestOpalVariableRFCavity USE_GTEST_MAIN)
