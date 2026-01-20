file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

# Tests for Boost/GSL replacement utilities
add_opalx_test(TestBiMap USE_GTEST_MAIN NO_MPI)
add_opalx_test(TestCubicSpline USE_GTEST_MAIN NO_MPI)
add_opalx_test(TestFFT USE_GTEST_MAIN NO_MPI)
add_opalx_test(TestRandom USE_GTEST_MAIN NO_MPI)
add_opalx_test(TestHistogram USE_GTEST_MAIN NO_MPI)
add_opalx_test(TestQuasiRandom USE_GTEST_MAIN NO_MPI)
add_opalx_test(TestGSLCompat USE_GTEST_MAIN NO_MPI)
add_opalx_test(TestGSLComplex USE_GTEST_MAIN NO_MPI)
add_opalx_test(TestGSLMatrix USE_GTEST_MAIN NO_MPI)
add_opalx_test(TestGSLBLAS USE_GTEST_MAIN NO_MPI)
add_opalx_test(TestGSLLinalg USE_GTEST_MAIN NO_MPI)
add_opalx_test(TestGSLEigen USE_GTEST_MAIN NO_MPI)
add_opalx_test(TestGSLIntegration USE_GTEST_MAIN NO_MPI)

