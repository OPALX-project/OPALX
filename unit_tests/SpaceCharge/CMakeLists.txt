file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestSpaceChargeConfig USE_GTEST_MAIN)
# Temporarily disabled: CUDA 13 on GH200 reports cudaErrorInvalidResourceHandle
# during FFT2D5 initialization. Re-enable after the failure is resolved.
# add_opalx_test(TestFFT2D5Algorithm USE_GTEST_MAIN RUN_SERIAL)
# add_opalx_test(TestFFT2D5Numerics USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestSpaceChargeSolver USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestParticleBinTraversal USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestCartesianPIC3DComponents USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestCartesianPIC3DAlgorithms USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestCartesianPIC3DFields USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestPoissonSolver USE_GTEST_MAIN RUN_SERIAL)
