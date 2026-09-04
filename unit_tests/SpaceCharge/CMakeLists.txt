file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestSpaceChargeConfig USE_GTEST_MAIN)
add_opalx_test(TestFFT2D5Algorithm USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestSpaceChargeSolver USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestParticleBinTraversal USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestCartesianPICComponents USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestCartesianPICFields USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestPoissonSolver USE_GTEST_MAIN RUN_SERIAL)

set(_fft2d5_regression_dir "${CMAKE_CURRENT_BINARY_DIR}/FFT2D5EndToEnd")
file(MAKE_DIRECTORY "${_fft2d5_regression_dir}")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/data/FFT2D5EndToEnd.in"
  "${_fft2d5_regression_dir}/FFT2D5EndToEnd.in" COPYONLY)

if(BUILD_TESTING)
  add_test(
    NAME FFT2D5EndToEnd
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:opalx_exe> FFT2D5EndToEnd.in --info 2)
  set_tests_properties(
    FFT2D5EndToEnd
    PROPERTIES WORKING_DIRECTORY "${_fft2d5_regression_dir}"
               LABELS "unit;regression;space-charge;fft2d5"
               RUN_SERIAL TRUE
               TIMEOUT 120
               PASS_REGULAR_EXPRESSION "Total FieldSolver calls: 16")
endif()
