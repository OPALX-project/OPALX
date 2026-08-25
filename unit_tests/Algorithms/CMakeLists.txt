file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestPolynomialTimeDependence USE_GTEST_MAIN)
add_opalx_test(TestSinusoidalTimeDependence USE_GTEST_MAIN)
add_opalx_test(TestSplineTimeDependence USE_GTEST_MAIN)
add_opalx_test(TestQuaternion USE_GTEST_MAIN)
add_opalx_test(TestCoordinateSystemTrafo USE_GTEST_MAIN)
add_opalx_test(TestIndexMap USE_GTEST_MAIN)
add_opalx_test(TestOrbitThreader USE_GTEST_MAIN)
add_opalx_test(TestDistributionMoments USE_GTEST_MAIN)
add_opalx_test(TestStepSizeConfig USE_GTEST_MAIN)

if(BUILD_TESTING)
  add_test(
    NAME TestCheckpointMonitorRestart
    COMMAND
      ${MPIEXEC_EXECUTABLE}
      ${MPIEXEC_NUMPROC_FLAG}
      1
      ${MPIEXEC_PREFLAGS}
      ${CMAKE_COMMAND}
      "-DOPALX_EXE=$<TARGET_FILE:opalx_exe>"
      "-DINPUT_FILE=${CMAKE_CURRENT_SOURCE_DIR}/CheckpointMonitorRestart.in"
      "-DTEST_WORKDIR=${CMAKE_CURRENT_BINARY_DIR}/TestCheckpointMonitorRestart.run"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/TestCheckpointMonitorRestart.cmake"
  )
  set_tests_properties(
    TestCheckpointMonitorRestart
    PROPERTIES
      LABELS "integration;checkpoint;monitor"
      RUN_SERIAL TRUE
      TIMEOUT 60
  )
endif()
