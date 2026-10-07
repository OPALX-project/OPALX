file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestSpaceChargeConfig USE_GTEST_MAIN)
add_opalx_test(TestFFT2D5Algorithm USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestFFT2D5Numerics USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestSpaceChargeSolver USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestParticleBinTraversal USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestCartesianPIC3DComponents USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestCartesianPIC3DAlgorithms USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestCartesianPIC3DFields USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestPoissonSolver USE_GTEST_MAIN RUN_SERIAL)

if(OPALX_ENABLE_BH)
  add_opalx_test(TestBarnesHutAlgorithm USE_GTEST_MAIN RUN_SERIAL TIMEOUT 900)
  # Same binary on 4 ranks: Barnes-Hut migrates particles, so the multi-rank paths
  # (copy-back growth/shrink, halos, redistribution) only run here.
  if(DEFINED MPIEXEC_EXECUTABLE AND BUILD_TESTING)
    add_test(NAME TestBarnesHutAlgorithm_4ranks
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                     $<TARGET_FILE:TestBarnesHutAlgorithm>)
    set_tests_properties(TestBarnesHutAlgorithm_4ranks PROPERTIES
                         TIMEOUT 900 RUN_SERIAL TRUE LABELS "unit;nbody")
  endif()
endif()
