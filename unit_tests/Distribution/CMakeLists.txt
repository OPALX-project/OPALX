file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestGaussian USE_GTEST_MAIN)
add_opalx_test(TestMultiVariateGaussian USE_GTEST_MAIN)
add_opalx_test(TestFlatTop USE_GTEST_MAIN)
add_opalx_test(TestOpalFlatTop USE_GTEST_MAIN)
add_opalx_test(TestFromFile USE_GTEST_MAIN)
add_opalx_test(TestEmittedFromFile USE_GTEST_MAIN)
if(BUILD_TESTING AND TARGET TestEmittedFromFile)
  set_tests_properties(TestEmittedFromFile PROPERTIES
    TIMEOUT 30 RUN_SERIAL TRUE LABELS "unit;beambeam;regression"
    ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  if(MPIEXEC_EXECUTABLE)
    add_test(NAME TestBeamBeamTimedBirthsMPI4
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
        ${MPIEXEC_PREFLAGS} $<TARGET_FILE:TestEmittedFromFile> ${MPIEXEC_POSTFLAGS}
        --gtest_filter=EmittedFromFileTest.TwelveTimedBirthsPreserveIdentityAndBoundarySteps)
    set_tests_properties(TestBeamBeamTimedBirthsMPI4 PROPERTIES
      DISABLED TRUE
      TIMEOUT 30 PROCESSORS 4 RUN_SERIAL TRUE
      LABELS "unit;beambeam;regression;mpi"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  endif()
endif()
add_opalx_test(TestUniform USE_GTEST_MAIN)
