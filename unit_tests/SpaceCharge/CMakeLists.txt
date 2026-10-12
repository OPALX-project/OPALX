file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestSpaceChargeConfig USE_GTEST_MAIN)
add_opalx_test(TestFFT2D5Algorithm USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestFFT2D5Numerics USE_GTEST_MAIN RUN_SERIAL)
# Separate processes ensure a CUDA launch failure in one implementation does not
# prevent the other implementations from running. Two TUs exercise header instantiations.
add_opalx_test(TestSelfFieldClear USE_GTEST_MAIN RUN_SERIAL TIMEOUT 30
  SOURCES TestSelfFieldClear.cpp SelfFieldClearTestOps.cpp
  ARGS --gtest_filter=SelfFieldClearTest.DeepCopy)
if(BUILD_TESTING AND TARGET TestSelfFieldClear)
  foreach(_case IN ITEMS OriginalAssignment FencedAssignment)
    add_test(NAME TestSelfFieldClear${_case}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:TestSelfFieldClear> --gtest_filter=SelfFieldClearTest.${_case}
        ${MPIEXEC_POSTFLAGS})
    set_tests_properties(TestSelfFieldClear${_case} PROPERTIES
      TIMEOUT 30 RUN_SERIAL TRUE LABELS "unit;spacecharge;cuda-diagnostic"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  endforeach()
  set_tests_properties(TestSelfFieldClear PROPERTIES
    LABELS "unit;spacecharge;cuda-diagnostic"
    ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
endif()
add_opalx_test(TestSpaceChargeSolver USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestParticleBinTraversal USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestCartesianPIC3DComponents USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestCartesianPIC3DAlgorithms USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestCartesianPIC3DFields USE_GTEST_MAIN RUN_SERIAL)
add_opalx_test(TestBeamBeamFieldServices USE_GTEST_MAIN RUN_SERIAL)
if(BUILD_TESTING AND TARGET TestBeamBeamFieldServices)
  set_tests_properties(TestBeamBeamFieldServices PROPERTIES
    TIMEOUT 30 LABELS "unit;beambeam;regression"
    ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  # Exercise empty-local witness ownership and the non-slab mirror case in CI,
  # not just in manually launched validation. Reuse the same gtest executable.
  if(MPIEXEC_EXECUTABLE)
    foreach(_ranks IN ITEMS 2 4)
      add_test(NAME TestBeamBeamFieldServicesMPI${_ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_ranks}
          ${MPIEXEC_PREFLAGS} $<TARGET_FILE:TestBeamBeamFieldServices>
          ${MPIEXEC_POSTFLAGS})
      set_tests_properties(TestBeamBeamFieldServicesMPI${_ranks} PROPERTIES
        TIMEOUT 30 PROCESSORS ${_ranks} RUN_SERIAL TRUE
        LABELS "unit;beambeam;regression;mpi"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
      if(_ranks EQUAL 4)
        set_tests_properties(TestBeamBeamFieldServicesMPI${_ranks} PROPERTIES DISABLED TRUE)
      endif()
    endforeach()
  endif()
endif()
add_opalx_test(TestPoissonSolver USE_GTEST_MAIN RUN_SERIAL)
