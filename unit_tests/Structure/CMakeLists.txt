file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestBinningCmd USE_GTEST_MAIN)
add_opalx_test(TestBinConfigWriter USE_GTEST_MAIN)
add_opalx_test(TestDirichletPlaneWriter USE_GTEST_MAIN)
add_opalx_test(TestLaser USE_GTEST_MAIN)
add_opalx_test(TestBeam USE_GTEST_MAIN)
add_opalx_test(TestBeamBeamDiagnosticsWriter USE_GTEST_MAIN RUN_SERIAL TIMEOUT 30)
if(BUILD_TESTING AND TARGET TestBeamBeamDiagnosticsWriter)
  set_tests_properties(TestBeamBeamDiagnosticsWriter PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
endif()
if(BUILD_TESTING AND TARGET TestBeamBeamDiagnosticsWriter AND MPIEXEC_EXECUTABLE)
  add_test(NAME TestBeamBeamDiagnosticsWriterMPI2
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:TestBeamBeamDiagnosticsWriter> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(TestBeamBeamDiagnosticsWriterMPI2 PROPERTIES
    TIMEOUT 30 PROCESSORS 2 RUN_SERIAL TRUE LABELS "unit;beambeam;mpi"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
endif()
add_opalx_test(TestLossDataSink USE_GTEST_MAIN)
add_opalx_test(TestCheckpointFile USE_GTEST_MAIN)
add_opalx_test(TestMeshGenerator USE_GTEST_MAIN)
