add_opalx_test(TestMacro LABELS parser)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(BUILD_TESTING AND Python3_Interpreter_FOUND)
  # Launch each OPALX invocation separately, so an expected MPI_Abort on bad
  # input cannot terminate the Python test runner as well.
  set(_parse_only_launcher)
  if(MPIEXEC_EXECUTABLE)
    set(_parse_only_launcher ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                             ${MPIEXEC_PREFLAGS})
  endif()

  # Keep the script with the build artifacts when CI moves them between jobs.
  configure_file(test_parse_only.py test_parse_only.py COPYONLY)
  configure_file(run_parse_only.cmake.in run_parse_only.cmake @ONLY)
  # Resolve CMake on the test host: CSCS's build-time uenv may be absent there.
  add_test(NAME ParseOnlyCLI
    COMMAND cmake "-DOPALX_EXECUTABLE=$<TARGET_FILE:opalx_exe>"
            -P "${CMAKE_CURRENT_BINARY_DIR}/run_parse_only.cmake")
  set_tests_properties(ParseOnlyCLI PROPERTIES TIMEOUT 180 LABELS "integration;parser")
endif()
