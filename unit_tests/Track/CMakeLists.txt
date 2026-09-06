# Parser/runtime tests need only Python's standard library, not analysis packages.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
  set(_cof_mpi_args)
  if(MPIEXEC_EXECUTABLE)
    list(APPEND _cof_mpi_args --mpiexec ${MPIEXEC_EXECUTABLE})
  endif()
  add_test(NAME TestCofCommand
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_cof_command.py
            --opalx $<TARGET_FILE:opalx_exe> ${_cof_mpi_args})
  set_tests_properties(TestCofCommand PROPERTIES TIMEOUT 240 RUN_SERIAL TRUE LABELS "integration;cof")
else()
  message(STATUS "Python3 interpreter unavailable: COF command integration test not registered")
endif()
