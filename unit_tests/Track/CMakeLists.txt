file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestDecayChain USE_GTEST_MAIN)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(BUILD_TESTING AND Python3_Interpreter_FOUND)
  set(_cof_mpi_flags)
  foreach(_flag IN LISTS MPIEXEC_PREFLAGS)
    list(APPEND _cof_mpi_flags "--mpi-preflag=${_flag}")
  endforeach()
  foreach(_flag IN LISTS MPIEXEC_POSTFLAGS)
    list(APPEND _cof_mpi_flags "--mpi-postflag=${_flag}")
  endforeach()
  foreach(_ranks IN ITEMS 1 2)
    add_test(NAME TestCofHousekeeping${_ranks}Rank
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_cof_housekeeping.py
        --executable $<TARGET_FILE:opalx_exe> --mpiexec ${MPIEXEC_EXECUTABLE}
        --numproc-flag=${MPIEXEC_NUMPROC_FLAG} --ranks ${_ranks} ${_cof_mpi_flags})
    set_tests_properties(TestCofHousekeeping${_ranks}Rank PROPERTIES
      TIMEOUT 300 PROCESSORS ${_ranks} LABELS "integration;cof")
  endforeach()
endif()
