file(RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Adding unit tests found in ${_relPath}")

add_opalx_test(TestDecayChain USE_GTEST_MAIN)

if(BUILD_TESTING)
  set(_cof_config "${CMAKE_CURRENT_BINARY_DIR}/CofHousekeeping-$<CONFIG>.cmake")
  file(GENERATE OUTPUT "${_cof_config}" CONTENT
"set(COF_EXECUTABLE [==[$<TARGET_FILE:opalx_exe>]==])
set(COF_MPIEXEC [==[${MPIEXEC_EXECUTABLE}]==])
set(COF_NUMPROC_FLAG [==[${MPIEXEC_NUMPROC_FLAG}]==])
set(COF_MPI_PREFLAGS [==[${MPIEXEC_PREFLAGS}]==])
set(COF_MPI_POSTFLAGS [==[${MPIEXEC_POSTFLAGS}]==])
set(COF_WORK_ROOT [==[${CMAKE_CURRENT_BINARY_DIR}]==])
")
  foreach(_ranks IN ITEMS 1 2)
    add_test(NAME TestCofHousekeeping${_ranks}Rank
      COMMAND ${CMAKE_COMMAND} "-DCOF_CONFIG=${_cof_config}" "-DCOF_RANKS=${_ranks}"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/CheckCofHousekeeping.cmake)
    set_tests_properties(TestCofHousekeeping${_ranks}Rank PROPERTIES
      TIMEOUT 300 PROCESSORS ${_ranks} LABELS "integration;cof")
  endforeach()
endif()
