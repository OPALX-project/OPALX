cmake_minimum_required(VERSION 3.20)

foreach(required_variable IN ITEMS OPALX_EXE INPUT_FILE TEST_WORKDIR)
  if(NOT DEFINED ${required_variable} OR "${${required_variable}}" STREQUAL "")
    message(FATAL_ERROR "${required_variable} must be provided")
  endif()
endforeach()

if(NOT EXISTS "${OPALX_EXE}")
  message(FATAL_ERROR "OPALX executable does not exist: ${OPALX_EXE}")
endif()
if(NOT EXISTS "${INPUT_FILE}")
  message(FATAL_ERROR "Regression input does not exist: ${INPUT_FILE}")
endif()
if("${TEST_WORKDIR}" STREQUAL "/" OR "${TEST_WORKDIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "Refusing unsafe test working directory: ${TEST_WORKDIR}")
endif()

file(REMOVE_RECURSE "${TEST_WORKDIR}")
file(MAKE_DIRECTORY "${TEST_WORKDIR}")
file(COPY "${INPUT_FILE}" DESTINATION "${TEST_WORKDIR}")
get_filename_component(input_name "${INPUT_FILE}" NAME)

function(run_opalx phase)
  execute_process(
    COMMAND "${OPALX_EXE}" ${ARGN} "${input_name}"
    WORKING_DIRECTORY "${TEST_WORKDIR}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr
  )
  set(output "${stdout}\n${stderr}")
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${phase} failed with exit code ${result}:\n${output}")
  endif()
  set(${phase}_OUTPUT "${output}" PARENT_SCOPE)
endfunction()

run_opalx(FRESH --info 2)

set(checkpoint_file "${TEST_WORKDIR}/CheckpointMonitorRestart_checkpoint.h5")
set(m1_file "${TEST_WORKDIR}/M1.h5")
set(m2_file "${TEST_WORKDIR}/M2.h5")

if(NOT EXISTS "${checkpoint_file}")
  message(FATAL_ERROR "Fresh run did not create the checkpoint")
endif()
if(NOT EXISTS "${m1_file}")
  message(FATAL_ERROR "Fresh run did not reach M1")
endif()
if(EXISTS "${m2_file}")
  message(FATAL_ERROR "Fresh run unexpectedly reached M2")
endif()

run_opalx(RESTART --info 2 --restart "${checkpoint_file}")

string(FIND "${RESTART_OUTPUT}" "Starting track with dt" tracking_position)
string(FIND "${RESTART_OUTPUT}" "Save 'M1.h5'" first_m1_save_position)
string(FIND "${RESTART_OUTPUT}" "Save 'M2.h5'" first_m2_save_position)

if(tracking_position EQUAL -1)
  message(FATAL_ERROR "Restart output does not contain the tracking-loop marker")
endif()
if(first_m1_save_position EQUAL -1)
  message(FATAL_ERROR "Restart did not physically reach M1")
endif()
if(first_m1_save_position LESS tracking_position)
  message(FATAL_ERROR "M1 wrote output during restart preparation")
endif()
if(NOT first_m2_save_position EQUAL -1 OR EXISTS "${m2_file}")
  message(FATAL_ERROR "Unreachable M2 wrote output during restart preparation")
endif()
