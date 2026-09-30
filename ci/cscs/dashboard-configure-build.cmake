cmake_minimum_required(VERSION 3.20)

set(DASHBOARD_PROJECT "OPALX")

if(NOT DEFINED BUILD_TYPE)
  set(BUILD_TYPE Debug)
endif()

if(NOT DEFINED CDASH_LABEL)
  set(CDASH_LABEL "branch")
endif()

if(NOT DEFINED BUILD_DIR)
  message(FATAL_ERROR "BUILD_DIR must be defined")
endif()

# The full build name here gets overwritten in test phase, so leave this blank for now
set(TEST_INFO "")

# --- CDash metadata ---
set(CTEST_SITE "${CTEST_SITE}")
set(CTEST_BUILD_CONFIGURATION ${BUILD_TYPE})
set(CTEST_BUILD_NAME "${CDASH_LABEL}-${BUILD_ARCH}-${BUILD_TYPE}-${TEST_INFO}")

set(CTEST_SOURCE_DIRECTORY "$ENV{CI_PROJECT_DIR}")
set(CTEST_BINARY_DIRECTORY "${BUILD_DIR}")
set(CTEST_CMAKE_GENERATOR "Ninja")
set(CTEST_GROUP "Experimental")

# --- start a new build in CDash ---
ctest_start(Experimental GROUP "${CTEST_GROUP}")

# --- Initialize base configure command as a CMake LIST ---
set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND}")
string(APPEND CTEST_CONFIGURE_COMMAND " -S${CTEST_SOURCE_DIRECTORY}")
string(APPEND CTEST_CONFIGURE_COMMAND " -B${CTEST_BINARY_DIRECTORY}")
string(APPEND CTEST_CONFIGURE_COMMAND " -G${CTEST_CMAKE_GENERATOR}")
string(APPEND CTEST_CONFIGURE_COMMAND " --preset=${PRESET}")
string(APPEND CTEST_CONFIGURE_COMMAND " -DCMAKE_BUILD_TYPE=${BUILD_TYPE}")
string(APPEND CTEST_CONFIGURE_COMMAND " -DCMAKE_BUILD_RPATH_USE_ORIGIN=ON")
string(APPEND CTEST_CONFIGURE_COMMAND " -DOPALX_USE_STANDARD_FOLDERS=ON")

# --- Append remaining static flags ---
string(APPEND CTEST_CONFIGURE_COMMAND " -DIPPL_ENABLE_SOLVERS=ON")
string(APPEND CTEST_CONFIGURE_COMMAND " -DIPPL_MARK_FAILING_TESTS=ON")
if(DEFINED Kokkos_ARCH_FLAG)
  string(APPEND CTEST_CONFIGURE_COMMAND " -D${Kokkos_ARCH_FLAG}=ON")
endif()

# ---------------------------------
# cmake-format: off
# ---------------------------------
# --- Forward variables cleanly ---
set(VARS_TO_FORWARD
  # Opal options  
  OPALX_PLATFORMS
  OPALX_OPENMP_THREADS
  OPALX_ENABLE_SCRIPTS
  # Ippl 
  IPPL_GIT_TAG
  Heffte_VERSION
  Kokkos_VERSION
  Kokkos_ENABLE_DEBUG_BOUNDS_CHECK
  # srun/mpi
  MPIEXEC_EXECUTABLE
  MPIEXEC_PREFLAGS      
  MPIEXEC_MAX_NUMPROCS
)
# ---------------------------------
# cmake-format: on
# ---------------------------------

foreach(VAR IN LISTS VARS_TO_FORWARD)
  if(DEFINED ${VAR})
    set(VAL "${${VAR}}")
    if("${VAL}" MATCHES ";")
      # Force CMake to treat VAR it as a literal string cache entry to preserve semicolons
      string(APPEND CTEST_CONFIGURE_COMMAND " -D${VAR}:STRING=${VAL}")
      # Expose VAR to the local scope so it is passed through
      set(${VAR} "${VAL}")
    else()
      # Standard scalar variable (no semicolons)
      string(APPEND CTEST_CONFIGURE_COMMAND " -D${VAR}=${VAL}")
    endif()
  endif()
endforeach()

# --- Output our configure command for debug purposes---
message("Final CTest configure command: ${CTEST_CONFIGURE_COMMAND}")

# --- configure & build ---
ctest_configure(RETURN_VALUE configure_result)
# --- submit configure results immediately, so they reach CDash even if the
# build later hangs and the job is killed by the SLURM timelimit ---
ctest_submit(PARTS Configure)

# Optional GH200 include probe. Run before compilation so its report survives a
# build failure; diagnostic failures must not replace the normal build result.
if(OPALX_DIAG_DESUL_HEADERS AND configure_result EQUAL 0)
  find_program(DIAGNOSTIC_PYTHON NAMES python3)
  if(DIAGNOSTIC_PYTHON)
    execute_process(
      COMMAND "${DIAGNOSTIC_PYTHON}"
              "${CTEST_SOURCE_DIRECTORY}/ci/cscs/diagnose_desul_headers.py"
              "${CTEST_BINARY_DIRECTORY}"
      RESULT_VARIABLE diagnostic_result
      OUTPUT_VARIABLE diagnostic_output
      ERROR_VARIABLE diagnostic_error
      TIMEOUT 150
    )
    set(diagnostic_report "${diagnostic_output}\n${diagnostic_error}\nDiagnostic exit status: ${diagnostic_result}\n")
  else()
    set(diagnostic_report "Desul header diagnostic unavailable: python3 was not found.\n")
  endif()
  set(diagnostic_file "${CTEST_BINARY_DIRECTORY}/desul-header-diagnostic.log")
  file(WRITE "${diagnostic_file}" "${diagnostic_report}")
  message("${diagnostic_report}")
  list(APPEND CTEST_NOTES_FILES "${diagnostic_file}")
  ctest_submit(
    PARTS Notes
    RETURN_VALUE diagnostic_submit_result
    CAPTURE_CMAKE_ERROR diagnostic_submit_error
  )
endif()

ctest_build(RETURN_VALUE build_result)

# --- fail if any test failed ---
# configure results were already submitted above; on failure also submit the
# build results so errors are visible on CDash immediately (the test phase
# submits everything else when it succeeds)
if(configure_result OR build_result)
  # submit configure + build results
  ctest_submit()
  # make sure to fail the build if configure or build failed
  message(FATAL_ERROR "CTest reported configure/build failures")
endif()

string(ASCII 27 ESC)
set(BLUE "${ESC}[34m")
set(RESET "${ESC}[0m")
message("${BLUE}# ---------------------------------${RESET}")
message("${BLUE}To view ALL configure/build/test results and error logs visit: ${RESET}")
message("${BLUE}https://my.cdash.org/index.php?project=${DASHBOARD_PROJECT}${RESET}")
message("${BLUE}For this PR visit: ${RESET}")
message(
  "${BLUE}https://my.cdash.org/index.php?project=${DASHBOARD_PROJECT}&filtercount=1&showfilters=1&field1=buildname&compare1=63&value1=${CDASH_LABEL}${RESET}"
)
message("${BLUE}# ---------------------------------${RESET}")
