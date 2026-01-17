
# -----------------------------------------------------------------------------
# OPALXOptions.cmake
# ~~~
# High-level OPALX configuration options (use access)
# Handles platform/backend selection for OPALX: SERIAL, OPENMP, CUDA, or CUDA + OPENMP.
#
# Responsibilities:
#   - Normalize and validate the value
#   - Enable relevant Kokkos, IPPL options
#   - Selecting default CMake build type
# ~~~
# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------
# Build type (user-facing)
# -----------------------------------------------------------------------------
set(BUILD_TYPE "Release" CACHE STRING "Build type: Debug, Release, RelWithDebInfo, MinSizeRel")
set(CMAKE_BUILD_TYPE "${BUILD_TYPE}" CACHE STRING "" FORCE)

# -----------------------------------------------------------------------------
# Unit Test (user-facing)
# -----------------------------------------------------------------------------

option(OPALX_ENABLE_UNIT_TESTS "Enable unit tests using GoogleTest" OFF)
message(STATUS "🔧 OPALX unit tests: ${OPALX_ENABLE_UNIT_TESTS}")

# -----------------------------------------------------------------------------
# Platform selection (user-facing → forwarded to OPALX)
# -----------------------------------------------------------------------------

set(PLATFORMS "" CACHE STRING
    "Execution backends: SERIAL, OPENMP, CUDA, HIP")
set_property(CACHE PLATFORMS PROPERTY STRINGS SERIAL OPENMP CUDA HIP)

# Only set default if user didn't specify anything
if(NOT PLATFORMS)
    set(PLATFORMS "SERIAL" CACHE STRING "" FORCE)
endif()

# Forward to OPALX internal variable
set(OPALX_PLATFORMS "${PLATFORMS}" CACHE STRING "" FORCE)

# -----------------------------------------------------------------------------
# platforms we do support
# -----------------------------------------------------------------------------

set(OPALX_SUPPORTED_PLATFORMS "SERIAL;OPENMP;CUDA;HIP")

# === Normalize to uppercase ===
string(TOUPPER "${OPALX_PLATFORMS}" OPALX_PLATFORMS)

# === Declare a HIP profiler option ===
if("HIP" IN_LIST OPALX_PLATFORMS)
  option(OPALX_ENABLE_HIP_PROFILER "Enable HIP Systems Profiler" OFF)
endif()

# -----------------------------------------------------------------------------
# Sanity check for known platforms
# -----------------------------------------------------------------------------

set(unhandled_platforms_ ${OPALX_PLATFORMS})
foreach(platform ${OPALX_SUPPORTED_PLATFORMS})
  if(platform IN_LIST OPALX_PLATFORMS)
    list(REMOVE_ITEM unhandled_platforms_ ${platform})
  endif()
endforeach()

if(NOT unhandled_platforms_ STREQUAL "")
  message(FATAL_ERROR "Unknown or unsupported OPALX_PLATFORMS requested: '${unhandled_platforms_}'")
endif()

if("HIP" IN_LIST OPALX_PLATFORMS AND "CUDA" IN_LIST OPALX_PLATFORMS)
  message(FATAL_ERROR "CUDA and HIP should not both be present in OPALX_PLATFORMS")
endif()

# -----------------------------------------------------------------------------
# Profiler section
# -----------------------------------------------------------------------------

if(OPALX_ENABLE_HIP_PROFILER)
  if("HIP" IN_LIST OPALX_PLATFORMS)
    message(STATUS "🧩 Enabling HIP Profiler and KOKKOS profiliing")
    add_compile_definitions(-DOPALX_ENABLE_HIP_PROFILER)
  else()
    message(FATAL_ERROR "Cannot enable HIP Systems Profiler since platform is not HIP")
  endif()
endif()

if(OPALX_ENABLE_NSYS_PROFILER)
  if("CUDA" IN_LIST OPALX_PLATFORMS)
    message(STATUS "Enabling Nsys Profiler")
    add_compile_definitions(-DOPALX_ENABLE_NSYS_PROFILER)
  else()
    message(FATAL_ERROR "Cannot enable Nvidia Nsys Profiler since platform is not CUDA")
  endif()
endif()

# -----------------------------------------------------------------------------
# Kokkos architecture selection
# -----------------------------------------------------------------------------

set(ARCH "" CACHE STRING
    "Kokkos architecture: AMPERE80, PASCAL61, VOLTA70, BDW, SKX")
set_property(CACHE ARCH PROPERTY STRINGS AMPERE80 PASCAL61 VOLTA70 BDW SKX)

# Clear all architecture flags
foreach(a AMPERE80 PASCAL61 VOLTA70 BDW SKX)
    set(Kokkos_ARCH_${a} OFF CACHE BOOL "" FORCE)
endforeach()

# Enable selected arch (if provided)
if(ARCH)
    set(Kokkos_ARCH_${ARCH} ON CACHE BOOL "" FORCE)
endif()

# -----------------------------------------------------------------------------
# IPPL default features
# -----------------------------------------------------------------------------

set(IPPL_ENABLE_FFT     ON  CACHE BOOL "" FORCE)
set(IPPL_ENABLE_SOLVERS ON  CACHE BOOL "" FORCE)
set(IPPL_ENABLE_ALPINE  OFF CACHE BOOL "" FORCE)
set(IPPL_ENABLE_TESTS   OFF CACHE BOOL "" FORCE)
set(IPPL_PLATFORMS "${OPALX_PLATFORMS}" CACHE STRING "" FORCE)

# -----------------------------------------------------------------------------
# Other OPALX options
# -----------------------------------------------------------------------------

option(BUILD_SHARED_LIBS "Build OPALX as a shared library" OFF)
option(OPALX_ENABLE_EXAMPLES "Enable building the Example module" OFF)
option(OPALX_ENABLE_TESTS "Build integration tests in test/ directory" OFF)
option(OPALX_ENABLE_COVERAGE "Enable code coverage" OFF)
option(OPALX_ENABLE_NSYS_PROFILER "Enable Nvidia Nsys Profiler" OFF)
option(OPALX_ENABLE_SANITIZER "Enable sanitizer(s)" OFF)
option(OPALX_USE_ALTERNATIVE_VARIANT
       "Use modified variant implementation (required for CUDA 12.2 + GCC 12.3.0)" OFF)
option(OPALX_USE_STANDARD_FOLDERS "Put all generated binaries in bin/lib folders" OFF)
option(OPALX_SKIP_FAILING_TESTS "Do not build/test tests that are currently marked as failing" OFF)
option(OPALX_ENABLE_SCRIPTS "Generate job script templates for some benchmarks/tests" OFF)
option(OPALX_FIELD_DEBUG "Enable FieldSolver field-dump debug code" OFF)


if(OPALX_FIELD_DEBUG)
  message(STATUS "⚠️  OPALX_FIELD_DEBUG enabled — field dumps will be emitted during simulation")
endif()

# "Build OPALX as a shared library (ON) or static library (OFF)" OFF) 
if(OPALX_DYL)
    set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
    message(WARNING "OPALX_DYL is deprecated; use -DBUILD_SHARED_LIBS=ON instead.")
endif()