
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
    "Execution backends: SERIAL, OPENMP, CUDA, HIP, SYCL")
set_property(CACHE PLATFORMS PROPERTY STRINGS SERIAL OPENMP CUDA HIP SYCL)

# Only set default if user didn't specify anything
if(NOT PLATFORMS)
    set(PLATFORMS "SERIAL" CACHE STRING "" FORCE)
endif()

# Forward to OPALX internal variable (internal, derived from user's PLATFORMS)
set(OPALX_PLATFORMS "${PLATFORMS}" CACHE STRING "" FORCE)

# -----------------------------------------------------------------------------
# platforms we do support (internal, set based on user selections)
# -----------------------------------------------------------------------------

set(OPALX_SUPPORTED_PLATFORMS "SERIAL;OPENMP;CUDA;HIP;SYCL")

# === Normalize to uppercase ===
string(TOUPPER "${OPALX_PLATFORMS}" OPALX_PLATFORMS)

# -----------------------------------------------------------------------------
# HIP profiler (user-facing)
#
# Always declare the option so it shows up in cmake-gui/cmake -L,
# but enforce that it can only be enabled for the HIP backend.
# -----------------------------------------------------------------------------
option(OPALX_ENABLE_HIP_PROFILER "Enable HIP Systems Profiler (HIP backend only)" OFF)

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
# Internal: Device compilation flag detection (auto-derived from PLATFORMS)
# 
# Automatically detects if the user selected a device backend (CUDA, HIP, SYCL)
# and sets OPALX_DEVICE_COMPILATION compile definition accordingly.
# This allows code to distinguish device-specific implementations from 
# host-only code. Users do NOT set this directly; it is derived from their 
# PLATFORMS selection.
# -----------------------------------------------------------------------------

set(HOST_ONLY_PLATFORMS "SERIAL;OPENMP;THREADS;HPX")
set(OPALX_HAS_DEVICE_BACKEND FALSE)  # internal: derived variable
foreach(platform ${OPALX_PLATFORMS})
  if(NOT platform IN_LIST HOST_ONLY_PLATFORMS)
    set(OPALX_HAS_DEVICE_BACKEND TRUE)
    break()
  endif()
endforeach()

if(OPALX_HAS_DEVICE_BACKEND)
  add_compile_definitions(OPALX_DEVICE_COMPILATION)
  colour_message(STATUS ${Green} "✅ Device backend detected - defining OPALX_DEVICE_COMPILATION")
else()
  colour_message(STATUS ${Cyan} "ℹ️  Host-only backend (SERIAL/OPENMP) - OPALX_DEVICE_COMPILATION not defined")
endif()

# -----------------------------------------------------------------------------
# Profiler section
# -----------------------------------------------------------------------------

message(STATUS "🔧 HIP profiler (OPALX_ENABLE_HIP_PROFILER): ${OPALX_ENABLE_HIP_PROFILER}")
if(OPALX_ENABLE_HIP_PROFILER)
  if("HIP" IN_LIST OPALX_PLATFORMS)
    message(FATAL_ERROR
      "OPALX_ENABLE_HIP_PROFILER was enabled, but HIP profiling support is currently not implemented in OPALX. "
      "This option is reserved for future work. Please configure with -DOPALX_ENABLE_HIP_PROFILER=OFF for now.")
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
    "Kokkos architecture (CMAKE_CUDA_ARCHITECTURES is auto-configured from this): HOPPER90, AMPERE80, AMPERE86, AMPERE87, VOLTA70, VOLTA72, PASCAL61, PASCAL60, ADA102, BDW, SKX, CNL, ZEN, ZEN2, ZEN3, ARM80, ARM81, INTEL_GEN9, INTEL_GEN11, INTEL_GEN12LP, INTEL_DG1, INTEL_XEHP, INTEL_PVC")
set_property(CACHE ARCH PROPERTY STRINGS "" HOPPER90 AMPERE80 AMPERE86 AMPERE87 VOLTA70 VOLTA72 PASCAL61 PASCAL60 ADA102 BDW SKX CNL ZEN ZEN2 ZEN3 ARM80 ARM81 INTEL_GEN9 INTEL_GEN11 INTEL_GEN12LP INTEL_DG1 INTEL_XEHP INTEL_PVC)

# Clear all architecture flags (include new ones for future-proofing)
foreach(a HOPPER90 AMPERE80 AMPERE86 AMPERE87 VOLTA70 VOLTA72 PASCAL61 PASCAL60 ADA102 BDW SKX CNL ZEN ZEN2 ZEN3 ARM80 ARM81 INTEL_GEN9 INTEL_GEN11 INTEL_GEN12LP INTEL_DG1 INTEL_XEHP INTEL_PVC)
    set(Kokkos_ARCH_${a} OFF CACHE BOOL "" FORCE)
endforeach()

# Enable selected arch (if provided)
if(ARCH)
    set(Kokkos_ARCH_${ARCH} ON CACHE BOOL "" FORCE)
endif()

# -----------------------------------------------------------------------------
# Auto-configure CMAKE_CUDA_ARCHITECTURES from selected ARCH.
# For now, this only works with CUDAs architecture naming convention (which
# is "microarchitecture + compute capability"). At some point, we may want to
# extend this to all architectures that are supported by Kokkos.
# See https://kokkos.org/kokkos-core-wiki/get-started/configuration-guide.html#gpu-architectures
# Note: For non-CUDA backends (HIP, SYCL), the auto-configuration does not apply.
# Users must explicitly specify CMAKE_HIP_ARCHITECTURES or equivalent via command-line.
# -----------------------------------------------------------------------------

if("CUDA" IN_LIST OPALX_PLATFORMS)
    include(ArchitectureMapping)
    configure_cuda_architectures_from_kokkos_arch()
endif()

# Validate architecture configuration for non-CUDA backends
if("HIP" IN_LIST OPALX_PLATFORMS)
    if(NOT DEFINED CMAKE_HIP_ARCHITECTURES OR CMAKE_HIP_ARCHITECTURES STREQUAL "")
        message(FATAL_ERROR 
            "HIP backend selected but CMAKE_HIP_ARCHITECTURES is not set. "
            "Please explicitly specify the target HIP GPU architecture(s) via:\n"
            " -DCMAKE_HIP_ARCHITECTURES=<arch>.\n")
    endif()
    colour_message(STATUS ${Green} "✅ HIP backend: CMAKE_HIP_ARCHITECTURES=${CMAKE_HIP_ARCHITECTURES}")
endif()

if("SYCL" IN_LIST OPALX_PLATFORMS)
    # SYCL architecture configuration is compiler-dependent (Intel DPC++/icpx is
    # the supported toolchain). Use -DARCH=INTEL_PVC (or another INTEL_* arch
    # listed above) to enable the corresponding Kokkos_ARCH_INTEL_* flag.
    colour_message(STATUS ${Yellow} "⚠️  SYCL backend: Architecture configuration is compiler-dependent.")
    colour_message(STATUS ${Yellow} "   Ensure your SYCL compiler has GPU support enabled via environment or compiler flags.")
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
option(OPALX_USE_KOKKOS_MATH_CONSTANTS "Use Kokkos mathematical constants in Physics.h" ON)
option(OPALX_USE_INSTALLED_HDF5 "Use system-installed HDF5 instead of building from source" OFF)
option(OPALX_USE_INSTALLED_H5HUT "Use system-installed H5HUT instead of building from source" OFF)
option(OPALX_USE_INSTALLED_GTEST "Use system-installed GoogleTest instead of building from source" OFF)

# -----------------------------------------------------------------------------
# Compile-definition style options (keep consistent with other global compile definitions)
# -----------------------------------------------------------------------------
if(OPALX_USE_KOKKOS_MATH_CONSTANTS)
  add_compile_definitions(OPALX_USE_KOKKOS_MATH_CONSTANTS)
  colour_message(STATUS ${Green} "✅ Kokkos math constants enabled (OPALX_USE_KOKKOS_MATH_CONSTANTS)")
else()
  colour_message(STATUS ${Cyan} "ℹ️  Kokkos math constants disabled (using literal constants in Physics.h)")
endif()

# Keep high-signal toggles visible in configure output.
message(STATUS "🔧 Skip failing tests (OPALX_SKIP_FAILING_TESTS): ${OPALX_SKIP_FAILING_TESTS}")

# "Build OPALX as a shared library (ON) or static library (OFF)" OFF) 
if(OPALX_DYL)
    set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
    message(WARNING "OPALX_DYL is deprecated; use -DBUILD_SHARED_LIBS=ON instead.")
endif()
