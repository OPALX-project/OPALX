#
# High-level OPALX configuration options
#

# ---------------------------------------------------------
# Build type (user-facing)
# ---------------------------------------------------------
set(BUILD_TYPE "Release" CACHE STRING "Build type: Debug, Release, RelWithDebInfo, MinSizeRel")
set(CMAKE_BUILD_TYPE "${BUILD_TYPE}" CACHE STRING "" FORCE)

# ---------------------------------------------------------
# Platform selection (user-facing → forwarded to OPALX)
# ---------------------------------------------------------
set(PLATFORMS "SERIAL" CACHE STRING
    "Execution backends: SERIAL, OPENMP, CUDA, HIP")
set_property(CACHE PLATFORMS PROPERTY STRINGS SERIAL OPENMP CUDA HIP)

# Forward to the internal variable handled by Platforms.cmake
set(OPALX_PLATFORMS "${PLATFORMS}" CACHE STRING "" FORCE)

# ---------------------------------------------------------
# Kokkos architecture selection
# ---------------------------------------------------------
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

# --- IPPL default features ---
set(IPPL_ENABLE_FFT     ON  CACHE BOOL "" FORCE)
set(IPPL_ENABLE_SOLVERS ON  CACHE BOOL "" FORCE)
set(IPPL_ENABLE_ALPINE  OFF CACHE BOOL "" FORCE)
set(IPPL_ENABLE_TESTS   OFF CACHE BOOL "" FORCE)

# --- OPALX options ---
option(BUILD_SHARED_LIBS "Build OPALX as a shared library" OFF)
option(OPALX_ENABLE_UNIT_TESTS "Enable unit tests using GoogleTest" OFF)
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

# "Build OPALX as a shared library (ON) or static library (OFF)" OFF) 
if(OPALX_DYL)
    set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
    message(WARNING "OPALX_DYL is deprecated; use -DBUILD_SHARED_LIBS=ON instead.")
endif()