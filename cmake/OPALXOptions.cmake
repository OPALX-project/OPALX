
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
# Honour OPALX_PLATFORMS if it was already set (e.g. via preset or -D)
if(NOT PLATFORMS AND OPALX_PLATFORMS)
    set(PLATFORMS "${OPALX_PLATFORMS}" CACHE STRING "" FORCE)
endif()
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
#
# ARCH switches on Kokkos_ARCH_<ARCH>. Kokkos_ARCH_<name> flags given directly
# (the alps-* presets set Kokkos_ARCH_HOPPER90 or Kokkos_ARCH_AMD_GFX942_APU)
# reach Kokkos unchanged. Kokkos ignores names it does not know and then picks
# the architecture by autodetection, so ARCH is checked against the names of
# Kokkos 5.2 (kokkos/cmake/kokkos_arch.cmake).
# -----------------------------------------------------------------------------

set(OPALX_KOKKOS_ARCHITECTURES
    # NVIDIA GPUs
    MAXWELL50 MAXWELL52 MAXWELL53 PASCAL60 PASCAL61 VOLTA70 VOLTA72 TURING75
    AMPERE80 AMPERE86 AMPERE87 ADA89 HOPPER90
    BLACKWELL100 BLACKWELL103 BLACKWELL120 BLACKWELL121
    # AMD GPUs
    AMD_GFX906 AMD_GFX908 AMD_GFX90A AMD_GFX940 AMD_GFX942 AMD_GFX942_APU AMD_GFX950
    AMD_GFX1030 AMD_GFX1100 AMD_GFX1101 AMD_GFX1103 AMD_GFX1151 AMD_GFX1152 AMD_GFX1201
    VEGA906 VEGA908 VEGA90A NAVI1030
    # Intel GPUs
    INTEL_GEN INTEL_GEN9 INTEL_GEN11 INTEL_GEN12LP INTEL_DG1 INTEL_DG2 INTEL_XEHP INTEL_PVC
    # CPUs
    NATIVE AMDAVX ZEN ZEN2 ZEN3 ZEN4 ZEN5 SNB HSW BDW SKL SKX ICL ICX SPR KNC KNL
    ARMV80 ARMV81 ARMV84 ARMV84_SVE ARMV8_THUNDERX ARMV8_THUNDERX2 A64FX ARMV9_GRACE
    POWER8 POWER9 RISCV_SG2042 RISCV_RVA22V RISCV_U74MC)

set(ARCH "" CACHE STRING
    "Kokkos architecture. CUDA architectures auto-configure CMAKE_CUDA_ARCHITECTURES; AMD_GFX architectures auto-configure CMAKE_HIP_ARCHITECTURES for HIP.")
set_property(CACHE ARCH PROPERTY STRINGS "" ${OPALX_KOKKOS_ARCHITECTURES})

# Honour Kokkos_ARCH_FLAG if ARCH was not set (e.g. via CI -D)
if(NOT ARCH AND Kokkos_ARCH_FLAG)
    string(REGEX REPLACE "^Kokkos_ARCH_" "" ARCH "${Kokkos_ARCH_FLAG}")
    set(ARCH "${ARCH}" CACHE STRING "" FORCE)
endif()

# A raw compiler target such as gfx90a is never a Kokkos name and stops here. An
# unknown upper-case name only warns, because a newer installed Kokkos may know it.
if(ARCH AND NOT ARCH IN_LIST OPALX_KOKKOS_ARCHITECTURES)
    set(_opalx_arch_doc "https://kokkos.org/kokkos-core-wiki/get-started/configuration-guide.html")
    if(ARCH MATCHES "^[A-Z][A-Z0-9_]*$")
        message(WARNING
            "ARCH=${ARCH} is not a Kokkos 5.2 architecture. It is passed on as "
            "Kokkos_ARCH_${ARCH}=ON, which Kokkos ignores if it does not know the name. "
            "See ${_opalx_arch_doc}")
    else()
        set(_opalx_arch_hint "")
        if(ARCH MATCHES "^gfx([0-9A-Fa-f]+)$")
            string(TOUPPER "${CMAKE_MATCH_1}" _opalx_gfx_id)
            set(_opalx_arch_hint " For ${ARCH} use ARCH=AMD_GFX${_opalx_gfx_id}.")
        endif()
        message(FATAL_ERROR
            "ARCH='${ARCH}' is not a Kokkos architecture name such as HOPPER90, AMD_GFX90A "
            "or ZEN3.${_opalx_arch_hint} See ${_opalx_arch_doc}")
    endif()
endif()

# ARCH switches its own flag on and the flag of an earlier ARCH off again; other
# Kokkos_ARCH_ flags are left as they were given.
if(OPALX_KOKKOS_ARCH_ENABLED AND NOT OPALX_KOKKOS_ARCH_ENABLED STREQUAL "${ARCH}")
    set(Kokkos_ARCH_${OPALX_KOKKOS_ARCH_ENABLED} OFF CACHE BOOL "" FORCE)
endif()
if(ARCH)
    set(Kokkos_ARCH_${ARCH} ON CACHE BOOL "" FORCE)
endif()
set(OPALX_KOKKOS_ARCH_ENABLED "${ARCH}" CACHE INTERNAL "Kokkos_ARCH_ flag switched on through ARCH")

# -----------------------------------------------------------------------------
# Compiler GPU targets derived from ARCH
#
# CUDA: CMAKE_CUDA_ARCHITECTURES is the compute capability of ARCH (HOPPER90 -> 90).
# HIP:  AMD_GFX<id> and AMD_GFX<id>_APU set CMAKE_HIP_ARCHITECTURES=gfx<id>
#       (AMD_GFX942_APU -> gfx942), unless the user set it, e.g. to gfx942:xnack+.
# SYCL: the target comes from the compiler setup (see below).
# -----------------------------------------------------------------------------

include(ArchitectureMapping)

if("CUDA" IN_LIST OPALX_PLATFORMS)
    configure_cuda_architectures_from_kokkos_arch()
endif()

if("HIP" IN_LIST OPALX_PLATFORMS)
    opalx_arch_is_nvidia(_opalx_arch_is_nvidia)
    if(_opalx_arch_is_nvidia OR ARCH MATCHES "^INTEL_")
        message(FATAL_ERROR
            "❌ HIP is enabled in OPALX_PLATFORMS, but ARCH=${ARCH} is not an AMD GPU "
            "architecture. Use an AMD ARCH such as AMD_GFX90A or AMD_GFX942_APU, or select "
            "a matching backend with -DPLATFORMS=CUDA/SYCL.")
    endif()

    # A value derived from an earlier ARCH follows ARCH; any other value is the user's.
    if(ARCH MATCHES "^AMD_GFX([0-9A-Fa-f]+)(_APU)?$")
        string(TOLOWER "gfx${CMAKE_MATCH_1}" _opalx_hip_architecture)
        if(NOT CMAKE_HIP_ARCHITECTURES
           OR CMAKE_HIP_ARCHITECTURES STREQUAL "${OPALX_DERIVED_HIP_ARCHITECTURES}")
            set(CMAKE_HIP_ARCHITECTURES "${_opalx_hip_architecture}" CACHE STRING
                "HIP GPU architecture (auto-derived from ARCH=${ARCH})" FORCE)
            set(OPALX_DERIVED_HIP_ARCHITECTURES "${_opalx_hip_architecture}" CACHE INTERNAL
                "CMAKE_HIP_ARCHITECTURES derived from ARCH")
        endif()
    endif()

    if(NOT DEFINED CMAKE_HIP_ARCHITECTURES OR CMAKE_HIP_ARCHITECTURES STREQUAL "")
        message(FATAL_ERROR 
            "HIP backend selected but CMAKE_HIP_ARCHITECTURES is not set. "
            "Please explicitly specify the target HIP GPU architecture(s) via:\n"
            " -DARCH=AMD_GFX1152 or -DCMAKE_HIP_ARCHITECTURES=gfx1152.\n")
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
set(IPPL_PLATFORMS "${OPALX_PLATFORMS}" CACHE STRING "" FORCE)
set(IPPL_USE_STANDARD_FOLDERS ${OPALX_USE_STANDARD_FOLDERS})

# -----------------------------------------------------------------------------
# Barnes-Hut space-charge solver (FIELDSOLVER TYPE=BH)
#
# Builds IPPL's gridless N-body module (cstone + ryoanji) and the OPALX
# BarnesHutAlgorithm. The module selects its backend from IPPL_PLATFORMS
# (CUDA -> GPU, otherwise CPU/OpenMP); SYCL is not supported by cstone/ryoanji.
# -----------------------------------------------------------------------------
option(OPALX_ENABLE_BH "Enable the Barnes-Hut space-charge solver (IPPL NBody module)" OFF)
message(STATUS "🔧 Barnes-Hut space-charge solver (OPALX_ENABLE_BH): ${OPALX_ENABLE_BH}")
if(OPALX_ENABLE_BH)
  if("SYCL" IN_LIST OPALX_PLATFORMS)
    message(FATAL_ERROR "OPALX_ENABLE_BH is not supported with the SYCL platform.")
  endif()
  if("HIP" IN_LIST OPALX_PLATFORMS)
    message(FATAL_ERROR "OPALX_ENABLE_BH is untested with HIP; use CUDA or a host platform.")
  endif()
  set(IPPL_ENABLE_NBODY ON CACHE BOOL "" FORCE)
endif()

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

# "Build OPALX as a shared library (ON) or static library (OFF)" OFF) 
if(OPALX_DYL)
    set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
    message(WARNING "OPALX_DYL is deprecated; use -DBUILD_SHARED_LIBS=ON instead.")
endif()
