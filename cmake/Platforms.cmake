# -----------------------------------------------------------------------------
# Platforms.cmake
# ~~~
#
# Handles platform/backend selection for OPALX: SERIAL, OPENMP, CUDA, or CUDA + OPENMP.
#
# Responsibilities:
#   - Set OPALX_PLATFORMS to default (SERIAL) if unset
#   - Normalize and validate the value
#   - Enable relevant Kokkos/Heffte options
#
# Not responsible for:
#   - Selecting default CMake build type    → ProjectSetup.cmake
#
# ~~~
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# platforms we do support
# -----------------------------------------------------------------------------
set(OPALX_SUPPORTED_PLATFORMS "SERIAL;OPENMP;CUDA;HIP")

# === Default to SERIAL if OPALX_PLATFORMS not set ===
if(NOT OPALX_PLATFORMS)
  set(OPALX_PLATFORMS "SERIAL")
  message(STATUS "No OPALX_PLATFORMS specified — defaulting to SERIAL")
endif()

# === Normalize to uppercase ===
string(TOUPPER "${OPALX_PLATFORMS}" OPALX_PLATFORMS)

# === Declare a HIP profiler option ===
if("HIP" IN_LIST OPALX_PLATFORMS)
  option(OPALX_ENABLE_HIP_PROFILER "Enable HIP Systems Profiler" OFF)
endif()

if(NOT "SERIAL" IN_LIST OPALX_PLATFORMS AND NOT "OPENMP" IN_LIST OPALX_PLATFORMS)
  list(APPEND OPALX_PLATFORMS "SERIAL")
  message(STATUS "Appending SERIAL to OPALX_PLATFORMS as no HOST execution space set")
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
# ---------------------------------------------------------------------
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