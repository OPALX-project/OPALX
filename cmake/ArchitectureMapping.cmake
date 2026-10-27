# -----------------------------------------------------------------------------
# ArchitectureMapping.cmake
#
# Automatically extracts CUDA compute capability from Kokkos architecture names
# (e.g., HOPPER90, AMPERE80) for use with CMAKE_CUDA_ARCHITECTURES and tools
# that need CUDA architecture info.
#
# This works for any Kokkos architecture with a numeric suffix by extracting
# the trailing digits, so no hardcoded mapping table is needed.
#
# Example:
#   configure_cuda_architectures_from_kokkos_arch()
#   # If ARCH=HOPPER90, sets CMAKE_CUDA_ARCHITECTURES=90
#   # If ARCH=AMPERE80, sets CMAKE_CUDA_ARCHITECTURES=80
#   # If ARCH=BDW or AMD_GFX90A, stops with a backend mismatch error
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Function: opalx_arch_is_nvidia(<out-var>)
#
# Sets <out-var> to TRUE if ARCH names an NVIDIA GPU. Matching the generation
# prefix instead of full names keeps future compute capabilities working.
# -----------------------------------------------------------------------------
function(opalx_arch_is_nvidia OUT_VAR)
    set(_is_nvidia FALSE)
    if(ARCH MATCHES "^(MAXWELL|PASCAL|VOLTA|TURING|AMPERE|ADA|HOPPER|BLACKWELL)[0-9]+$")
        set(_is_nvidia TRUE)
    endif()
    set(${OUT_VAR} ${_is_nvidia} PARENT_SCOPE)
endfunction()

# -----------------------------------------------------------------------------
# Function: configure_cuda_architectures_from_kokkos_arch
#
# Automatically sets CMAKE_CUDA_ARCHITECTURES based on the ARCH cache variable
# when CUDA is in the enabled platforms. Uses regex to extract numeric suffix
# from architecture names (e.g., HOPPER90 -> 90, AMPERE80 -> 80).
#
# This function:
#   - Checks if CUDA is enabled in OPALX_PLATFORMS
#   - Extracts trailing digits from the ARCH value using regex
#   - Stops with an error if ARCH is not an NVIDIA GPU (AMD, Intel or CPU name)
#   - Sets CMAKE_CUDA_ARCHITECTURES from the numeric suffix
# -----------------------------------------------------------------------------
function(configure_cuda_architectures_from_kokkos_arch)
    # Check if CUDA is enabled
    if(NOT DEFINED OPALX_PLATFORMS)
        message(STATUS "ℹ️  OPALX_PLATFORMS not defined; CUDA architecture auto-configuration skipped")
        return()
    endif()

    # Only auto-configure if CUDA is in the enabled platforms
    if(NOT "CUDA" IN_LIST OPALX_PLATFORMS)
        message(STATUS "ℹ️  CUDA not in OPALX_PLATFORMS; CMAKE_CUDA_ARCHITECTURES not auto-configured")
        return()
    endif()

    if(NOT DEFINED ARCH)
        message(STATUS "ℹ️  ARCH not defined; CMAKE_CUDA_ARCHITECTURES will not be auto-configured")
        return()
    endif()

    if(NOT ARCH)
        message(STATUS "ℹ️  ARCH is empty; CMAKE_CUDA_ARCHITECTURES will not be auto-configured")
        return()
    endif()

    # Guard against silently deriving a bogus compute capability from a
    # non-NVIDIA architecture: AMD_GFX1152 with PLATFORMS=CUDA would set
    # CMAKE_CUDA_ARCHITECTURES=1152 and a CPU name like ZEN3 would set 3.
    opalx_arch_is_nvidia(_opalx_arch_is_nvidia)
    if(NOT _opalx_arch_is_nvidia)
        message(FATAL_ERROR
            "❌ CUDA is enabled in OPALX_PLATFORMS, but ARCH=${ARCH} is not an NVIDIA GPU "
            "architecture. Combining PLATFORMS=CUDA with an AMD, CPU, or Intel architecture "
            "is a backend mismatch. Use an NVIDIA ARCH such as HOPPER90, AMPERE80, or "
            "BLACKWELL120, or select a matching backend with -DPLATFORMS=HIP/OPENMP.")
    endif()

    # Extract trailing digits from ARCH using regex
    # This works for HOPPER90, AMPERE80, VOLTA70, etc.
    string(REGEX MATCH "[0-9]+$" CUDA_CAPABILITY "${ARCH}")

    if(CUDA_CAPABILITY)
        set(CMAKE_CUDA_ARCHITECTURES "${CUDA_CAPABILITY}" CACHE STRING
            "CUDA compute capability (auto-extracted from ARCH=${ARCH})" FORCE)
        message(STATUS "✅ Extracted CUDA compute capability from ARCH=${ARCH} -> CMAKE_CUDA_ARCHITECTURES=${CUDA_CAPABILITY}")
    else()
        message(FATAL_ERROR 
            "❌ CUDA is enabled in OPALX_PLATFORMS, but ARCH=${ARCH} has no numeric suffix to indicate compute capability.\n"
            "   When CUDA is enabled, you must specify a GPU architecture with compute capability like HOPPER90, AMPERE80, VOLTA70, etc.\n"
            "   If you don't want CUDA, use e.g. -DPLATFORMS=OPENMP.")
    endif()
endfunction()
