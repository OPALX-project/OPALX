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
#   # If ARCH=BDW (CPU only), skips CUDA configuration
# -----------------------------------------------------------------------------

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
#   - Sets CMAKE_CUDA_ARCHITECTURES if ARCH has a numeric suffix
#   - Handles CPU-only architectures gracefully
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
