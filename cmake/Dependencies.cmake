# -----------------------------------------------------------------------------
# Dependencies.cmake
# ~~~
#
# Resolves third-party libraries: Kokkos and Heffte.
#
# Not responsible for:
#   - Selecting platform backends            → Platforms.cmake
#   - Enabling compiler flags                → CompilerOptions.cmake
#   - Version variables or target creation   → Version.cmake / src/
#
# ~~~
# -----------------------------------------------------------------------------
set(FETCHCONTENT_BASE_DIR "${PROJECT_BINARY_DIR}/_deps")
set(FETCHCONTENT_UPDATES_DISCONNECTED ON) # opt out of auto-updates
set(FETCHCONTENT_QUIET ON)

include(FetchContent)

# ------------------------------------------------------------------------------
# MPI
# ------------------------------------------------------------------------------
find_package(MPI COMPONENTS CXX REQUIRED)
colour_message(STATUS ${Green} "✅ MPI found ${MPI_CXX_VERSION}")

# ------------------------------------------------------------------------------
# CUDA
# ------------------------------------------------------------------------------
if("CUDA" IN_LIST IPPL_PLATFORMS)
  find_package(CUDAToolkit REQUIRED)
  colour_message(STATUS ${Green}
                 "✅ CUDA platform requested and CUDAToolkit found ${CUDAToolkit_VERSION}")
endif()

# ------------------------------------------------------------------------------
# OpenMP
# ------------------------------------------------------------------------------
if("OPENMP" IN_LIST IPPL_PLATFORMS)
  find_package(OpenMP REQUIRED)
  colour_message(STATUS ${Green} "✅ OpenMP platform requested OpenMP found ${OPENMP_VERSION}")
endif()

# ------------------------------------------------------------------------------
# Utility function to clear a list of vars one by one
# ------------------------------------------------------------------------------
function(unset_vars)
  foreach(VAR IN LISTS ARGN)
    unset(${VAR} PARENT_SCOPE)
  endforeach()
endfunction()

# ------------------------------------------------------------------------------
# Utility function to get git tag/sha/version from version string
# ------------------------------------------------------------------------------
function(extract_git_label VERSION_STRING RESULT_VAR)
  if("${${VERSION_STRING}}" MATCHES "^git\\.(.+)$")
    set(${RESULT_VAR} "${CMAKE_MATCH_1}" PARENT_SCOPE)
  else()
    unset(${RESULT_VAR} PARENT_SCOPE)
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Utility function to get git tags from a repo before downloading (unused currently)
# ------------------------------------------------------------------------------
function(get_git_tags GIT_REPOSITORY RESULT_VAR)
  message("Fetching git tags for repo ${GIT_REPOSITORY}")
  execute_process(
    COMMAND git -c versionsort.suffix=- ls-remote --tags --sort=v:refname ${GIT_REPOSITORY}
    COMMAND cut --delimiter=/ --fields=3
    COMMAND grep -Po "^[\\d.]+$"
    OUTPUT_VARIABLE GIT_TAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Convert the output string into a CMake list
  string(REPLACE "\n" ";" GIT_TAGS_LIST "${GIT_TAGS}")
  set(${RESULT_VAR} "${GIT_TAGS_LIST}" PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# IPPL library
# ------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if ("${IPPL_PLATFORMS}" STREQUAL "CUDA")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wno-deprecated-gpu-targets")
endif()

# Disable compile time assert (used by IPPL)
add_definitions (-DNOCTAssert)

# Allow user to specify branch/tag, default to master
set(IPPL_GIT_TAG "master" CACHE STRING "Branch or tag for IPPL (default: master)")
message(STATUS "Fetching IPPL branch/tag: ${IPPL_GIT_TAG}")

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose build type" FORCE)
endif()
message(STATUS "Build type is: ${CMAKE_BUILD_TYPE}")

FetchContent_Declare(
    IPPL
    GIT_REPOSITORY https://github.com/IPPL-framework/ippl.git
    GIT_TAG ${IPPL_GIT_TAG}
    GIT_SHALLOW TRUE
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

FetchContent_MakeAvailable(IPPL)
message(STATUS "IPPL include path: ${IPPL_SOURCE_DIR}/src")

# set(IPPL_INCLUDE_DIR "${IPPL_SOURCE_DIR}/src")
set(IPPL_LIBRARY ippl)

message(STATUS "Found IPPL_DIR: ${IPPL_DIR}")
if(NOT IPPL_VERSION)
    set(IPPL_VERSION "3.2.0")
    message(STATUS "Defaulting to IPPL-${IPPL_VERSION}")
endif()

# ------------------------------------------------------------------------------
# HDF5
# ------------------------------------------------------------------------------
if(OPALX_USE_INSTALLED_HDF5)

    message(STATUS "⚙ Using system-installed HDF5")

    # Require parallel HDF5 with C++ components disabled (you only use C libs)
    find_package(HDF5 REQUIRED COMPONENTS C HL)

    if(NOT HDF5_FOUND)
        message(FATAL_ERROR "System HDF5 requested but not found.")
    endif()

    if (HDF5_VERSION VERSION_LESS "1.10.0")
    message(FATAL_ERROR
        "System HDF5 version ${HDF5_VERSION} is too old. "
        "Required version is >= 1.10.0")
    endif()

    # Normalize your variable names
    set(HDF5_INCLUDE_DIR ${HDF5_INCLUDE_DIRS})
    set(HDF5_LIBRARIES   ${HDF5_LIBRARIES})

    message(STATUS "✔ Found system HDF5 version: ${HDF5_VERSION}")
    message(STATUS "✔ HDF5 include dir: ${HDF5_INCLUDE_DIR}")

else()
    message(STATUS "⚙ Building HDF5 from source (FetchContent)")
    set(HDF5_ENABLE_PARALLEL ON)
    set(HDF5_BUILD_HL_LIB OFF CACHE BOOL “” FORCE) # Disable high-level APIs for thread safety
    set(HDF5_BUILD_EXAMPLES OFF CACHE BOOL “” FORCE) # Disable examples
    set(HDF5_BUILD_TOOLS OFF CACHE BOOL “” FORCE) # Disable tools
    set(HDF5_ENABLE_THREADSAFE OFF CACHE BOOL “” FORCE)
    set(HDF5_TEST_PARALLEL OFF)
    set(HDF5_VERSION "1.14.6")
    set(HDF5_VERSION_MAJOR "1.14")

    set(HDF5_VERSION "1.14.6")
    set(HDF5_ENABLE_FLOAT16 OFF CACHE BOOL "Disable half-precision floats" FORCE)
    FetchContent_Declare(
        HDF5
        URL https://github.com/HDFGroup/hdf5/releases/download/hdf5_${HDF5_VERSION}/hdf5-${HDF5_VERSION}.tar.gz
        URL_HASH "SHA256=e4defbac30f50d64e1556374aa49e574417c9e72c6b1de7a4ff88c4b1bea6e9b"
    )

    # Now FetchContent_MakeAvailable will see the option
    FetchContent_MakeAvailable(HDF5)
    set(HDF5_FOUND TRUE)

    if (TARGET hdf5-shared)
        install(TARGETS hdf5-shared EXPORT ipplTargets DESTINATION lib)
        add_library(hdf5::hdf5 ALIAS hdf5-shared)
    elseif(TARGET hdf5-static)
        install(TARGETS hdf5-static EXPORT ipplTargets DESTINATION lib)
        add_library(hdf5::hdf5 ALIAS hdf5-static)
    endif()

    set(HDF5_LIBRARIES hdf5::hdf5)
    message ("HDF5_FOUND and dir are ${HDF5_FOUND} ${HDF5_BINARY_DIR} and ${HDF5_LIBRARIES} and ${HDF5_VERSION}")
endif()

message(STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")
message(STATUS "HDF5 include dir: ${HDF5_INCLUDE_DIR}")

# ------------------------------------------------------------------------------
# H5hut
# ------------------------------------------------------------------------------
if(OPALX_USE_INSTALLED_H5HUT)
    message(STATUS "⚙ Using system-installed H5hut")

    # If you know the module set these env variables
    if(DEFINED ENV{H5HUT_INCLUDE_DIR} AND DEFINED ENV{H5HUT_LIBRARY_DIR})
        set(H5HUT_INCLUDE_DIRS $ENV{H5HUT_INCLUDE_DIR})
        set(H5HUT_LIBRARIES    $ENV{H5HUT_LIBRARY_DIR}/libh5hut.so)
        message(STATUS "✔ Using H5hut include: ${H5HUT_INCLUDE_DIRS}")
        message(STATUS "✔ Using H5hut library: ${H5HUT_LIBRARIES}")
    endif()

    # Prefer config mode first
    find_package(H5hut QUIET CONFIG)

    if(NOT H5hut_FOUND)
        # fallback to module mode if someone uses FindH5hut.cmake
        find_package(H5hut REQUIRED MODULE)
    endif()

    if(NOT H5hut_FOUND)
        message(FATAL_ERROR "System H5hut requested, but not found.")
    endif()

    message(STATUS "✔ Found system H5Hut")

    # Normalize variables for downstream use
    set(H5HUT_INCLUDE_DIR ${H5HUT_INCLUDE_DIRS} CACHE PATH "H5hut include dir")
    set(H5HUT_LIBRARY     ${H5HUT_LIBRARIES}    CACHE FILEPATH "H5hut library")
else()
    message(STATUS "⚙ Building H5Hut from source (FetchContent)")
    set(H5hut_VERSION cmake)
    set(H5hut_GIT https://github.com/eth-cscs/h5hut.git)
    set(H5hut_WITH_MPI ON)
    set(fetch_string
      GIT_TAG ${H5hut_VERSION}
      GIT_REPOSITORY ${H5hut_GIT})

    # Invoke cmake fetch/find
    FetchContent_Declare(H5hut ${fetch_string})
    FetchContent_MakeAvailable(H5hut)

    # Check that kokkos actually has the platform backends that we need
    if (H5hut_FOUND)
      message(STATUS "H5hut ${H5hut_VERSION} found externally")
    else()
      message(STATUS "H5hut ${H5hut_VERSION} building from source in ${H5hut_SOURCE_DIR}")
      set(H5hut_FOUND ON)
    endif()

    # The H5Hut CMake project itself creates the target "H5hut"
    # and sets these variables:
    #     H5HUT_INCLUDE_DIRS
    #     H5HUT_LIBRARIES
    set(H5HUT_INCLUDE_DIR ${H5HUT_INCLUDE_DIRS})
endif()

# Normalize exported variables used by your project
message(STATUS "H5HUT include dir: ${H5HUT_INCLUDE_DIR}")

# ------------------------------------------------------------------------------
# GoogleTest
# ------------------------------------------------------------------------------
if(OPALX_ENABLE_UNIT_TESTS)
    if(OPALX_USE_INSTALLED_GTEST)
        message(STATUS "⚙ Using system-installed GoogleTest")

        find_package(GTest QUIET)

        # fallback if GTest_FOUND not set, but dirs exist
        if(NOT GTest_FOUND AND DEFINED GTEST_LIBRARIES AND DEFINED GTEST_INCLUDE_DIRS)
            set(GTest_FOUND TRUE)
        endif()

        if(NOT GTest_FOUND)
            message(FATAL_ERROR "OPALX_USE_INSTALLED_GTEST=ON but system GTest not found.")
        endif()

        # Normalize variables
        set(GTest_INCLUDE_DIRS ${GTEST_INCLUDE_DIRS})
        set(GTest_LIBRARIES   ${GTEST_LIBRARIES})

        message(STATUS "✔ Found system GTest: ${GTest_LIBRARIES}")

        # Create imported target for gtest framework
        if(NOT TARGET GTest::gtest)
            add_library(GTest::gtest STATIC IMPORTED GLOBAL)
            set_target_properties(GTest::gtest PROPERTIES
                IMPORTED_LOCATION "${GTEST_LIBRARIES}"
                INTERFACE_INCLUDE_DIRECTORIES "${GTEST_INCLUDE_DIRS}"
            )
        endif()

        # Imported target for gtest_main (depends on gtest!)
        get_filename_component(GTEST_LIB_DIR "${GTEST_LIBRARIES}" DIRECTORY)
        set(GTEST_MAIN_LIB "${GTEST_LIB_DIR}/libgtest_main.a")

        if(EXISTS "${GTEST_MAIN_LIB}")
            if(NOT TARGET GTest::gtest_main)
                add_library(GTest::gtest_main STATIC IMPORTED GLOBAL)
                set_target_properties(GTest::gtest_main PROPERTIES
                    IMPORTED_LOCATION "${GTEST_MAIN_LIB}"
                    INTERFACE_LINK_LIBRARIES "GTest::gtest"   # crucial!
                    INTERFACE_INCLUDE_DIRECTORIES "${GTEST_INCLUDE_DIRS}"
                )
            endif()
            message(STATUS "✔ Found libgtest_main.a: ${GTEST_MAIN_LIB}")
        else()
            message(STATUS "⚠ libgtest_main.a not found; tests will need their own main()")
        endif()



    else()
          message(STATUS "⚙ Building GoogleTest from source (FetchContent)")

          FetchContent_Declare(GTest GIT_REPOSITORY "https://github.com/google/googletest"
                              GIT_TAG "v1.16.0" GIT_SHALLOW ON)

          # For Windows: force shared crt, ignored on linux
          set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

          # Turn off GTest install/tests in the subproject
          set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
          set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
          set(BUILD_GTEST ON CACHE BOOL "" FORCE)

          FetchContent_MakeAvailable(GTest)
          message(STATUS "✅ GoogleTest built from source (${GTest_VERSION})")
    endif()
endif()
# ------------------------------------------------------------------------------
