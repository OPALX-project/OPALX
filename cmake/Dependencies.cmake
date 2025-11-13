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

# set(HDF5_INCLUDE_DIRS ${HDF5_BINARY_DIR}/src)
set(HDF5_LIBRARIES hdf5::hdf5)
message ("HDF5_FOUND and dir are ${HDF5_FOUND} ${HDF5_BINARY_DIR} and ${HDF5_LIBRARIES} and ${HDF5_VERSION}")

# ------------------------------------------------------------------------------
# H5hut
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Boost library
# ------------------------------------------------------------------------------

include(ExternalProject)

set(BOOST_VERSION "1.84.0")
set(BOOST_INCLUDE_LIBRARIES filesystem system optional regex iostreams)
set(BOOST_ENABLE_CMAKE ON)

message(STATUS "Downloading and extracting boost library sources. This will take some time...")
include(FetchContent)
Set(FETCHCONTENT_QUIET FALSE) # Needed to print downloading progress
FetchContent_Declare(
    Boost
    #URL https://github.com/boostorg/boost/releases/download/boost-${BOOST_VERSION}/boost-${BOOST_VERSION}.7z # downloading a zip release speeds up the download
    URL https://boostorg.jfrog.io/artifactory/main/release/${BOOST_VERSION}/source/boost_${BOOST_VERSION//./_}.tar.gz
    USES_TERMINAL_DOWNLOAD TRUE
    GIT_PROGRESS TRUE
    DOWNLOAD_NO_EXTRACT FALSE
)
FetchContent_MakeAvailable(Boost)
set(Boost_LIBRARIES Boost::filesystem Boost::optional Boost::iostreams)
message(STATUS "Boost include dir: ${Boost_INCLUDE_DIR}")
message(STATUS "Boost libraries: ${Boost_LIBRARIES}")

# Header-only Boost interface for uBLAS and other headers
add_library(boost_HEADERS INTERFACE)
target_include_directories(boost_HEADERS INTERFACE
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/numeric/ublas/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/serialization/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/bimap/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/multi_index/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/variant/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/type_index/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/geometry/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/assign/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/algorithm/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/numeric/conversion/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/numeric/odeint/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/units/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/format/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/spirit/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/phoenix/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/proto/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/endian/include
    ${CMAKE_BINARY_DIR}/_deps/boost-src/libs/math/include
)

# ------------------------------------------------------------------------------
# GSL
# ------------------------------------------------------------------------------

set(GSL_VERSION 2.7)
set(GSL_TAR "gsl-${GSL_VERSION}.tar.gz")
set(GSL_URL "http://ftp.gnu.org/gnu/gsl/${GSL_TAR}")
set(GSL_INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/gsl)
set(GSL_SRC_DIR ${CMAKE_BINARY_DIR}/_deps/gsl-src)

include(ExternalProject)

ExternalProject_Add(gsl_external
    URL ${GSL_URL}
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    PREFIX ${GSL_INSTALL_DIR}
    CONFIGURE_COMMAND
        <SOURCE_DIR>/configure
            --prefix=<INSTALL_DIR>
            --disable-shared
            --enable-static
            CFLAGS=-fPIC
    BUILD_COMMAND make -j${NPROC}
    INSTALL_COMMAND make install
    BUILD_IN_SOURCE TRUE
)

ExternalProject_Get_Property(gsl_external INSTALL_DIR)
set(GSL_ROOT ${INSTALL_DIR})
set(GSL_INCLUDE_DIR ${GSL_ROOT}/include)
set(GSL_LIBRARIES
    ${GSL_ROOT}/lib/libgsl.a
    ${GSL_ROOT}/lib/libgslcblas.a
)
message(STATUS "GSL include dir: ${GSL_INCLUDE_DIR}")
message(STATUS "GSL libraries: ${GSL_LIBRARIES}")

# ------------------------------------------------------------------------------
# GoogleTest
# ------------------------------------------------------------------------------
if(OPALX_ENABLE_UNIT_TESTS)
  find_package(GTest)
  if(NOT GTest_FOUND)
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
