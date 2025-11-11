# ------------------------------------------------------------------------------
# CMake version requirement
# ------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.24)

# ------------------------------------------------------------------------------
# Policies - use the latest of everything
# ------------------------------------------------------------------------------
cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})

# but disable c++20 module scanning until we have better support for it this seems to cause some
# problems when using gcc@14
cmake_policy(SET CMP0155 OLD)

# ------------------------------------------------------------------------------
# Project declaration and version
# ------------------------------------------------------------------------------
project(OPALX VERSION 1.0 LANGUAGES CXX)

# ------------------------------------------------------------------------------
# Module path, setup before including any cmake modules
# ------------------------------------------------------------------------------
list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" "${PROJECT_SOURCE_DIR}/CMakeModules")

# ------------------------------------------------------------------------------
# Project setup
# ------------------------------------------------------------------------------
include(ProjectSetup)

# ------------------------------------------------------------------------------
# Primary OPALX options
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# Setup Output directories
# ------------------------------------------------------------------------------
if(OPALX_USE_STANDARD_FOLDERS)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin
      CACHE PATH "Single Directory for all Executables.")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
      CACHE PATH "Single Directory for all Libraries")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
      CACHE PATH "Single Directory for all static libraries.")
else()
  unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY CACHE)
  unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY CACHE)
  unset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY CACHE)
endif()

# ------------------------------------------------------------------------------
# cmake modules
# ------------------------------------------------------------------------------
# cmake supplied modules
include(GNUInstallDirs)
include(FetchContent)
include(CMakePackageConfigHelpers)
include(ExternalProject)
include(CTest)

# custom modules
include(Messages)
include(CompilerOptions)
include(Platforms)
include(Dependencies)
include(FailingTests)

# ------------------------------------------------------------------------------
# Deprecated vars
# ------------------------------------------------------------------------------
if(DEFINED USE_ALTERNATIVE_VARIANT)
  colour_message(
    WARNING ${Red}
    "USE_ALTERNATIVE_VARIANT is deprecated. Please set OPALX_USE_ALTERNATIVE_VARIANT instead.")
  set(OPALX_USE_ALTERNATIVE_VARIANT ${USE_ALTERNATIVE_VARIANT} CACHE BOOL "" FORCE)
endif()

# ------------------------------------------------------------------------------
# Debug: This tells the compiler to replace occurrences of ${<DIR>} with <x> in debug info and error
# messages.
# ------------------------------------------------------------------------------
add_compile_options(
  $<$<COMPILE_LANG_AND_ID:C,Clang,AppleClang,GNU>:-ffile-prefix-map=${CMAKE_SOURCE_DIR}=.>
  $<$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang,GNU>:-ffile-prefix-map=${CMAKE_SOURCE_DIR}=.>)

if(DEFINED FETCHCONTENT_BASE_DIR)
  add_compile_options(
    $<$<COMPILE_LANG_AND_ID:C,Clang,AppleClang,GNU>:-ffile-prefix-map=${FETCHCONTENT_BASE_DIR}=.3p>
    $<$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang,GNU>:-ffile-prefix-map=${FETCHCONTENT_BASE_DIR}=.3p>
  )
endif()

# ------------------------------------------------------------------------------
# Define sources for project
# ------------------------------------------------------------------------------
add_subdirectory(src)

include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${OPALX_INCLUDE_DIR}
)

if (BUILD_OPAL_UNIT_TESTS)
    add_subdirectory (tests)
endif ()

# ------------------------------------------------------------------------------
# Platform-specific options
# ------------------------------------------------------------------------------

message (STATUS "Project Name ${PROJECT_NAME}")
# Print relevant paths
message(STATUS "IPPL_INC: ${IPPL_INCLUDE_DIR}")
message(STATUS "IPPL_LIB: ${IPPL_LIBRARY}")
message (STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")
message (STATUS "HDF5 include dir: ${HDF5_INCLUDE_DIR}")
message (STATUS "H5HUT libraries: ${H5HUT_LIBRARIES}")
message (STATUS "H5HUT include dir: ${H5HUT_INCLUDE_DIR}")
message (STATUS "gsl libraries: ${GSL_LIBRARIES}")
message (STATUS "gsl include dirs: ${GSL_INCLUDE_DIR}")
message (STATUS "Boost headers: ${boost_HEADERS}")
message (STATUS "Boost include dirs: ${BOOST_INCLUDE_DIR}")
