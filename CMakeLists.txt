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
if(POLICY CMP0155)
    cmake_policy(SET CMP0155 OLD)
endif()

# ------------------------------------------------------------------------------
# Project declaration and version
# ------------------------------------------------------------------------------
project(OPALX VERSION 0.9 LANGUAGES C CXX)

# ------------------------------------------------------------------------------
# Module path, setup before including any cmake modules
# ------------------------------------------------------------------------------
list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" "${PROJECT_SOURCE_DIR}/CMakeModules")

# ------------------------------------------------------------------------------
# load message macros (for colors), needs to be done here to get the order right
# ------------------------------------------------------------------------------
include(Messages)

# ------------------------------------------------------------------------------
# Primary OPALX options
# ------------------------------------------------------------------------------
include(OPALXOptions)

# ------------------------------------------------------------------------------
# Project setup
# ------------------------------------------------------------------------------
include(ProjectSetup)

# ------------------------------------------------------------------------------
# Setup Output directories
# ------------------------------------------------------------------------------
if(OPALX_USE_STANDARD_FOLDERS)
  message(STATUS "🔧 OPALX_USE_STANDARD_FOLDERS=ON: using ${PROJECT_BINARY_DIR}/bin and ${PROJECT_BINARY_DIR}/lib for all targets")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin
      CACHE PATH "Single Directory for all Executables.")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
      CACHE PATH "Single Directory for all Libraries")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
      CACHE PATH "Single Directory for all static libraries.")
else()
  message(STATUS "🔧 OPALX_USE_STANDARD_FOLDERS=OFF: using CMake default output directories")
  unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY CACHE)
  unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY CACHE)
  unset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY CACHE)
endif()

# ------------------------------------------------------------------------------                                                                                                                                    
# Make OPAL revision
# Run gen_OPALXrevision script at configure time
# ------------------------------------------------------------------------------
execute_process(
    COMMAND ${CMAKE_SOURCE_DIR}/tools/gen_OPALXrevision
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    RESULT_VARIABLE REVISION_RESULT
)

if(NOT REVISION_RESULT EQUAL 0)
    message(FATAL_ERROR "gen_OPALXrevision failed with exit code ${REVISION_RESULT}")
else()
    colour_message(STATUS ${Green} "✅ gen_OPALXrevision completed successfully.")
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

if (OPALX_ENABLE_UNIT_TESTS)
    add_subdirectory (unit_tests)
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
