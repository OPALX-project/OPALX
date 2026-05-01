# -----------------------------------------------------------------------------
# src/CMakeLists.txt
# cmake-format: off
#
# This file builds the main OPALX target. It:
#   - Generates BuildInfo.h from current build system info
#   - Gathers sources and headers
#   - Adds all component subdirectories
#   - Links to external libraries
#   - Installs headers and the library
#
# cmake-format: on
# -----------------------------------------------------------------------------
message(STATUS "📄 Generating BuildInfo.h")

# Make embedded metadata optional (reduces rebuild churn).
option(OPALX_EMBED_BUILD_METADATA "Embed user/machine/date in BuildInfo.h" OFF)
message(STATUS "🔧 OPALX_EMBED_BUILD_METADATA: ${OPALX_EMBED_BUILD_METADATA}")

if(OPALX_EMBED_BUILD_METADATA)
  message(STATUS "🧩 Embedding user/machine/date into BuildInfo.h (this may increase rebuild churn)")

  execute_process(COMMAND date OUTPUT_VARIABLE OPALX_COMPILE_DATE OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND uname -a OUTPUT_VARIABLE OPALX_COMPILE_MACHINE
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND whoami OUTPUT_VARIABLE OPALX_COMPILE_USER OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  message(STATUS "🧩 Not embedding build metadata into BuildInfo.h")
  set(OPALX_COMPILE_DATE "unknown")
  set(OPALX_COMPILE_MACHINE "unknown")
  set(OPALX_COMPILE_USER "unknown")
endif()

set(OPALX_COMPILE_ARCH "${CMAKE_SYSTEM_PROCESSOR}")
set(OPALX_COMPILE_LINE "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
set(OPALX_COMPILE_OPTIONS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
set(OPALX_BUGREPORT "https://github.com/OPALX-project/OPALX/issues")
set(OPALX_VERSION_MAJOR_VALUE "${PROJECT_VERSION_MAJOR}")
set(OPALX_VERSION_MINOR_VALUE "${PROJECT_VERSION_MINOR}")
if(PROJECT_VERSION_PATCH)
  set(OPALX_VERSION_PATCH_VALUE "${PROJECT_VERSION_PATCH}")
else()
  set(OPALX_VERSION_PATCH_VALUE "0")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/BuildInfo.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/BuildInfo.h @ONLY)


# ------------------------------------------------------------------
# Build the main executable target: opalx
# ------------------------------------------------------------------

add_library(opalx
        AbsBeamline/MultipoleTConfig.h)
target_sources(opalx PRIVATE changes.cpp)

target_compile_features(opalx PUBLIC cxx_std_20)

target_compile_options(
  opalx
  PUBLIC $<$<COMPILE_LANGUAGE:CXX>:
         $<$<CONFIG:Debug>:-O0
         -g>
         $<$<CONFIG:RelWithDebInfo>:-O3
         -g>
         $<$<CONFIG:Release>:-O3>
         >)

get_target_property(_opalx_type opalx TYPE)
if(_opalx_type)
  string(REPLACE "_LIBRARY" "" _opalx_type "${_opalx_type}") # -> STATIC / SHARED / OBJECT / INTERFACE
else()
  set(_opalx_type "UNKNOWN")
endif()

# Query PIC state (may inherit from CMAKE_POSITION_INDEPENDENT_CODE)
get_target_property(_opalx_pic opalx POSITION_INDEPENDENT_CODE)
if(NOT DEFINED _opalx_pic)
  set(_opalx_pic "DEFAULT(${CMAKE_POSITION_INDEPENDENT_CODE})")
endif()

message(STATUS "🔧 Target 'opalx': type=${_opalx_type}, POSITION_INDEPENDENT_CODE=${_opalx_pic}")

target_include_directories(
  opalx PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
              $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/opalx>
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR})


target_sources(opalx PUBLIC changes.cpp)
#
add_subdirectory (AbsBeamline)
add_subdirectory (AbstractObjects)
add_subdirectory (Algorithms)
add_subdirectory (Attributes)
add_subdirectory (BasicActions)
add_subdirectory (BeamlineCore)
add_subdirectory (BeamlineGeometry)
add_subdirectory (Beamlines)
add_subdirectory (Channels)
add_subdirectory (Distribution)
add_subdirectory (Elements)
add_subdirectory (Expressions)
add_subdirectory (Fields)
add_subdirectory (Filters)
add_subdirectory (Lines)
add_subdirectory (OpalConfigure)
add_subdirectory (OpalParser)
add_subdirectory (PartBunch)
add_subdirectory (Physics)
add_subdirectory (Processes)
add_subdirectory (Steppers)
add_subdirectory (Structure)
add_subdirectory (Tables)
add_subdirectory (Track)
add_subdirectory (Utilities)
add_subdirectory (ValueDefinitions)

include(${PROJECT_SOURCE_DIR}/cmake/PlatformOptions.cmake)

target_link_libraries(opalx PUBLIC
    MPI::MPI_CXX
    IPPL::ippl
    H5hut
    ${HDF5_LIBRARIES}
    ZLIB::ZLIB
)

target_include_directories(opalx PRIVATE
  ${HDF5_INCLUDE_DIR}
  ${H5HUT_INCLUDE_DIR}
  ZLIB::ZLIB
)

# CLI executable
add_executable(opalx_exe Main.cpp)
target_link_libraries(opalx_exe PUBLIC opalx)
target_include_directories(opalx_exe PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# this alias should be created after all target properties have been set
add_library(OPALX::opalx ALIAS opalx)

include(InstallOPALX)
