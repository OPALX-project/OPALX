# -----------------------------------------------------------------------------
# src/CMakeLists.txt
# cmake-format: off
#
# This file builds the main OPALX target. It:
#   - Generates OPALXVersions.h from current build system info
#   - Gathers sources and headers
#   - Adds all component subdirectories
#   - Links to external libraries
#   - Installs headers and the library
#
# cmake-format: on
# -----------------------------------------------------------------------------
message(STATUS "📄 Generating OPALXVersions.h")

# Make embedded metadata optional (reduces rebuild churn).
option(OPALX_EMBED_BUILD_METADATA "Embed user/machine/date in OPALXVersions.h" OFF)

if(OPALX_EMBED_BUILD_METADATA)

  execute_process(COMMAND date OUTPUT_VARIABLE OPALX_COMPILE_DATE OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND uname -a OUTPUT_VARIABLE OPALX_COMPILE_MACHINE
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND whoami OUTPUT_VARIABLE OPALX_COMPILE_USER OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  set(OPALX_COMPILE_DATE "unknown")
  set(OPALX_COMPILE_MACHINE "unknown")
  set(OPALX_COMPILE_USER "unknown")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/OPALXVersions.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/OPALXVersions.h @ONLY)


# ------------------------------------------------------------------
# Build the main executable target: opalx
# ------------------------------------------------------------------

add_executable(opalx)
target_sources(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Main.cpp)


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

target_sources(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Main.cpp)

target_include_directories(
  opalx PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
              $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/opalx>
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR})


#add_opal_sources (changes.cpp)
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
add_subdirectory (Steppers)
add_subdirectory (Structure)
add_subdirectory (Tables)
add_subdirectory (Track)
add_subdirectory (Utilities)
add_subdirectory (ValueDefinitions)
add_subdirectory (MemoryManagement)

include(${PROJECT_SOURCE_DIR}/cmake/PlatformOptions.cmake)

target_link_libraries(opalx PUBLIC
    MPI::MPI_CXX
    ippl
    ${OPAL_FFTW_LIBS}
    ${GSL_LIBRARY}
    ${GSL_CBLAS_LIBRARY}
    ${H5Hut_LIBRARY}
    ${HDF5_LIBRARIES}
    ${Boost_LIBRARIES})

target_include_directories(opalx PRIVATE
  ${FFTW_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR}
  ${H5HUT_INCLUDE_DIR}
  ${GSL_INCLUDE_DIR}
  "${ippl_SOURCE_DIR}/src"
  ${BOOST_INCLUDE_DIR}
  ${BOOST_SRC_DIR}/boost_${BOOST_VERSION_UNDERSCORE}  # raw source headers
)

include(InstallOPALX)