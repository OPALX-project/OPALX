OPTION (OPTP_USED_FROM_OPAL "No OPT-PILOT standalone build" ON)

IF (NOT OPTP_USED_FROM_OPAL)
    cmake_minimum_required (VERSION 2.8.2)
    set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)
    project (opt-pilot)
    set (PACKAGE \"OPT-PILOT\")
    set (PACKAGE_VERSION \"0.1.4\")
ENDIF (NOT OPTP_USED_FROM_OPAL)

#-----------------------------------------------------------------------------
# config.h
IF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    set (GIT_REVISION 0)
    EXECUTE_PROCESS(
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND git describe
            OUTPUT_VARIABLE GIT_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    set (GIT_REVISION "#define GIT_VERSION=\"${GIT_OUT}\"")
ELSE(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    set (GIT_REVISION "")
ENDIF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")

set (BUILD_DATE 0)
EXECUTE_PROCESS(
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMAND date
        OUTPUT_VARIABLE DATE_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
    )
set (BUILD_DATE \"${DATE_OUT}\")


configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

#-----------------------------------------------------------------------------
# dependent packages
find_package (MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

#-----------------------------------------------------------------------------

set(BOOST_ROOT $ENV{BOOST_HOME})
set(Boost_USE_STATIC_LIBS ON)
set(BOOST_NO_SYSTEM_PATHS TRUE)
set(BOOST_INCLUDEDIR $ENV{BOOST_HOME})
set(BOOST_LIBRARYDIR $ENV{BOOST_HOME}}/lib)
find_package (Boost 1.55.0 REQUIRED chrono filesystem regex serialization system iostreams)
if (Boost_INCLUDE_DIRS)
    message (STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
    include_directories (${Boost_INCLUDE_DIRS})
endif ()

#-----------------------------------------------------------------------------
# options

OPTION (ARCH_BG "Set architecture to BlueGene" OFF)

#-----------------------------------------------------------------------------
# Common sources and libs

# add_sources macro implemented by Christof Kraus
# for all common source files (Expressions and Utils)

set (_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
macro (add_sources)
    file (RELATIVE_PATH _relPath "${_SRC_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND OPTP_SRCS "${CMAKE_SOURCE_DIR}/opt-pilot/${_relPath}/${_src}")
        else()
            list (APPEND OPTP_SRCS "${CMAKE_SOURCE_DIR}/opt-pilot/${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SRCS to parent directory
        set (OPTP_SRCS ${OPTP_SRCS} PARENT_SCOPE)
    endif()
endmacro()

set (OPTP_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)

add_subdirectory (Expression)
add_subdirectory (Util)
add_subdirectory (extlib)


#-----------------------------------------------------------------------------
# Build liboptp for linking to OPAL
#
set (CMAKE_CXX_FLAGS
    "${IPPL_CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS}"
    )
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )
add_library ( optp ${OPTP_SRCS} ${PILOT_OPAL_SRCS} )
# Set OPTP variables everywhere (childs and parents)
set (OPTP_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}" CACHE INTERNAL "" FORCE)
set (OPTP_LIBRARY "optp" CACHE INTERNAL "" FORCE)
set (OPTP_LIBS
    pthread
    m
    ${Boost_LIBRARIES}
    wfgHypervolume
    CACHE INTERNAL "" FORCE
)

IF (NOT OPTP_USED_FROM_OPAL)
    #-----------------------------------------------------------------------------
    # Executables
    add_subdirectory (Drivers)
ENDIF (NOT OPTP_USED_FROM_OPAL)

#-----------------------------------------------------------------------------
# Testing

IF (BUILD_OPAL_UNIT_TESTS)
    enable_testing()
    add_subdirectory(Tests)
ENDIF (BUILD_OPAL_UNIT_TESTS)

#-----------------------------------------------------------------------------
# Doxygen
find_package(Doxygen)
if (DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(optpilot-doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif(DOXYGEN_FOUND)