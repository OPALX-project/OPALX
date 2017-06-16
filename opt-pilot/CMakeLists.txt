IF (NOT OPTP_USED_FROM_OPAL)

cmake_minimum_required (VERSION 2.8.2)

set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project (opt-pilot)

set (PACKAGE \"OPT-PILOT\")
set (PACKAGE_VERSION \"0.1.4\")

ENDIF (NOT OPTP_USED_FROM_OPAL)

#-----------------------------------------------------------------------------
# config.h
IF(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set (GIT_REVISION 0)
    EXECUTE_PROCESS(
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMAND git describe
            OUTPUT_VARIABLE GIT_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        set (GIT_REVISION \"${GIT_OUT}\")
ELSE(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set (GIT_REVISION 0)
ENDIF(EXISTS "${CMAKE_SOURCE_DIR}/.git")

set (BUILD_DATE 0)
EXECUTE_PROCESS(
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMAND date
        OUTPUT_VARIABLE DATE_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set (BUILD_DATE \"${DATE_OUT}\")


configure_file(config.h.in ${CMAKE_SOURCE_DIR}/config.h)

#-----------------------------------------------------------------------------
# Flags
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    SET (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
    )
endif (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

message (STATUS "Build type is: " ${CMAKE_BUILD_TYPE})

# Select flags.
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -Wno-reorder -Wno-enum-compare")
set (CMAKE_CXX_FLAGS_RELEASE "-O2 -Wno-reorder -Wno-enum-compare")
set (CMAKE_CXX_FLAGS_DEBUG  "-O0 -g -Wno-reorder -Wno-enum-compare")

message (STATUS "Building on ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} on $ENV{HOSTNAME}")


#-----------------------------------------------------------------------------
# dependent packages
find_package (MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

#-----------------------------------------------------------------------------

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

IF (NOT OPTP_USED_FROM_OPAL)
find_package (OPAL REQUIRED)
find_package (IPPL REQUIRED)
find_package (H5hut REQUIRED)
find_package (HDF5 REQUIRED)
find_package (GSL REQUIRED)
find_package (Trilinos REQUIRED HINTS $ENV{TRILINOS_PREFIX} $ENV{TRILINOS_DIR} $ENV{TRILINOS})
ENDIF (NOT OPTP_USED_FROM_OPAL)

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

if (Boost_FILESYSTEM_FOUND)
    add_definitions(-DBOOST_FILESYSTEM)
endif (Boost_FILESYSTEM_FOUND)


### Uldis Locans 2014.09.05 ###
OPTION (ENABLE_DKS "Enable DKS" OFF)

### Load DKS package ###
IF (ENABLE_DKS)
   FIND_PACKAGE(DKS 1.1.1 REQUIRED HINTS $ENV{DKS_PREFIX} $ENV{DKS_DIR} $ENV{DKS})
   MESSAGE (STATUS "DKS Version ${DKS_VERSION}")
   MESSAGE (STATUS "Found DKS library: ${DKS_LIBRARY}")
   MESSAGE (STATUS "Found DKS include dir: ${DKS_INCLUDE_DIR}")
   SET (DKS_PACKAGE_VERSION "\"${DKS_VERSION}\"")
   SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDKS_MPI")

   ### OpenCL compiler flags ###
   SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lOpenCL -pthread -DDKS_OPENCL")

   ### CUDA compiler flags ###
   SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lcudart -lcufft -lcublas -lnvToolsExt -DDKS_CUDA")

   ### if any accelerator enabled set flag to use DKS ###
   SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPPL_DKS -DIPPL_DKS_CUDA -DOPAL_DKS")
ENDIF (ENABLE_DKS)
### END DKS part ###


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
            list (APPEND OPTP_SRCS "${CMAKE_SOURCE_DIR}/${_relPath}/${_src}")
        else()
            list (APPEND OPTP_SRCS "${CMAKE_SOURCE_DIR}/${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SRCS to parent directory
        set (OPTP_SRCS ${OPTP_SRCS} PARENT_SCOPE)
    endif()
endmacro()

set (OPTP_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)

SET (OPTP_USED_FROM_OPAL "TRUE")

add_subdirectory (Expression)
add_subdirectory (Util)
add_subdirectory (extlib)

set (PILOT_LIBS
     pthread
     m
     boost_chrono
     boost_regex
     boost_system
     boost_serialization
#    nlopt
     Problems
)

IF (NOT OPTP_USED_FROM_OPAL)
#-----------------------------------------------------------------------------
# Executables
add_subdirectory (Drivers)

#-----------------------------------------------------------------------------
# Testing
enable_testing()
add_subdirectory(Tests)
ENDIF (NOT OPTP_USED_FROM_OPAL)
#-----------------------------------------------------------------------------
# Doxygen
find_package(Doxygen)
if (DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif(DOXYGEN_FOUND)
