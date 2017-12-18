CMAKE_MINIMUM_REQUIRED (VERSION 2.8.10)
PROJECT (OPAL)
SET (OPAL_VERSION_MAJOR 1)
SET (OPAL_VERSION_MINOR 9.0)
set (PACKAGE \"opal\")
set (PACKAGE_BUGREPORT \"opal@lists.psi.ch\")
set (PACKAGE_NAME \"OPAL\")
set (PACKAGE_TARNAME \"opal\")
set (PACKAGE_VERSION "${OPAL_VERSION_MAJOR}.${OPAL_VERSION_MINOR}")
set (PACKAGE_VERSION_STR "\"${PACKAGE_VERSION}\"")
STRING (REGEX REPLACE "\\.([0-9])\\." ".0\\1." PACKAGE_VERSION ${PACKAGE_VERSION})
STRING (REGEX REPLACE "\\.([0-9])$" ".0\\1" PACKAGE_VERSION ${PACKAGE_VERSION})
STRING (REGEX REPLACE "\\." "" PACKAGE_VERSION ${PACKAGE_VERSION})

IF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    SET (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
	)
ENDIF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)


CMAKE_HOST_SYSTEM_INFORMATION(RESULT HOSTNAME QUERY HOSTNAME)

MESSAGE (STATUS "Build type is: " ${CMAKE_BUILD_TYPE})
MESSAGE (STATUS "Host OS System: ${CMAKE_HOST_SYSTEM}")
MESSAGE (STATUS "Hostname:  ${HOSTNAME}")

# Select flags.
SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ")
SET (CMAKE_CXX_FLAGS_RELEASE "-O3")
SET (CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

# Resolve all library dependencies
SET (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

find_package (MPI REQUIRED)
message (STATUS "The C++ compiler identification is: ${CMAKE_CXX_COMPILER_ID}")
message (STATUS "The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
message (STATUS "The MPI C++ compiler is: ${MPI_CXX_COMPILER}")
message (STATUS "The underlying C++ compiler is: ${CMAKE_CXX_COMPILER}")

# compiler dependent flags
if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -diag-disable 383 -diag-disable 981 -std=c++11")

elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wsign-compare -Wunused-variable -Warray-bounds -DIPPL_RESTRICT_BUG -std=c++11")

elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    if (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.7.3")
        message (FATAL_ERROR "To build OPAL you need GCC version 4.7.3 or greater")
    endif ()
    # Warning flags, -Wno-unused-local-typedefs is needed for boost versions below 1.59
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-local-typedefs -Werror=unused-variable -std=c++11")
else ()
    message (STATUS "Unknown C++ compiler. Please use the GNU or Intel compiler, IF you are having problems.")
endif ()

FIND_PACKAGE (H5Hut REQUIRED)
FIND_PACKAGE (HDF5 REQUIRED)
FIND_PACKAGE (GSL REQUIRED)

set (BOOSTROOT $ENV{BOOST_DIR})
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package (Boost 1.55.0 REQUIRED COMPONENTS filesystem system iostreams regex)
if (Boost_INCLUDE_DIRS)
    message (STATUS "Found boost include dir: ${Boost_INCLUDE_DIR}")
    message (STATUS "Found boost library dir: ${Boost_LIBRARY_DIR}")
    message (STATUS "Found boost libraries: ${Boost_LIBRARIES}")
    include_directories (${Boost_INCLUDE_DIRS})
endif ()


### Uldis Locans 2014.09.05 ###
OPTION (ENABLE_DKS "Enable DKS" OFF)

### Load DKS package ###
IF (ENABLE_DKS)
    FIND_PACKAGE(DKS 1.1.1 REQUIRED HINTS $ENV{DKS_PREFIX} $ENV{DKS_DIR} $ENV{DKS})
    set (DKS_VERSION_OPAL \"${DKS_VERSION}\")
    MESSAGE (STATUS "Found DKS version: ${DKS_VERSION}")
    MESSAGE (STATUS "Found DKS library: ${DKS_LIBRARY}")
    MESSAGE (STATUS "Found DKS include dir: ${DKS_INCLUDE_DIR}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDKS_MPI")

    ### OpenCL compiler flags ###
    #SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lOpenCL -pthread -DDKS_OPENCL")

    ### CUDA compiler flags ###
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDKS_CUDA")

    ### if any accelerator enabled set flag to use DKS ###
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPPL_DKS -DIPPL_DKS_CUDA -DOPAL_DKS")
ENDIF (ENABLE_DKS)
### END DKS part ###

# Handle options
OPTION (ENABLE_AMR "Enable AMReX based AMR solver" OFF)
OPTION (ENABLE_AMR_MG_SOLVER "Enable AMR MG solver" OFF)

OPTION (ENABLE_SAAMG_SOLVER "Enable iteartive SA-AMG-PCG self field solver" OFF)

OPTION (DBG_SCALARFIELD "Enable dump of scalar field rho_m" OFF)
OPTION (DBG_STENCIL "Dump stencil of iterative solver in matlab format" OFF)
OPTION (AMR_YT_DUMP "Dump AMR field data to yt readable format" OFF)
OPTION (AMR_PYTHON_DUMP "Dump AMR grids and particles" OFF)
OPTION (BUILD_LIBOPAL "Enable building the OPAL library" OFF)

OPTION (NO_FIELD_ASSIGN_OPTIMIZATION "Disable compiler optimization of IPPL field assignment" OFF)

OPTION (BUILD_OPAL_UNIT_TESTS "Unit tests" OFF)

IF (BUILD_OPAL_UNIT_TESTS)
    FIND_PACKAGE (GTest REQUIRED)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_UNIT_TESTS")
ENDIF (BUILD_OPAL_UNIT_TESTS)

IF (USE_H5HUT2)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_H5HUT2")
ENDIF (USE_H5HUT2)

IF (ENABLE_AMR)
    MESSAGE (STATUS "Enable AMR_SOLVER    " ${ENABLE_AMR})
    ENABLE_LANGUAGE (Fortran)

    find_package (AMReX CONFIG REQUIRED
        HINTS $ENV{AMREX_HOME}/cmake $ENV{AMREX_DIR}/cmake $ENV{AMREX_PREFIX}/cmake)

    MESSAGE (STATUS "Found AMReX include dir: ${AMREX_INCLUDE_DIR}")
    MESSAGE (STATUS "Found AMReX library dir: ${AMREX_LIBRARY_DIR}")
    MESSAGE (STATUS "Found AMReX libraries: ${AMREX_LIBRARIES}")
    echo_amrex_config_options()
    
    add_definitions(${AMREX_DEFINES})
    
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_AMR")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable -Wno-unused-but-set-variable -Wno-maybe-uninitialized")

ENDIF (ENABLE_AMR)

# remove digits from hostname: edison03 -> edison
STRING(REGEX MATCH "[^0-9]*" HOSTNAME_BASE "${HOSTNAME}")

# In case of Edison and Cori the find_package(Trilinos) gives the wrong
# Trilinos_TPL_LIBRARIES, instead of cray-tpsl/16.07.1 it uses
# cray-tpsl/16.06.1 libraries -> Linker error.
# We can avoid this issue by not going into the if-statement
IF (ENABLE_SAAMG_SOLVER)
    MESSAGE (STATUS "Enable SAAMG_SOLVER  " ${ENABLE_SAAMG_SOLVER})
    IF (NOT (
                ${HOSTNAME_BASE} MATCHES "edison" OR
                ${HOSTNAME_BASE} MATCHES "cori"))
        
        find_package (Trilinos REQUIRED HINTS $ENV{TRILINOS_PREFIX} $ENV{TRILINOS_DIR} $ENV{TRILINOS})

        message (STATUS "Found Trilinos: ${Trilinos_DIR}")
        message (STATUS "      Trilinos version: ${Trilinos_VERSION}")
        message (STATUS "      Trilinos package list: ${Trilinos_PACKAGE_LIST}")
        message (STATUS "      Trilinos libraries: ${Trilinos_LIBRARIES}")
        message (STATUS "      Trilinos TPL libraries: ${Trilinos_TPL_LIBRARIES}")

        if (NOT ${MPI_CXX_COMPILER} STREQUAL ${Trilinos_CXX_COMPILER} )
            message (STATUS "Compiler mismatch:")
            message (STATUS "        Trilinos was compiled with:")
            message (STATUS "            ${Trilinos_C_COMPILER}")
            message (STATUS "            ${Trilinos_CXX_COMPILER}")
            message (STATUS "            ${Trilinos_Fortran_COMPILER}")
            message (STATUS "        You are using:")
            message (STATUS "            ${CMAKE_C_COMPILER}")
            message (STATUS "            ${CMAKE_CXX_COMPILER}")
            message (STATUS "            ${CMAKE_Fortran_COMPILER}")
            message (FATAL  "")
        ENDIF ()
    ENDIF ()
ENDIF ()

IF (DBG_SCALARFIELD)
    MESSAGE (STATUS "\nWrite scalar rho_m field is enabled ")
    SET (CMAKE_CXX_FLAGS "-DDBG_SCALARFIELD ${CMAKE_CXX_FLAGS}")
ENDIF (DBG_SCALARFIELD)

IF (DBG_STENCIL)
    MESSAGE (STATUS "\nWrite stencil is enabled ")
    SET (CMAKE_CXX_FLAGS "-DDBG_STENCIL ${CMAKE_CXX_FLAGS}")
ENDIF (DBG_STENCIL)

IF (AMR_YT_DUMP AND ENABLE_AMR)
    MESSAGE (STATUS "\nWrite AMR field data in yt format is enabled ")
    SET (CMAKE_CXX_FLAGS "-DAMR_YT_DUMP ${CMAKE_CXX_FLAGS}")
ENDIF (AMR_YT_DUMP AND ENABLE_AMR)

IF (AMR_PYTHON_DUMP AND ENABLE_AMR)
    MESSAGE (STATUS "\nWrite AMR grids and particles is enabled ")
    SET (CMAKE_CXX_FLAGS "-DAMR_PYTHON_DUMP ${CMAKE_CXX_FLAGS}")
ENDIF (AMR_PYTHON_DUMP AND ENABLE_AMR)

IF (BUILD_LIBOPAL)
    MESSAGE (STATUS "\nBuilding OPAL library")
ENDIF (BUILD_LIBOPAL)

IF (NO_FIELD_ASSIGN_OPTIMIZATION)
    MESSAGE (STATUS "\nCompiler will not try to optimize field assignment")
    SET (CMAKE_CXX_FLAGS
        "-DdontOPTIMIZE_FIELD_ASSIGNMENT ${CMAKE_CXX_FLAGS}"
	)
ENDIF (NO_FIELD_ASSIGN_OPTIMIZATION)

INCLUDE(FindGit)
IF(GIT_FOUND)
    SET (OPAL_WC_REVISION 0)
    EXECUTE_PROCESS(
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMAND git rev-parse HEAD
        OUTPUT_VARIABLE GIT_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    SET (OPAL_WC_REVISION \"${GIT_OUT}\")
ELSE()
    MESSAGE (STATUS "could not find git")
    SET (OPAL_WC_REVISION \"\")
ENDIF()

SET (IPPL_USED_FROM_OPAL "TRUE")

ADD_SUBDIRECTORY (ippl)
ADD_SUBDIRECTORY (opt-pilot)

ADD_SUBDIRECTORY (src)
ADD_SUBDIRECTORY (tools)

IF (BUILD_OPAL_UNIT_TESTS)
    ADD_SUBDIRECTORY(tests)
ENDIF (BUILD_OPAL_UNIT_TESTS)

# make variables visible in other CMakeLists files
set(OPAL_CXX_FLAGS ${OPAL_CXX_FLAGS}
    CACHE INTERNAL "" FORCE
    )

set(OPAL_LIBS ${OPAL_LIBS}
    CACHE INTERNAL "" FORCE
    )

# configure file
configure_file (
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake
    )

install (
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake
    DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/${PROJECT_NAME}"
    RENAME ${PROJECT_NAME}Config.cmake
    )

# vi: set et ts=4 sw=4 sts=4:
# Local Variables:
# mode:cmake
# cmake-tab-width: 4
# indent-tabs-mode:nil
# End:
