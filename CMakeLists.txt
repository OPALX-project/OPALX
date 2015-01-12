CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (OPAL)
SET (OPAL_VERSION_MAJOR 1)
SET (OPAL_VERSION_MINOR 3.0)


IF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    SET (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
    )
ENDIF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

MESSAGE (STATUS "Build type is: " ${CMAKE_BUILD_TYPE})
MESSAGE (STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")

# Select flags.
SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ")
SET (CMAKE_CXX_FLAGS_RELEASE "-O3")
SET (CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

# Resolve all library dependencies
SET (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

FIND_PACKAGE (MPI REQUIRED)

FIND_PACKAGE (IPPL REQUIRED HINTS $ENV{IPPL_PREFIX} $ENV{IPPL_DIR} $ENV{IPPL})
MESSAGE (STATUS "Found IPPL library: ${IPPL_LIBRARY}")
MESSAGE (STATUS "Found IPPL include dir: ${IPPL_INCLUDE_DIR}")

FIND_PACKAGE (H5Hut REQUIRED)
FIND_PACKAGE (HDF5 REQUIRED)
FIND_PACKAGE (GSL REQUIRED)

set (BOOSTROOT $ENV{BOOST_DIR})
find_package (Boost 1.55.0 REQUIRED)
if (Boost_INCLUDE_DIRS)
    message (STATUS "Found boost include dir: ${Boost_INCLUDE_DIR}")
    include_directories (${Boost_INCLUDE_DIRS})
endif ()

# Handle options
OPTION (ENABLE_AMR_SOLVER "Enable BoxLib based AMR solver" OFF)

OPTION (ENABLE_SAAMG_SOLVER "Enable iteartive SA-AMG-PCG self field solver" OFF)

OPTION (DBG_SCALARFIELD "Enable dump of scalar field rho_m" OFF)
OPTION (DBG_STENCIL "Dump stencil of iterative solver in matlab format" OFF)
OPTION (BUILD_LIBOPAL "Enable building the OPAL library" OFF)

OPTION (NOCPLUSPLUS11_FOREACH "Disable C++11 foreach support" OFF)
OPTION (NOCPLUSPLUS11_NULLPTR "Disable C++11 nullptr support" OFF)

OPTION (NO_FIELD_ASSIGN_OPTIMIZATION "Disable compiler optimization of IPPL field assignment" OFF)

if (ENABLE_AMR_SOLVER)
   message (STATUS "----> Enable AMR")

   find_package (CCSE REQUIRED HINTS $ENV{CCSE_INCLUDE_DIRS} $ENV{CCSE_LIBRARY_DIR})
   include ($ENV{BOXLIB_HOME}/Tools/CMake/CCSEOptions.cmake)

   set (BL_SPACEDIM 3 CACHE INT "Dimension of BoxLib build")
   set (ENABLE_MPI 1 CACHE INT "Enable build with MPI")
   set (ENABLE_OpenMP 0 CACHE INT "Enable build with OpenMP")
   set (BL_PRECISION "DOUBLE" CACHE INT "Precision of BoxLib build")
   set (BL_USE_PARTICLES 1 CACHE INT "Include Particles classes in BoxLib build")
   set (ENABLE_PROFILING 0 CACHE INT "Include profiling information in BoxLib build")
   set (ENABLE_BACKTRACE 1 CACHE INT "Include backtrace information in BoxLib build")
   set (CMAKE_CXX_FLAGS "-Wno-sign-compare -DHAVE_AMR_SOLVER -DTULIN ${CMAKE_CXX_FLAGS}")
   set (CCSE_INCLUDE_DIRS $ENV{CCSE_INCLUDE_DIRS})
   set (CCSE_LIBRARY_DIR $ENV{CCSE_LIBRARY_DIR})

   set(CCSE_LIBRARIES box_camrdata;cboxlib;cfboxlib;fboxlib)
   foreach (L ${CCSE_LIBRARIES})
            find_library(CCSE_LIBRARY
                         NAMES ${L}
                         HINTS ${CCSE_LIBRARY_DIR}
                         NO_DEFAULT_PATH)

            if ( NOT CCSE_LIBRARY )
                message(SEND_ERROR "Cannot locate CCSE library: ${L}")
            endif()

    endforeach()
    message (STATUS "----> Enable AMR done")
endif (ENABLE_AMR_SOLVER)

if (ENABLE_SAAMG_SOLVER OR ENABLE_AMR_SOLVER)
    find_package (Trilinos REQUIRED HINTS $ENV{TRILINOS_PREFIX} $ENV{TRILINOS_DIR} $ENV{TRILINOS})

    message (STATUS "Found Trilinos: ${Trilinos_DIR}")
    message (STATUS "      Trilinos version: ${Trilinos_VERSION}")
    message (STATUS "      Trilinos package list: ${Trilinos_PACKAGE_LIST}")
    message (STATUS "      Trilinos libraries: ${Trilinos_LIBRARIES}")
    message (STATUS "      Trilinos TPL libraries: ${Trilinos_TPL_LIBRARIES}")

    if (NOT ${CMAKE_CXX_COMPILER} STREQUAL ${Trilinos_CXX_COMPILER} )
        message (STATUS "Compiler mismatch:")
        message (STATUS "        Trilinos was compiled with:")
        message (STATUS "            ${Trilinos_C_COMPILER}")
        message (STATUS "            ${Trilinos_CXX_COMPILER}")
        message (STATUS "            ${Trilinos_Fortran_COMPILER}")
        message (STATUS "        You are using:")
        message (STATUS "            ${CMAKE_C_COMPILER}")
        message (STATUS "            ${CMAKE_CXX_COMPILER}")
        message (STATUS "            ${CMAKE_Fortran_COMPILER}")
    endif ()
endif ()

IF (DBG_SCALARFIELD)
   MESSAGE (STATUS "\nWrite scalar rho_m field is enabled ")
   SET (CMAKE_CXX_FLAGS "-DDBG_SCALARFIELD ${CMAKE_CXX_FLAGS}")
ENDIF (DBG_SCALARFIELD)

IF (DBG_STENCIL)
   MESSAGE (STATUS "\nWrite stencil is enabled ")
   SET (CMAKE_CXX_FLAGS "-DDBG_STENCIL ${CMAKE_CXX_FLAGS}")
ENDIF (DBG_STENCIL)

IF (BUILD_LIBOPAL)
   MESSAGE (STATUS "\nBuilding OPAL library")
ENDIF (BUILD_LIBOPAL)

IF (NOCPLUSPLUS11_FOREACH)
    MESSAGE (STATUS "\nBuilding OPAL without c++11 foreach support")
    SET (CMAKE_CXX_FLAGS
        "-DOPAL_NOCPLUSPLUS11_FOREACH ${CMAKE_CXX_FLAGS}"
    )
ENDIF (NOCPLUSPLUS11_FOREACH)

IF (NOCPLUSPLUS11_NULLPTR)
    MESSAGE (STATUS "\nBuilding OPAL without c++11 nullptr support")
    SET (CMAKE_CXX_FLAGS
        "-DOPAL_NOCPLUSPLUS11_NULLPTR ${CMAKE_CXX_FLAGS}"
    )
ENDIF (NOCPLUSPLUS11_NULLPTR)

IF (NO_FIELD_ASSIGN_OPTIMIZATION)
    MESSAGE (STATUS "\nCompiler will not try to optimize field assignment")
    SET (CMAKE_CXX_FLAGS
        "-DdontOPTIMIZE_FIELD_ASSIGNMENT ${CMAKE_CXX_FLAGS}"
    )
ENDIF (NO_FIELD_ASSIGN_OPTIMIZATION)

# compiler dependent flags

SET (MPI_UNDERLYING_COMPILER "")

EXECUTE_PROCESS (COMMAND ${CMAKE_CXX_COMPILER} -show
                OUTPUT_VARIABLE MPI_COMPILER_OUTPUT
                ERROR_VARIABLE MPI_COMPILER_ERROR)
IF ("${MPI_COMPILER_ERROR}" STREQUAL "")
    STRING (REGEX REPLACE "([A-Za-z.0-9/]*) .*$" "\\1" MPI_UNDERLYING_COMPILER ${MPI_COMPILER_OUTPUT})
    STRING (REGEX REPLACE "/" ";" MPI_UNDERLYING_COMPILER_SPLIT ${MPI_UNDERLYING_COMPILER})
    LIST (REVERSE MPI_UNDERLYING_COMPILER_SPLIT)
    LIST (GET MPI_UNDERLYING_COMPILER_SPLIT 0 MPI_UNDERLYING_COMPILER_NAME)
ENDIF ("${MPI_COMPILER_ERROR}" STREQUAL "")

IF ("${MPI_UNDERLYING_COMPILER_NAME}" STREQUAL "icpc")
    # using intel compiler
    MESSAGE (STATUS "The underlying compiler is: icpc")
    SET (CMAKE_CXX_FLAGS
         "-diag-disable 383 -diag-disable 981 $ENV{CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${OTHER_CMAKE_CXX_FLAGS}")

ELSEIF ("${MPI_UNDERLYING_COMPILER_NAME}" STREQUAL "clang++")
    # using clang compiler
    MESSAGE (STATUS "The underlying compiler is: clang++")
    SET (CMAKE_CXX_FLAGS
         "-Wall -Werror -Wsign-compare -Wunused-variable -Warray-bounds -DIPPL_RESTRICT_BUG $ENV{CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${OTHER_CMAKE_CXX_FLAGS}")

ELSEIF ("${MPI_UNDERLYING_COMPILER_NAME}" STREQUAL "g++")
    # using gnu compiler
    MESSAGE (STATUS "The underlying compiler is: g++")
    EXECUTE_PROCESS (COMMAND ${MPI_UNDERLYING_COMPILER} --version
                    OUTPUT_VARIABLE GCC_VERSION_OUTPUT
                    ERROR_VARIABLE GCC_VERSION_ERROR)
    STRING (REGEX REPLACE ".*([0-9]\\.[0-9]\\.[0-9]).*" "\\1" GCC_VERSION ${GCC_VERSION_OUTPUT})
    MESSAGE (STATUS "Your gcc version is: ${GCC_VERSION}")

    IF (${GCC_VERSION} VERSION_LESS "4.5.0")
        MESSAGE (FATAL_ERROR "To build OPAL you need gcc version 4.5.0 or greater")

    ELSEIF (${GCC_VERSION} VERSION_LESS "4.6.0")
        SET (CMAKE_CXX_FLAGS
             "-DOPAL_NOCPLUSPLUS11_NULLPTR -DOPAL_NOCPLUSPLUS11_FOREACH $ENV{CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${OTHER_CMAKE_CXX_FLAGS}")

    ELSEIF (${GCC_VERSION} VERSION_LESS "4.7.0")
        SET (CMAKE_CXX_FLAGS
    	"-Wall -Werror -Wno-reorder -std=c++0x $ENV{CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${OTHER_CMAKE_CXX_FLAGS}")

    ELSEIF (${GCC_VERSION} VERSION_EQUAL "4.7.0")
        SET (CMAKE_CXX_FLAGS
            "-DdontOPTIMIZE_FIELD_ASSIGNMENT $ENV{CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${OTHER_CMAKE_CXX_FLAGS}")

    ELSEIF (${GCC_VERSION} VERSION_GREATER "4.7")
        SET (CMAKE_CXX_FLAGS
        "-Wall -Werror -Wno-reorder -Wno-unused-local-typedefs -Werror=unused-variable -std=c++11 $ENV{CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${OTHER_CMAKE_CXX_FLAGS}")
    ENDIF ()

ELSE ()
    # using unknown compiler
    MESSAGE (STATUS "${CMAKE_CXX_COMPILER} is using a compiler (${MPI_UNDERLYING_COMPILER}) we were not thinking of!
             Please use the gnu compiler or the intel compiler IF you are having problems.")
ENDIF ()


MESSAGE (STATUS "Compiling with ${CMAKE_CXX_FLAGS}")


INCLUDE(FindSubversion)
IF(Subversion_SVN_EXECUTABLE)

    EXECUTE_PROCESS(COMMAND ${Subversion_SVN_EXECUTABLE} info ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE OPAL_WC_INFO
        ERROR_VARIABLE Subversion_svn_info_error
        RESULT_VARIABLE Subversion_svn_info_result
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    IF(${Subversion_svn_info_result} EQUAL 0)
        Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} OPAL)
    ELSE()
        SET (OPAL_WC_REVISION 0)
        EXECUTE_PROCESS(
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            COMMAND git rev-parse HEAD
            OUTPUT_VARIABLE GIT_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
            )
        SET (OPAL_WC_REVISION \"${GIT_OUT}\")
    ENDIF()
ELSE()
    MESSAGE (STATUS "could not find subversion")
    SET (OPAL_WC_REVISION \"\")
ENDIF()

ADD_SUBDIRECTORY (classic/5.0/src)
ADD_SUBDIRECTORY (src)
