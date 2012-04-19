CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (OPAL)
SET (OPAL_VERSION_MAJOR 1)
SET (OPAL_VERSION_MINOR 1.9)

SET (CMAKE_CXX_COMPILER ${MPI_COMPILER})

IF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    SET (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
    )
ENDIF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

MESSAGE (STATUS "Build type is: " ${CMAKE_BUILD_TYPE})

# Select flags.
SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
SET (CMAKE_CXX_FLAGS_RELEASE "-O2")
SET (CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")


# Resolve all library dependencies
SET (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

FIND_PACKAGE (IPPL REQUIRED HINTS $ENV{IPPL_PREFIX})
MESSAGE (STATUS "Found IPPL library dir: ${IPPL_LIBRARY_DIR}")
MESSAGE (STATUS "Found IPPL include dir: ${IPPL_INCLUDE_DIR}")

FIND_PACKAGE (H5Hut REQUIRED)
FIND_PACKAGE (HDF5 REQUIRED)
FIND_PACKAGE (GSL REQUIRED)


# Handle options
OPTION (ENABLE_ML_SOLVER "Enable iteartive SA-AMG-PCG self field solver" OFF)

OPTION (DBG_SCALARFIELD "Enable dump of scalar field rho_m" OFF)
OPTION (DBG_STENCIL "Dump stencil of iterative solver in matlab format" OFF)
OPTION (BUILD_LIBOPAL "Enable building the OPAL library" OFF)

OPTION (NOCPLUSPLUS11_FOREACH "Disable C++11 foreach support" OFF)
OPTION (NOCPLUSPLUS11_NULLPTR "Disable C++11 nullptr support" OFF)

OPTION (NO_FIELD_ASSIGN_OPTIMIZATION "Disable compiler optimization of IPPL field assignment" OFF)

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

IF (ENABLE_ML_SOLVER)
    # Get Trilinos as one entity
    FIND_PACKAGE (Trilinos PATHS $ENV{TRILINOS_INCLUDE_PATH})

    IF (NOT Trilinos_FOUND)
        MESSAGE (FATAL_ERROR "Could not find Trilinos!")
    ENDIF (NOT Trilinos_FOUND)
    # Echo trilinos build info just for fun
    MESSAGE (STATUS "\nFound Trilinos!  Here are the details: ")
    MESSAGE (STATUS "   Trilinos_DIR = ${Trilinos_DIR}")
    MESSAGE (STATUS "   Trilinos_VERSION = ${Trilinos_VERSION}")
    MESSAGE (STATUS "   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
    MESSAGE (STATUS "   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
    MESSAGE (STATUS "   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
    MESSAGE (STATUS "   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
    MESSAGE (STATUS "   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
    MESSAGE (STATUS "   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
    MESSAGE (STATUS "   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
    MESSAGE (STATUS "   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
    MESSAGE (STATUS "   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
    MESSAGE (STATUS "End of Trilinos details\n")

    # Make sure to use same compilers and flags as Trilinos
    #SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
    #SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
    #SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

    #SET(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
    #SET(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
    #SET(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")
ENDIF (ENABLE_ML_SOLVER)


# compiler dependent flags
STRING (REGEX REPLACE ".*/([A-Za-z]*)$" "\\1" COMPILER_NAME ${CMAKE_CXX_COMPILER})
MESSAGE (STATUS "Your compiler is: ${COMPILER_NAME}")

IF ("${COMPILER_NAME}" STREQUAL "mpicxx")

    SET (MPI_UNDERLYING_COMPILER
        "")

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
        MESSAGE (STATUS "The underlying compiler of ${COMPILER_NAME} is: icpc")
        SET (OTHER_CMAKE_CXX_FLAGS
             "-diag-disable 383 -diag-disable 981 ${OTHER_CMAKE_CXX_FLAGS}")

    ELSEIF ("${MPI_UNDERLYING_COMPILER_NAME}" STREQUAL "clang++")
        # using clang compiler
        MESSAGE (STATUS "The underlying compiler of ${COMPILER_NAME} is: clang++")
        SET (OTHER_CMAKE_CXX_FLAGS
             "-Wsign-compare -Wunused-variable -Warray-bounds -DIPPL_RESTRICT_BUG ${OTHER_CMAKE_CXX_FLAGS}")

    ELSEIF ("${MPI_UNDERLYING_COMPILER_NAME}" STREQUAL "g++")
        # using gnu compiler
        MESSAGE (STATUS "The underlying compiler of ${COMPILER_NAME} is: g++")
        EXECUTE_PROCESS (COMMAND ${MPI_UNDERLYING_COMPILER} --version
                        OUTPUT_VARIABLE GCC_VERSION_OUTPUT
                        ERROR_VARIABLE GCC_VERSION_ERROR)
        STRING (REGEX REPLACE ".*([0-9]\\.[0-9]\\.[0-9]).*" "\\1" GCC_VERSION ${GCC_VERSION_OUTPUT})
        MESSAGE (STATUS "Your gcc version is: ${GCC_VERSION}")

        IF (${GCC_VERSION} VERSION_LESS "4.5.0")
            MESSAGE (FATAL_ERROR "To build OPAL you need gcc version 4.5.0 or greater")
        ELSEIF (${GCC_VERSION} VERSION_LESS "4.6.0")
            SET (OTHER_CMAKE_CXX_FLAGS
                 "-DOPAL_NOCPLUSPLUS11_NULLPTR -DOPAL_NOCPLUSPLUS11_FOREACH ${OTHER_CMAKE_CXX_FLAGS}")
        ELSEIF (${GCC_VERSION} VERSION_EQUAL "4.7.0")
            SET (OTHER_CMAKE_CXX_FLAGS
                "-DdontOPTIMIZE_FIELD_ASSIGNMENT ${OTHER_CMAKE_CXX_FLAGS}")
        ENDIF (${GCC_VERSION} VERSION_LESS "4.5.0")

    ELSE ("${MPI_UNDERLYING_COMPILER_NAME}" STREQUAL "icpc")
        # using unknown compiler
        MESSAGE (STATUS "${COMPILER_NAME} is using a compiler (${MPI_UNDERLYING_COMPILER}) we were not thinking of!
                 Please use the gnu compiler or the intel compiler IF you are having problems.")
    ENDIF ("${MPI_UNDERLYING_COMPILER_NAME}" STREQUAL "icpc")

ELSE ("${COMPILER_NAME}" STREQUAL "mpicxx")
    # using unknown mpi implementation
    MESSAGE (STATUS "You are using an unsupported MPI compiler: ${COMPILER_NAME}!
             Please use the OpenMPI or MPICH IF you are having problems.")
ENDIF ("${COMPILER_NAME}" STREQUAL "mpicxx")

SET (CMAKE_CXX_FLAGS
  "-Wall -Werror -std=c++0x ${CMAKE_CXX_FLAGS} ${OTHER_CMAKE_CXX_FLAGS}")

MESSAGE (STATUS "Compiling with ${CMAKE_CXX_FLAGS}")


ADD_SUBDIRECTORY (classic/5.0/src)
ADD_SUBDIRECTORY (src)

