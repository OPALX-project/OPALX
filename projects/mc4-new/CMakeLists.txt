CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (MC4)
SET (MC4_VERSION_MAJOR 1)
SET (MC4_VERSION_MINOR 1.9)

SET (CMAKE_CXX_COMPILER ${MPI_COMPILER})

IF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    SET (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
    )
ENDIF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

MESSAGE (STATUS "Build type is: " ${CMAKE_BUILD_TYPE})

# Select flags.
SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-Wno-deprecated -O2 -g")
SET (CMAKE_CXX_FLAGS_RELEASE "-Wno-deprecated -O2")
SET (CMAKE_CXX_FLAGS_DEBUG  "-Wno-deprecated -O0 -g")


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
OPTION (BUILD_LIBMC4 "Enable building the MC4 library" OFF)

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

IF (BUILD_LIBMC4)
   MESSAGE (STATUS "\nBuilding MC4 library")
ENDIF (BUILD_LIBMC4)

IF (NOCPLUSPLUS11_FOREACH)
    MESSAGE (STATUS "\nBuilding MC4 without c++11 foreach support")
    SET (CMAKE_CXX_FLAGS
        "-DMC4_NOCPLUSPLUS11_FOREACH ${CMAKE_CXX_FLAGS}"
    )
ENDIF (NOCPLUSPLUS11_FOREACH)

IF (NOCPLUSPLUS11_NULLPTR)
    MESSAGE (STATUS "\nBuilding MC4 without c++11 nullptr support")
    SET (CMAKE_CXX_FLAGS
        "-DMC4_NOCPLUSPLUS11_NULLPTR ${CMAKE_CXX_FLAGS}"
    )
ENDIF (NOCPLUSPLUS11_NULLPTR)

IF (NO_FIELD_ASSIGN_OPTIMIZATION)
    MESSAGE (STATUS "\nCompiler will not try to optimize field assignment")
    SET (CMAKE_CXX_FLAGS
        "-DdontOPTIMIZE_FIELD_ASSIGNMENT ${CMAKE_CXX_FLAGS}"
    )
ENDIF (NO_FIELD_ASSIGN_OPTIMIZATION)


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
            MESSAGE (FATAL_ERROR "To build MC4 you need gcc version 4.5.0 or greater")
        ELSEIF (${GCC_VERSION} VERSION_LESS "4.6.0")
            SET (OTHER_CMAKE_CXX_FLAGS
                 "-DMC4_NOCPLUSPLUS11_NULLPTR -DMC4_NOCPLUSPLUS11_FOREACH ${OTHER_CMAKE_CXX_FLAGS}")
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
    "-Wall -Werror -Wparentheses -std=c++0x -DUSENAMESPACE ${CMAKE_CXX_FLAGS} ${OTHER_CMAKE_CXX_FLAGS}")

#SET (CMAKE_CXX_FLAGS
#  "-Wall -Werror -std=c++0x ${CMAKE_CXX_FLAGS} ${OTHER_CMAKE_CXX_FLAGS}")

MESSAGE (STATUS "Compiling with ${CMAKE_CXX_FLAGS}")


ADD_SUBDIRECTORY (src)

