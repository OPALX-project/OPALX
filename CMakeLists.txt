CMAKE_MINIMUM_REQUIRED (VERSION 2.8.10)
PROJECT (OPAL)
SET (OPAL_VERSION_MAJOR 1)
SET (OPAL_VERSION_MINOR 3.99)

IF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    SET (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
	)
ENDIF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

MESSAGE (STATUS "Build type is: " ${CMAKE_BUILD_TYPE})

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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -diag-disable 383 -diag-disable 981")

elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wsign-compare -Wunused-variable -Warray-bounds -DIPPL_RESTRICT_BUG")

elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    if (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.7.3")
	message (FATAL_ERROR "To build OPAL you need GCC version 4.7.3 or greater")

    elseif (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "4.7.2")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wno-reorder -Wno-unused-local-typedefs -Werror=unused-variable -std=c++11")
    endif ()
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
find_package (Boost 1.55.0 REQUIRED COMPONENTS  regex filesystem system)
if (Boost_INCLUDE_DIRS)
    message (STATUS "Found boost include dir: ${Boost_INCLUDE_DIR}")
    message (STATUS "Found boost library dir: ${Boost_LIBRARY_DIR}")
    message (STATUS "Found boost libraries: ${Boost_LIBRARIES}")
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

ADD_SUBDIRECTORY (ippl)
ADD_SUBDIRECTORY (classic/5.0/src)
ADD_SUBDIRECTORY (src)

IF (BUILD_OPAL_UNIT_TESTS)
    FIND_PACKAGE (GTest REQUIRED)

    # touch CMakeLists to force cmake to redo the glob for new test source files
    # otherwise cmake does not detect new source files
    EXECUTE_PROCESS(COMMAND touch ${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt)
    ADD_SUBDIRECTORY(tests)

ENDIF (BUILD_OPAL_UNIT_TESTS)