CMAKE_MINIMUM_REQUIRED (VERSION 3.1)
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

MESSAGE ("Build type is: " ${CMAKE_BUILD_TYPE})
MESSAGE ("Host OS System: ${CMAKE_HOST_SYSTEM}")
MESSAGE ("Hostname:  ${HOSTNAME}")

# Select flags.
SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ")
SET (CMAKE_CXX_FLAGS_RELEASE "-O3")
SET (CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

# Resolve all library dependencies
SET (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

find_package (MPI REQUIRED)
message ("The C++ compiler identification is: ${CMAKE_CXX_COMPILER_ID}")
message ("The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
message ("The MPI C++ compiler is: ${MPI_CXX_COMPILER}")
message ("The underlying C++ compiler is: ${CMAKE_CXX_COMPILER}")

OPTION (ENABLE_OpenMP "Use hybrid parallelism MPI-OpenMP" OFF)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

add_compile_options (-Wall)

# compiler dependent flags
if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel" )
    add_compile_options (-diag-disable 383 -diag-disable 981)
    if (ENABLE_OpenMP)
	add_compile_options (-openmp)
    endif ()

elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    add_compile_options (-Werror -Wsign-compare -Wunused-variable -Warray-bounds)
    add_definitions(-DIPPL_RESTRICT_BUG)
    if (ENABLE_OpenMP)
	    add_compile_options (-fopenmp)
    endif ()

elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    if (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.7.3")
        message (FATAL_ERROR "To build OPAL you need GCC version 4.7.3 or greater")
    endif ()
    # Warning flags, -Wno-unused-local-typedefs is needed for boost versions below 1.59
    add_compile_options (-Wno-unused-local-typedefs -Werror=unused-variable)

    if (ENABLE_OpenMP)
        add_compile_options (-fopenmp)
    endif ()

else ()
    message ("Unknown C++ compiler. Please use the GNU or Intel compiler, IF you are having problems.")
endif ()

enable_language (Fortran OPTIONAL)
message ("Fortran compilers: ${CMAKE_Fortran_COMPILER_NAMES}")
if (CMAKE_Fortran_COMPILER)
    message ("The Fortran compiler identification is: ${CMAKE_Fortran_COMPILER_ID}")
    message ("The Fortran compiler version is: ${CMAKE_Fortran_COMPILER_VERSION}")
    message ("The Fortran compiler is: ${CMAKE_Fortran_COMPILER}")
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        set (Fortran_LIBS "-lgfortran")
    endif ()
else ()
    message ("No Fortran compiler found (optional)")
    if (ENABLE_AMR)
	message (FATAL_ERROR "AMR is enabled but no Fortran compiler found!")
    endif ()
endif()

find_package (H5Hut REQUIRED)
set (HDF5_USE_STATIC_LIBRARIES ON)
find_package (HDF5 REQUIRED)
message ("Found HDF5 libraries: ${HDF5_C_LIBRARIES}")
find_package (GSL REQUIRED)
message ("Found gsl libraries: ${GSL_LIBRARIES}")

set (BOOSTROOT $ENV{BOOST_DIR})
set (Boost_USE_STATIC_LIBS ON)
set (Boost_USE_MULTITHREADED OFF)
set (Boost_USE_STATIC_RUNTIME OFF)
find_package (Boost REQUIRED COMPONENTS filesystem system iostreams regex)
if (Boost_INCLUDE_DIRS)
    if (${Boost_VERSION} VERSION_LESS 106000)
	message (FATAL_ERROR
	       "Found Boost library is to old; required is version 1.60.0 or newer!")
    endif ()
    message ("Found Boost include dir: ${Boost_INCLUDE_DIR}")
    message ("Found Boost library dir: ${Boost_LIBRARY_DIR}")
    message ("Found Boost libraries: ${Boost_LIBRARIES}")
    include_directories (${Boost_INCLUDE_DIRS})
endif ()

OPTION (ENABLE_DKS "Enable DKS" OFF)
IF (ENABLE_DKS)
    find_package (dKS 1.1.1 REQUIRED HINTS $ENV{DKS_PREFIX} $ENV{DKS_DIR} $ENV{DKS})
    set (DKS_VERSION_OPAL \"${DKS_VERSION}\")
    MESSAGE ("Found DKS version: ${DKS_VERSION}")
    MESSAGE ("Found DKS library: ${DKS_LIBRARY}")
    MESSAGE ("Found DKS include dir: ${DKS_INCLUDE_DIR}")
    add_definitions(-DDKS_MPI)

    ### OpenCL compiler flags ###
    #SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lOpenCL -pthread -DDKS_OPENCL")

    ### CUDA compiler flags ###
    add_definitions(-DDKS_CUDA)

    ### if any accelerator enabled set flag to use DKS ###
    add_definitions(-DIPPL_DKS -DIPPL_DKS_CUDA -DOPAL_DKS)

ENDIF (ENABLE_DKS)
### END DKS part ###

# Handle options
OPTION (ENABLE_AMR "Enable AMReX based AMR solver" OFF)
OPTION (ENABLE_AMR_MG_SOLVER "Enable AMR MG solver" OFF)

OPTION (ENABLE_SAAMG_SOLVER "Enable iteartive SA-AMG-PCG self field solver" OFF)

OPTION (DBG_SCALARFIELD "Enable dump of scalar field rho_m" OFF)
OPTION (DBG_STENCIL "Dump stencil of iterative solver in matlab format" OFF)

OPTION (NO_FIELD_ASSIGN_OPTIMIZATION "Disable compiler optimization of IPPL field assignment" OFF)

OPTION (BUILD_OPAL_UNIT_TESTS "Unit tests" OFF)

IF (BUILD_OPAL_UNIT_TESTS)
    FIND_PACKAGE (GTest REQUIRED)
    ADD_DEFINITIONS(-DWITH_UNIT_TESTS)
ENDIF (BUILD_OPAL_UNIT_TESTS)

IF (ENABLE_AMR)
    MESSAGE ("Enable AMR_SOLVER    " ${ENABLE_AMR})
    ENABLE_LANGUAGE (Fortran)

    find_package (AMReX CONFIG REQUIRED
        HINTS $ENV{AMREX_HOME}/cmake $ENV{AMREX_DIR}/cmake $ENV{AMREX_PREFIX}/cmake)

    MESSAGE ("Found AMReX include dir: ${AMREX_INCLUDE_DIR}")
    MESSAGE ("Found AMReX library dir: ${AMREX_LIBRARY_DIR}")
    MESSAGE ("Found AMReX libraries: ${AMREX_LIBRARIES}")
    MESSAGE ("Found AMReX Fortran links: ${AMREX_EXTRA_Fortran_LINK_LINE}")
    echo_amrex_config_options()
    
    add_definitions(${AMREX_DEFINES})
    add_definitions(-DENABLE_AMR)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable -Wno-unused-but-set-variable -Wno-maybe-uninitialized")

ENDIF (ENABLE_AMR)

# remove digits from hostname: edison03 -> edison
STRING(REGEX MATCH "[^0-9]*" HOSTNAME_BASE "${HOSTNAME}")

# In case of Edison and Cori the find_package(Trilinos) gives the wrong
# Trilinos_TPL_LIBRARIES, instead of cray-tpsl/16.07.1 it uses
# cray-tpsl/16.06.1 libraries -> Linker error.
# We can avoid this issue by not going into the if-statement
IF (ENABLE_SAAMG_SOLVER OR ENABLE_AMR_MG_SOLVER)

    IF ( ENABLE_SAAMG_SOLVER )
        MESSAGE ("Enable SAAMG_SOLVER  " ${ENABLE_SAAMG_SOLVER})
    ELSE ( ENABLE_SAAMG_SOLVER )
        MESSAGE ("Enable AMR_MG_SOLVER " ${ENABLE_AMR_MG_SOLVER})
    ENDIF( ENABLE_SAAMG_SOLVER )
    
    IF (NOT (
                ${HOSTNAME_BASE} MATCHES "edison" OR
                ${HOSTNAME_BASE} MATCHES "cori"))
        
        find_package (Trilinos REQUIRED HINTS $ENV{TRILINOS_PREFIX} $ENV{TRILINOS_DIR} $ENV{TRILINOS})

        message ("Found Trilinos: ${Trilinos_DIR}")
        message ("      Trilinos version: ${Trilinos_VERSION}")
        message ("      Trilinos package list: ${Trilinos_PACKAGE_LIST}")
        message ("      Trilinos libraries: ${Trilinos_LIBRARIES}")
        message ("      Trilinos TPL libraries: ${Trilinos_TPL_LIBRARIES}")

        if (NOT ${MPI_CXX_COMPILER} STREQUAL ${Trilinos_CXX_COMPILER} )
            message ("Compiler mismatch:")
            message ("        Trilinos was compiled with:")
            message ("            ${Trilinos_C_COMPILER}")
            message ("            ${Trilinos_CXX_COMPILER}")
            message ("            ${Trilinos_Fortran_COMPILER}")
            message ("        You are using:")
            message ("            ${CMAKE_C_COMPILER}")
            message ("            ${CMAKE_CXX_COMPILER}")
            message ("            ${CMAKE_Fortran_COMPILER}")
	    message (FATAL_ERROR  "")
        ENDIF ()
    ENDIF ()
ENDIF ()

IF (DBG_SCALARFIELD)
    MESSAGE ("\nWrite scalar rho_m field is enabled ")
    ADD_DEFINITIONS(-DDBG_SCALARFIELD)
ENDIF (DBG_SCALARFIELD)

IF (DBG_STENCIL)
    MESSAGE ("\nWrite stencil is enabled ")
    ADD_DEFINITIONS(-DDBG_STENCIL)
ENDIF (DBG_STENCIL)

IF (NO_FIELD_ASSIGN_OPTIMIZATION)
    MESSAGE ("\nCompiler will not try to optimize field assignment")
    ADD_DEFINITIONS(-DdontOPTIMIZE_FIELD_ASSIGNMENT)
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
    MESSAGE ("could not find git")
    SET (OPAL_WC_REVISION \"\")
ENDIF()

SET (IPPL_USED_FROM_OPAL "TRUE")

ADD_SUBDIRECTORY (ippl)
ADD_SUBDIRECTORY (optimizer)
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
