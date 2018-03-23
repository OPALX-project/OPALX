cmake_minimum_required (VERSION 3.0.2)
project (OPAL)
set (OPAL_MAJOR_VERSION 1)
set (OPAL_MINOR_VERSION 6)
set (OPAL_PATCH_VERSION 1)
set (OPAL_PACKAGE_BUGREPORT opal@lists.psi.ch)
set (OPAL_PACKAGE_VERSION
  ${OPAL_MAJOR_VERSION}.${OPAL_MINOR_VERSION}.${OPAL_PATCH_VERSION})

# Select flags.
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")
set (CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

cmake_host_system_information (RESULT HOSTNAME QUERY HOSTNAME)

message ("Build type is: " ${CMAKE_BUILD_TYPE})
message ("Host OS System: ${CMAKE_HOST_SYSTEM}")
message ("Hostname:  ${HOSTNAME}")


find_package (MPI REQUIRED)
message ("The C++ compiler identification is: ${CMAKE_CXX_COMPILER_ID}")
message ("The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
message ("The MPI C++ compiler is: ${MPI_CXX_COMPILER}")
message ("The underlying C++ compiler is: ${CMAKE_CXX_COMPILER}")

# compiler dependent flags
if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel" )
    add_compile_options (-diag-disable 383 -diag-disable 981 -std=c++11)

elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    add_compile_options (-Wall -Werror -Wsign-compare -Wunused-variable -Warray-bounds -DIPPL_RESTRICT_BUG -std=c++11)

elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    if (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.7.3")
        message (FATAL_ERROR "To build OPAL you need GCC version 4.7.3 or greater")
    elseif (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "4.7.2")
       add_compile_options (-Wall -Werror -Wno-reorder -Wno-unused-local-typedefs -Werror=unused-variable -std=c++11)
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
endif ()

find_package (H5Hut REQUIRED)

set (HDF5_USE_STATIC_LIBRARIES ON)
find_package (HDF5 REQUIRED)

find_package (GSL REQUIRED)

set (BOOSTROOT $ENV{BOOST_DIR})
set (Boost_USE_STATIC_LIBS ON)
set (Boost_USE_MULTITHREADED OFF)
set (Boost_USE_STATIC_RUNTIME OFF)
find_package (Boost 1.55.0 REQUIRED COMPONENTS  regex filesystem system)
if (Boost_INCLUDE_DIRS)
    message ("Found boost include dir: ${Boost_INCLUDE_DIR}")
    message ("Found boost library dir: ${Boost_LIBRARY_DIR}")
    message ("Found boost libraries: ${Boost_LIBRARIES}")
    include_directories (${Boost_INCLUDE_DIRS})
endif ()


OPTION (ENABLE_DKS "Enable DKS" OFF)
if (ENABLE_DKS)
   find_package(DKS 1.1.2 REQUIRED HINTS $ENV{DKS_PREFIX} $ENV{DKS_DIR} $ENV{DKS})
   message ("DKS Version ${DKS_VERSION}")
   message ("Found DKS library: ${DKS_LIBRARY}")
   message ("Found DKS include dir: ${DKS_INCLUDE_DIR}")
   set (DKS_PACKAGE_VERSION "\"${DKS_VERSION}\"")
   add_compile_options (-DDKS_MPI)

   ### OpenCL compiler flags ###
   #add_compile_options (-lOpenCL -pthread -DDKS_OPENCL)

   ### CUDA compiler flags ###
   add_compile_options (-DDKS_CUDA)
   
   ### if any accelerator enabled set flag to use DKS ###
   add_compile_options (-DIPPL_DKS -DIPPL_DKS_CUDA -DOPAL_DKS)
endif (ENABLE_DKS)

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
    ENABLE_LANGUAGE (Fortran)

    set ( CCSE_DIR $ENV{BOXLIB_PREFIX} )
    set ( CCSE_LIBRARY_DIR $ENV{BOXLIB_LIBRARY_DIR} )

    find_package (CCSE REQUIRED)

    message ("Found BoxLib include dir: ${CCSE_INCLUDE_DIR}")
    message ("Found BoxLib library dir: ${CCSE_LIBRARY_DIR}")
    message ("Found BoxLib perl dir: ${CCSE_PERL_DIR}")

    set (BL_SPACEDIM 3 CACHE INT "Dimension of BoxLib build")
    set (ENABLE_MPI 1 CACHE INT "Enable build with MPI")
    set (ENABLE_OpenMP 0 CACHE INT "Enable build with OpenMP")
    set (BL_PRECISION "DOUBLE" CACHE INT "Precision of BoxLib build")
    set (BL_USE_PARTICLES 1 CACHE INT "Include Particles classes in BoxLib build")
    set (BL_DEBUG 0 CACHE INT "Include debugging information in BoxLib build")
    set (ENABLE_PROFILING 0 CACHE INT "Include profiling information in BoxLib build")
    set (ENABLE_BACKTRACE 0 CACHE INT "Include backtrace information in BoxLib build")

    # FIXME disable due to BoxLib
    add_compile_options (-Wno-unused-variable -Wno-maybe-uninitialized)

    add_compile_options (-DBL_PRECISION=${BL_PRECISION})
    add_compile_options (-DENABLE_MPI=${ENABLE_MPI})
    add_compile_options (-DENABLE_OpenMP=${ENABLE_OpenMP})
    add_compile_options (-DUSE_PARTICLES=${BL_USE_PARTICLES})
    add_compile_options (-DBL_DEBUG=${BL_DEBUG})
    add_compile_options (-DENABLE_BACKTRACE=${ENABLE_BACKTRACE})
    add_compile_options (-DENABLE_PROFILING=${ENABLE_PROFILING})

    message ("Use following BoxLib settings:")
    message ("  BL_SPACEDIM = ${BL_SPACEDIM}")
    message ("  BL_MACHINE = ${BL_MACHINE}")
    message ("  BL_PRECISION = ${BL_PRECISION}")
    message ("  ENABLE_MPI = ${ENABLE_MPI}")
    message ("  ENABLE_OpenMP = ${ENABLE_OpenMP}")
    message ("  BL_DEBUG = ${BL_DEBUG}")
    message ("  BL_USE_PARTICLES = ${BL_USE_PARTICLES}")
    message ("  ENABLE_BACKTRACE = ${ENABLE_BACKTRACE}")
    message ("  ENABLE_PROFILING = ${ENABLE_PROFILING}")
endif (ENABLE_AMR_SOLVER)

# remove digits from hostname: edison03 -> edison
STRING(REGEX MATCH "[^0-9]*" HOSTNAME_BASE "${HOSTNAME}")

# FIXME: This shouldn't be done like this
if (${HOSTNAME_BASE} MATCHES "daint")
   add_compile_options (-Wno-unused-variable)
endif(${HOSTNAME_BASE} MATCHES "daint")

# In case of Edison and Cori the find_package(Trilinos) gives the wrong
# Trilinos_TPL_LIBRARIES, instead of cray-tpsl/16.07.1 it uses
# cray-tpsl/16.06.1 libraries -> Linker error.
# We can avoid this issue by not going into the if-statement
if (NOT ${HOSTNAME_BASE} MATCHES "edison" AND NOT ${HOSTNAME_BASE} MATCHES "cori" AND (ENABLE_SAAMG_SOLVER OR ENABLE_AMR_SOLVER) )
    message ("Enable SAAMG_SOLVER  " ${ENABLE_SAAMG_SOLVER})
    message ("Enable AMR_SOLVER    " ${ENABLE_AMR_SOLVER})
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
    endif ()
endif ()

if (DBG_SCALARFIELD)
    message ("Write scalar rho_m field is enabled ")
    add_compile_options (-DDBG_SCALARFIELD)
endif (DBG_SCALARFIELD)

if (DBG_STENCIL)
    message ("Write stencil is enabled ")
    add_compile_options (-DDBG_STENCIL)
endif (DBG_STENCIL)

if (BUILD_LIBOPAL)
    message ("Building OPAL library")
endif (BUILD_LIBOPAL)

if (NOCPLUSPLUS11_FOREACH)
    message ("Building OPAL without c++11 foreach support")
    add_compile_options (-DOPAL_NOCPLUSPLUS11_FOREACH)
endif (NOCPLUSPLUS11_FOREACH)

if (NOCPLUSPLUS11_NULLPTR)
    message ("Building OPAL without c++11 nullptr support")
    add_compile_options (-DOPAL_NOCPLUSPLUS11_NULLPTR)
endif (NOCPLUSPLUS11_NULLPTR)

if (NO_FIELD_ASSIGN_OPTIMIZATION)
    message ("Compiler will not try to optimize field assignment")
    add_compile_options (-DdontOPTIMIZE_FIELD_ASSIGNMENT)
endif (NO_FIELD_ASSIGN_OPTIMIZATION)

INCLUDE(FindGit)
if (GIT_FOUND)
    set (OPAL_WC_REVISION 0)
    EXECUTE_PROCESS(
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMAND git rev-parse HEAD
        OUTPUT_VARIABLE GIT_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set (OPAL_WC_REVISION \"${GIT_OUT}\")
else()
    message ("could not find git")
    set (OPAL_WC_REVISION \"\")
endif()

set (IPPL_USED_FROM_OPAL "TRUE")
ADD_SUBDIRECTORY (ippl)
ADD_SUBDIRECTORY (src)

if (BUILD_OPAL_UNIT_TESTS)
     ADD_SUBDIRECTORY(tests)
endif (BUILD_OPAL_UNIT_TESTS)

# make variables visible in other CMakeKists files
set(OPAL_CXX_FLAGS ${OPAL_CXX_FLAGS}
         CACHE INTERNAL "" FORCE
)

set(OPAL_LIBS ${OPAL_LIBS}
         CACHE INTERNAL "" FORCE
)
