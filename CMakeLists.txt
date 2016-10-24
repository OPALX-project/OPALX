CMAKE_MINIMUM_REQUIRED (VERSION 2.8.10)
PROJECT (OPAL)
SET (OPAL_VERSION_MAJOR 1)
SET (OPAL_VERSION_MINOR 5.00.1)
set (PACKAGE \"opal\")
set (PACKAGE_BUGREPORT \"opal@lists.psi.ch\")
set (PACKAGE_NAME \"OPAL\")
set (PACKAGE_TARNAME \"opal\")
set (PACKAGE_VERSION "\"${OPAL_VERSION_MAJOR}.${OPAL_VERSION_MINOR}\"")
set (DKS_VERSION \"1.0.1\")

IF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    SET (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
	)
ENDIF (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

MESSAGE (STATUS "Build type is: " ${CMAKE_BUILD_TYPE})
MESSAGE (STATUS "Host OS System: ${CMAKE_HOST_SYSTEM}")
MESSAGE (STATUS "Hostname:  $ENV{HOSTNAME}")

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


### Uldis Locans 2014.09.05 ###
OPTION (ENABLE_OPENCL "Enable OpenCL" OFF)
OPTION (ENABLE_CUDA "Enable CUDA" OFF)
OPTION (ENABLE_OPENMP "Enable OpenMP + offload" OFF)

OPTION (USE_OPENCL "Use OpenCL" OFF)
OPTION (USE_CUDA "Use CUDA" OFF)
OPTION (USE_MIC "Use MIC" OFF)

### Load DKS package ###
IF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_MIC)
   FIND_PACKAGE(DKS REQUIRED HINTS $ENV{DKS_PREFIX} $ENV{DKS_DIR} $ENV{DKS})
   MESSAGE (STATUS "Found DKS library: ${DKS_LIBRARY}")
   MESSAGE (STATUS "Found DKS include dir: ${DKS_INCLUDE_DIR}")
   SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDKS_MPI")
ENDIF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_MIC)

### OpenCL compiler flags ###
IF (ENABLE_OPENCL)
	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lOpenCL -pthread -DDKS_OPENCL")
ENDIF (ENABLE_OPENCL)
### End OpenCL ###

### CUDA compiler flags ###
IF (ENABLE_CUDA)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lcudart -lcufft -lcublas -lnvToolsExt -DDKS_CUDA")
ENDIF (ENABLE_CUDA)
### End Cuda ###

### MIC compiler flags ###
IF (ENABLE_OPENMP)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp -DDKS_MIC")
ENDIF (ENABLE_OPENMP)
### END MIC ###

### if any accelerator enabled set flag to use DKS ###
IF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_OPENMP)
	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPPL_DKS -DOPAL_DKS")
ENDIF (ENABLE_OPENCL OR ENABLE_CUDA OR ENABLE_OPENMP)

### select API for DKS (temporary - will be moved to runtime) ###
IF (USE_OPENCL)
	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPPL_DKS_OPENCL")
ENDIF (USE_OPENCL)
IF (USE_CUDA)
	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPPL_DKS_CUDA")
ENDIF (USE_CUDA)
IF (USE_MIC)
	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPPL_DKS_MIC")
ENDIF (USE_MIC)
### END DKS part ###

# Handle options
OPTION (ENABLE_AMR_SOLVER "Enable BoxLib based AMR solver" OFF)

OPTION (ENABLE_SAAMG_SOLVER "Enable iteartive SA-AMG-PCG self field solver" OFF)

OPTION (DBG_SCALARFIELD "Enable dump of scalar field rho_m" OFF)
OPTION (DBG_STENCIL "Dump stencil of iterative solver in matlab format" OFF)
OPTION (BUILD_LIBOPAL "Enable building the OPAL library" OFF)

OPTION (NOCPLUSPLUS11_FOREACH "Disable C++11 foreach support" OFF)
OPTION (NOCPLUSPLUS11_NULLPTR "Disable C++11 nullptr support" OFF)

OPTION (NO_FIELD_ASSIGN_OPTIMIZATION "Disable compiler optimization of IPPL field assignment" OFF)

IF (USE_H5HUT2)
	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_H5HUT2")
ENDIF (USE_H5HUT2)

IF (ENABLE_AMR_SOLVER)
    ENABLE_LANGUAGE (Fortran)

    SET ( CCSE_DIR $ENV{BOXLIB_PREFIX} )
    SET ( CCSE_LIBRARY_DIR $ENV{BOXLIB_LIBRARY_DIR} )

    FIND_PACKAGE (CCSE REQUIRED)

    MESSAGE (STATUS "Found BoxLib include dir: ${CCSE_INCLUDE_DIR}")
    MESSAGE (STATUS "Found BoxLib library dir: ${CCSE_LIBRARY_DIR}")
    MESSAGE (STATUS "Found BoxLib perl dir: ${CCSE_PERL_DIR}")

    SET (BL_SPACEDIM 3 CACHE INT "Dimension of BoxLib build")
    SET (ENABLE_MPI 1 CACHE INT "Enable build with MPI")
    SET (ENABLE_OpenMP 0 CACHE INT "Enable build with OpenMP")
    SET (BL_PRECISION "DOUBLE" CACHE INT "Precision of BoxLib build")
    SET (BL_USE_PARTICLES 1 CACHE INT "Include Particles classes in BoxLib build")
    SET (BL_DEBUG 0 CACHE INT "Include debugging information in BoxLib build")
    SET (ENABLE_PROFILING 0 CACHE INT "Include profiling information in BoxLib build")
    SET (ENABLE_BACKTRACE 0 CACHE INT "Include backtrace information in BoxLib build")

    # FIXME disable due to BoxLib
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable -Wno-maybe-uninitialized")

    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBL_PRECISION=${BL_PRECISION}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_MPI=${ENABLE_MPI}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_OpenMP=${ENABLE_OpenMP}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_PARTICLES=${BL_USE_PARTICLES}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBL_DEBUG=${BL_DEBUG}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_BACKTRACE=${ENABLE_BACKTRACE}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_PROFILING=${ENABLE_PROFILING}")

    MESSAGE (STATUS "Use following BoxLib settings:")
    MESSAGE (STATUS "  BL_SPACEDIM = ${BL_SPACEDIM}")
    MESSAGE (STATUS "  BL_MACHINE = ${BL_MACHINE}")
    MESSAGE (STATUS "  BL_PRECISION = ${BL_PRECISION}")
    MESSAGE (STATUS "  ENABLE_MPI = ${ENABLE_MPI}")
    MESSAGE (STATUS "  ENABLE_OpenMP = ${ENABLE_OpenMP}")
    MESSAGE (STATUS "  BL_DEBUG = ${BL_DEBUG}")
    MESSAGE (STATUS "  BL_USE_PARTICLES = ${BL_USE_PARTICLES}")
    MESSAGE (STATUS "  ENABLE_BACKTRACE = ${ENABLE_BACKTRACE}")
    MESSAGE (STATUS "  ENABLE_PROFILING = ${ENABLE_PROFILING}")
ENDIF (ENABLE_AMR_SOLVER)

if (ENABLE_SAAMG_SOLVER OR ENABLE_AMR_SOLVER)
    MESSAGE ("Enable SAAMG_SOLVER  " ${ENABLE_SAAMG_SOLVER})
    MESSAGE ("Enable AMR_SOLVER    " ${ENABLE_AMR_SOLVER})
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

IF (AMR_DUMMY_SOLVE)
    MESSAGE (STATUS "\nSolve Poisson equation with rhs equal 1")
    SET (CMAKE_CXX_FLAGS "-DAMR_DUMMY_SOLVE ${CMAKE_CXX_FLAGS}")
ENDIF (AMR_DUMMY_SOLVE)

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

SET (IPPL_USED_FROM_OPAL "TRUE")
ADD_SUBDIRECTORY (ippl)
ADD_SUBDIRECTORY (src)

IF (BUILD_OPAL_UNIT_TESTS)
     ADD_SUBDIRECTORY(tests)
ENDIF (BUILD_OPAL_UNIT_TESTS)

# make variables visible in other CMakeKists files
set(OPAL_CXX_FLAGS ${OPAL_CXX_FLAGS}
         CACHE INTERNAL "" FORCE
)

set(OPAL_LIBS ${OPAL_LIBS}
         CACHE INTERNAL "" FORCE
)