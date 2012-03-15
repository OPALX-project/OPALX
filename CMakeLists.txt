cmake_minimum_required(VERSION 2.6)

project(OPAL)
set (OPAL_VERSION_MAJOR 1)
set (OPAL_VERSION_MINOR 1.7)

option (ENABLE_ML_SOLVER "Enable iteartive SA-AMG-PCG self field solver" OFF)

option (DBG_SCALARFIELD "Enable dump of scalar field rho_m" OFF)
option (DBG_STENCIL "Dump stencil of iterative solver in matlab format" OFF)

find_package (MPI REQUIRED)
#set (CMAKE_CXX_COMPILER ${MPI_COMPILER})

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")
find_package (IPPL REQUIRED)
find_package (H5Part REQUIRED)
find_package (HDF5 REQUIRED)
find_package (GSL REQUIRED)

if (DBG_SCALARFIELD)	
   MESSAGE("\nWrite scalar rho_m field is enabled ")
   set (CMAKE_CXX_FLAGS_RELEASE "-D DBG_SCALARFIELD") 
endif (DBG_SCALARFIELD)


if (DBG_STENCIL)
   MESSAGE("\nWrite stencil is enabled ")
   set (CMAKE_CXX_FLAGS_RELEASE "-D DBG_STENCIL") 
endif (DBG_STENCIL)


if (ENABLE_ML_SOLVER)
    # Get Trilinos as one entity
    find_package (Trilinos PATHS $ENV{TRILINOS_INCLUDE_PATH})

    IF(NOT Trilinos_FOUND)
        MESSAGE (FATAL_ERROR "Could not find Trilinos!")
    ENDIF()
    # Echo trilinos build info just for fun
    MESSAGE("\nFound Trilinos!  Here are the details: ")
    MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
    MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
    MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
    MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
    MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
    MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
    MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
    MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
    MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
    MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
    MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
    MESSAGE("End of Trilinos details\n")

    # Make sure to use same compilers and flags as Trilinos
    #SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
    SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
    SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

    #SET(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
    #SET(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
    #SET(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")
    
endif (ENABLE_ML_SOLVER)

set (CMAKE_CXX_FLAGS_RELEASE "-g -O2")
set (CMAKE_CXX_FLAGS_DEBUG  "-g -O2")

add_subdirectory (doom)
add_subdirectory (classic/5.0/src)
add_subdirectory (src)
