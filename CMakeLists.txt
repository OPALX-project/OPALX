cmake_minimum_required (VERSION 3.12)
cmake_policy(VERSION 3.12)

project (OPAL VERSION 2023.9.99)
set (PROJECT_BUGREPORT opal@lists.psi.ch)

string (REGEX REPLACE "\\.([0-9])\\." ".0\\1." PACKAGE_VERSION ${PROJECT_VERSION})
string (REGEX REPLACE "\\.([0-9])$" ".0\\1" PACKAGE_VERSION ${PACKAGE_VERSION})
string (REGEX REPLACE "\\." "" PACKAGE_VERSION ${PACKAGE_VERSION})

message (STATUS PACKAGE_VERSION=${PACKAGE_VERSION})

if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
        )
endif ()

cmake_host_system_information (RESULT HOSTNAME QUERY HOSTNAME)

message (STATUS "Build type is: " ${CMAKE_BUILD_TYPE})
message (STATUS "Host OS System: ${CMAKE_HOST_SYSTEM}")
message (STATUS "Hostname: ${HOSTNAME}")

# Select flags.
set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0 -g ")
set (CMAKE_CXX_FLAGS_RELEASE "-O0 -g")
set (CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

# Disable compile time assert (used by IPPL)
add_definitions (-DNOCTAssert)
add_compile_options (-Wno-nonportable-include-path)
add_compile_options (-ferror-limit=1)
add_compile_options (-Wno-error)
add_compile_options (-Wno-deprecated-declarations)
add_compile_options (-Wno-unused)
add_compile_options (-Wextra)
add_compile_options (-Wno-unused-parameter)
add_compile_options (-funroll-loops)
add_compile_options (-fstrict-aliasing)
add_compile_options (-DKOKKOS_DEPENDENCE)

string(TOUPPER "${CMAKE_BUILD_TYPE}" uppercase_CMAKE_BUILD_TYPE)
if (uppercase_CMAKE_BUILD_TYPE MATCHES "^(DEBUG|RELWITHDEBINFO)$")
    option (BUILD_FOR_PRODUCTION "Disable all assertion checks; Only use for production" OFF)
else()
    option (BUILD_FOR_PRODUCTION "Disable all assertion checks; Only use for production" ON)
endif()

# Resolve all library dependencies
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

find_package (MPI REQUIRED)
message (STATUS "The C++ compiler identification is: ${CMAKE_CXX_COMPILER_ID}")
message (STATUS "The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
message (STATUS "The MPI C++ compiler is: ${MPI_CXX_COMPILER}")
message (STATUS "The underlying C++ compiler is: ${CMAKE_CXX_COMPILER}")

option (ENABLE_OpWenMP "Use hybrid parallelism MPI-OpenMP" OFF)

find_package(Kokkos REQUIRED)
message (STATUS "Found Kokkos_DIR: ${Kokkos_DIR}")

if ("${Kokkos_DEVICES}" MATCHES "CUDA")
    set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -G --diag-suppress=20011")
endif()
find_library(KOKKOS_LIBRARY1 NAMES kokkoscore HINTS /Users/adelmann/OPAL/lib)
find_library(KOKKOS_LIBRARY2 NAMES kokkoscontainers HINTS /Users/adelmann/OPAL/lib)
find_library(KOKKOS_LIBRARY3 NAMES kokkossimd HINTS /Users/adelmann/OPAL/lib)
message (STATUS "Found Kokkos lib: ${KOKKOS_LIBRARY1}:${KOKKOS_LIBRARY2}:${KOKKOS_LIBRARY3}")

find_package(Heffte 2.2.0 REQUIRED)
message (STATUS "Found Heffte_DIR: ${Heffte_DIR}")  
find_library(HEFFTE_LIBRARY heffte HINTS /Users/adelmann/OPAL/lib)
message (STATUS "Found HEFFTE lib: ${HEFFTE_LIBRARY}")

message (STATUS "FFTW lib(s): ${OPAL_FFTW_LIBS}")  

message (STATUS "MPI lib(s): ${OPAL_MPI_LIBS}")

find_package (Ippl REQUIRED)
add_compile_options("-I/Users/adelmann/git/ippl/src")
find_library(IPPL_LIBRARY ippl HINTS ${IPPL_LIBRARY_DIR})
message (STATUS "Found IPPL lib: ${IPPL_LIBRAR}")

find_package (H5Hut REQUIRED)

if (USE_STATIC_LIBRARIES)
    set (HDF5_USE_STATIC_LIBRARIES ON)
endif ()

set (HDF5_PREFER_PARALLEL ON)
find_package (HDF5 REQUIRED)
if(UNIX AND NOT APPLE)
    set (HDF5_LIBRARIES ${HDF5_LIBRARIES} ${CMAKE_DL_LIBS} rt z)
endif ()
message (STATUS "Found HDF5 libraries: ${HDF5_LIBRARIES}")

find_package (GSL REQUIRED)
message (STATUS "Found gsl libraries: ${GSL_LIBRARIES}")

set (BOOSTROOT $ENV{BOOST_DIR})
if (USE_STATIC_LIBRARIES)
    set (Boost_USE_STATIC_LIBS ON)
endif ()
set (Boost_USE_MULTITHREADED OFF)
set (Boost_USE_STATIC_RUNTIME OFF)

find_package (Boost 1.66.0
    REQUIRED COMPONENTS filesystem iostreams regex serialization system timer chrono)
if(UNIX AND NOT APPLE)
    set (Boost_LIBRARIES ${Boost_LIBRARIES} rt)
endif ()
message (STATUS "Found Boost include dirs: ${Boost_INCLUDE_DIRS}")
message (STATUS "Found Boost library dirs: ${Boost_LIBRARY_DIRS}")
message (STATUS "Found Boost libraries: ${Boost_LIBRARIES}")
include_directories (SYSTEM ${Boost_INCLUDE_DIRS})

#
# Resolve all library dependencies
#

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

option (USE_STATIC_LIBRARIES "Link with static libraries if available" ON)
if (NOT ${USE_STATIC_LIBRARIES})
    add_compile_options (-fPIE -fPIC)
    add_link_options (-pie)
endif ()

#
# compiler dependent flags
#

message (STATUS "The underlying C++ compiler ID is: ${CMAKE_CXX_COMPILER_ID}")
if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel" )
    add_compile_options ("-diag-disable 383")
    add_compile_options ("-diag-disable 981")
    if (ENABLE_OpenMP)
        add_compile_options (-openmp)
        add_link_options (-openmp)
    endif ()

elseif (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    add_compile_options (-Warray-bounds-pointer-arithmetic)
    add_compile_options (-Wno-overloaded-virtual)
# Increase template depth for boost
    add_compile_options (-ftemplate-depth=1024)
    # get rid of linking warning with boost. Looks like that boost
    # is using this visibility setting.
    add_compile_options (-fvisibility=hidden)
    add_compile_options (-fvisibility-inlines-hidden)
    # Allow overloaded virtual functions (instances to be fixed in OPAL)
    #add_compile_options (-Wno-overloaded-virtual)
    if (ENABLE_OpenMP)
        add_compile_options (-fopenmp)
        add_link_options (-fopenmp)
    endif ()

elseif (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    # An error in Trilinos/12.18.1 causes
    # -Waggressive-loop-optimizations with -O3 optimization.
    # See issue 587 and corresponding merge request 415.
    # Clang does not support this flag!
    add_compile_options (-fno-aggressive-loop-optimizations)

    if (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.7.3")
        message (FATAL_ERROR "To build OPAL you need GCC version 4.7.3 or greater")
    endif ()
    # Warning flags

    if (ENABLE_OpenMP)
        add_compile_options (-fopenmp)
        add_link_options (-fopenmp)
    endif ()

    # Enables extra error checking in the form of precondition assertion, such
    # as bounds checking and null pointer checks when dereferencing smart pointers
    if (NOT ${BUILD_FOR_PRODUCTION})
        add_compile_options(-Wp,-D_GLIBCXX_ASSERTIONS)
    endif()
else ()
    message (STATUS "Unknown C++ compiler. Please use the GNU or Intel compiler, if you are having problems.")
endif ()

enable_language (Fortran OPTIONAL)

message (STATUS "Fortran compilers: ${CMAKE_Fortran_COMPILER_NAMES}")
if (CMAKE_Fortran_COMPILER)
    message (STATUS "The Fortran compiler identification is: ${CMAKE_Fortran_COMPILER_ID}")
    message (STATUS "The Fortran compiler version is: ${CMAKE_Fortran_COMPILER_VERSION}")
    message (STATUS "The Fortran compiler is: ${CMAKE_Fortran_COMPILER}")
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        set (Fortran_LIBS "-lgfortran")
    endif ()
    # CMake issue. Missing library libmpi_mpifh in MPI_Fortran_LINRARIES.
    # See issue 587 and corresponding merge request 415.
    set (Fortran_LIBS "${Fortran_LIBS} -lmpi_mpifh")
else ()
    message (STATUS "No Fortran compiler found (optional)")
    if (ENABLE_AMR)
        message (FATAL_ERROR "AMR is enabled but no Fortran compiler found!")
    endif ()
endif ()

#
# Handle options
#
option (BUILD_OPAL_UNIT_TESTS "Unit tests" OFF)
if (BUILD_OPAL_UNIT_TESTS)
    find_package (GTest REQUIRED HINT $ENV{GTEST_PREFIX})
    add_definitions (-DWITH_UNIT_TESTS)
endif ()

option (ENABLE_DOXYDOC "compile Doxygen documentation" OFF)
if (ENABLE_DOXYDOC)
    find_package(Doxygen REQUIRED)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY
	)
    add_custom_target (doxydoc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif ()

option (DBG_SCALARFIELD "Enable dump of scalar field rho_m" OFF)
if (DBG_SCALARFIELD)
    message (STATUS "Write scalar rho_m field is enabled ")
    add_definitions (-DDBG_SCALARFIELD)
endif ()

#
#
#

add_custom_target(always_rebuild ALL
    DEPENDS
    src/OPALrevision.h
    )

add_custom_command(
    OUTPUT src/OPALrevision.h
    COMMAND ./gen_OPALrevision
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS
    always_rebuild
    )

add_subdirectory (src)

include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${IPPL_INCLUDE_DIR}    
)

if (BUILD_OPAL_UNIT_TESTS)
    add_subdirectory (tests)
endif ()

# make variables visible in other CMakeLists files
set (OPAL_CXX_FLAGS ${OPAL_CXX_FLAGS}
    CACHE INTERNAL "" FORCE
    )

set (FFTW_LIBS ${FFTW_LIBS}
    CACHE INTERNAL "" FORCE
    )


set (OPAL_LIBS ${OPAL_LIBS}
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
