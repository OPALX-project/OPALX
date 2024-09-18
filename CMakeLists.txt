cmake_minimum_required (VERSION 3.20.3)
project (OPALX VERSION 2024.1.00 LANGUAGES C CXX)
include(FetchContent)
link_directories(LINK_DIRECTORIES)
set (PROJECT_BUGREPORT opal@lists.psi.ch)

set (OPALX_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
        )
endif ()

if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    add_definitions (-DNOPAssert)
endif ()

if (NOT "${IPPL_PLATFORMS}" STREQUAL "CUDA" AND ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    add_compile_options (-fsanitize=undefined,address)
    add_link_options (-fsanitize=undefined,address)
endif ()

add_compile_options (-Wall)
add_compile_options (-Wunused)
add_compile_options (-Wextra)
# add_compile_options (-Werror)
add_compile_options (-Wno-unused-parameter)
add_compile_options (-Wno-error=unused-variable)

# allow deprecated functions
add_compile_options (-Wno-deprecated-declarations)

option (USE_STATIC_LIBRARIES "Link with static libraries if available" ON)

message (STATUS "Build type is: " ${CMAKE_BUILD_TYPE})

# Select flags.
set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CUDA_EXTENSIONS OFF)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3")
set (CMAKE_CXX_FLAGS_RELEASE "-g -O3 ")
set (CMAKE_CXX_FLAGS_DEBUG  "-g -O0")

if (NOT "${IPPL_PLATFORMS}" STREQUAL "CUDA" AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Use sanitizer flags for GCC or Clang only
     set (CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address,undefined")
endif ()

message (STATUS "CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")

# Disable compile time assert (used by IPPL)
add_definitions (-DNOCTAssert)

#add_compile_options (-ferror-limit=1)
add_compile_options (-Wno-error=uninitialized)
add_compile_options (-Wno-error=restrict)
add_compile_options (-Wno-error=unused-function)
add_compile_options (-Wno-unused-parameter)
add_compile_options (-Wno-error=unused-but-set-variable)
add_compile_options (-Wno-error=maybe-uninitialized)
add_compile_options (-Wno-error=reorder)
add_compile_options (-Wno-cpp)
add_compile_options (-funroll-loops)
add_compile_options (-fstrict-aliasing)
add_compile_options (-DKOKKOS_DEPENDENCE)
add_compile_options (-Wno-return-type)
add_compile_options (-gdwarf-4)                           # avoid dwarf errors on merlin

# Resolve all library dependencies
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

find_package (MPI REQUIRED)
message (STATUS "The C++ compiler identification is: ${CMAKE_CXX_COMPILER_ID}")
message (STATUS "The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
message (STATUS "The MPI C++ compiler is: ${MPI_CXX_COMPILER}")
message (STATUS "The underlying C++ compiler is: ${CMAKE_CXX_COMPILER}")

message (STATUS "Found IPPL_DIR: ${IPPL_DIR}")
if(NOT IPPL_VERSION)
    set(IPPL_VERSION "3.2.0")
    message(STATUS "Defaulting to IPPL-${IPPL_VERSION}")
endif()
FetchContent_Declare(ippl
    #URL https://github.com/IPPL-framework/ippl/archive/refs/tags/IPPL-${IPPL_VERSION}.tar.gz
    GIT_REPOSITORY https://github.com/IPPL-framework/ippl.git
    GIT_SHALLOW True
    DOWNLOAD_EXTRACT_TIMESTAMP True
)
FetchContent_MakeAvailable(ippl)
include_directories("${IPPL_SOURCE_DIR}/src")

message (STATUS "Found Ippl_INC: ${IPPL_INCLUDE_DIR}")
message (STATUS "Found Ippl_LIB: ${IPPL_LIBRARY}")

message (STATUS "FFTW lib(s): ${OPAL_FFTW_LIBS}")  
message (STATUS "MPI lib(s): ${OPAL_MPI_LIBS}")


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
# Create build targets
#
#add_custom_target(always_rebuild ALL
add_custom_target(
    DEPENDS
    src/OPALrevision.h
    )

add_custom_command(
    OUTPUT src/OPALrevision.h
    COMMAND ./gen_OPALrevision
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS
#    always_rebuild
    )

add_subdirectory (src)

include_directories (
    ${CMAKE_SOURCE_DIR}/src
    ${IPPL_INCLUDE_DIR}    
)

if (BUILD_OPAL_UNIT_TESTS)
    add_subdirectory (tests)
endif ()

#
# make variables visible in other CMakeLists files
#
set (OPAL_CXX_FLAGS ${OPAL_CXX_FLAGS}
    CACHE INTERNAL "" FORCE
    )

set (FFTW_LIBS ${FFTW_LIBS}
    CACHE INTERNAL "" FORCE
    )


set (OPAL_LIBS ${OPAL_LIBS}
    CACHE INTERNAL "" FORCE
    )

#
# configure files
#

message (STATUS "Project Name ${PROJECT_NAME}")
message (STATUS "CMAKE_CURRENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}")
message (STATUS "CMAKE_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}")

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
# mode: cmake
# cmake-tab-width: 4
# indent-tabs-mode: nil
# require-final-newline: nil
# End:

