# ------------------------------------------------------------------------------
# Determine MPI compilers first
# ------------------------------------------------------------------------------
find_package(MPI REQUIRED)

# Use MPI compilers for all projects
set(CMAKE_C_COMPILER ${MPI_C_COMPILER} CACHE FILEPATH "C compiler" FORCE)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE FILEPATH "C++ compiler" FORCE)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} CACHE FILEPATH "Fortran compiler" FORCE)

# ------------------------------------------------------------------------------
# Project
# ------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.20)
project(opalx LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Determine number of cores
include(ProcessorCount)
ProcessorCount(NPROC)
if(NPROC EQUAL 0)
  set(NPROC 2)
endif()

message (STATUS "The C++ compiler identification is: ${CMAKE_CXX_COMPILER_ID}")
message (STATUS "The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}")
message (STATUS "The MPI C++ compiler is: ${MPI_CXX_COMPILER}")
message (STATUS "The underlying C++ compiler is: ${CMAKE_CXX_COMPILER}")

# ------------------------------------------------------------------------------
# FetchContent module
# ------------------------------------------------------------------------------
include(FetchContent)

# ------------------------------------------------------------------------------
# FFTW
# ------------------------------------------------------------------------------

# FFTW uses autotools (not CMake), so we need an external project build:
include(ExternalProject)
ExternalProject_Add(fftw_external
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE            # avoids CMP0135 dev warning
  URL http://www.fftw.org/fftw-3.3.10.tar.gz
  CONFIGURE_COMMAND ./configure --enable-float --enable-shared --enable-threads --prefix=<INSTALL_DIR>
  BUILD_COMMAND make -j
  INSTALL_COMMAND make install
  BUILD_IN_SOURCE TRUE
)
ExternalProject_Get_Property(fftw_external INSTALL_DIR)
set(FFTW_ROOT ${INSTALL_DIR})
set(FFTW_LIBRARIES ${FFTW_ROOT}/lib/libfftw3f.so)
set(FFTW_INCLUDE_DIR ${FFTW_ROOT}/include)

# ------------------------------------------------------------------------------
# HDF5
# ------------------------------------------------------------------------------

set(HDF5_VERSION 1.10.8)
set(HDF5_PREFIX ${CMAKE_BINARY_DIR}/_deps/hdf5)

ExternalProject_Add(hdf5_external
  URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz
  PREFIX ${HDF5_PREFIX}
  CONFIGURE_COMMAND CC=mpicc CXX=mpicxx ./configure --prefix=<INSTALL_DIR> --disable-shared --enable-parallel --disable-cxx --enable-unsupported --with-pic
  BUILD_COMMAND make -j${NPROC}
  INSTALL_COMMAND make install
  BUILD_IN_SOURCE TRUE
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

ExternalProject_Get_Property(hdf5_external INSTALL_DIR)
set(HDF5_ROOT ${INSTALL_DIR})
set(HDF5_INCLUDE_DIR ${HDF5_ROOT}/include)
set(HDF5_LIBRARIES ${HDF5_ROOT}/lib/libhdf5.a)

# ------------------------------------------------------------------------------
# GSL
# ------------------------------------------------------------------------------

set(GSL_VERSION 2.7)
set(GSL_TAR "gsl-${GSL_VERSION}.tar.gz")
set(GSL_URL "http://ftp.gnu.org/gnu/gsl/${GSL_TAR}")
set(GSL_INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/gsl)
set(GSL_SRC_DIR ${CMAKE_BINARY_DIR}/_deps/gsl-src)

include(ExternalProject)

ExternalProject_Add(gsl_external
    URL ${GSL_URL}
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    PREFIX ${GSL_INSTALL_DIR}
    CONFIGURE_COMMAND
        <SOURCE_DIR>/configure
            --prefix=<INSTALL_DIR>
            --disable-shared
            --enable-static
            CFLAGS=-fPIC
    BUILD_COMMAND make -j${NPROC}
    INSTALL_COMMAND make install
    BUILD_IN_SOURCE TRUE
)

ExternalProject_Get_Property(gsl_external INSTALL_DIR)
set(GSL_ROOT ${INSTALL_DIR})
set(GSL_INCLUDE_DIR ${GSL_ROOT}/include)
set(GSL_LIBRARIES ${GSL_ROOT}/lib/libgsl.a ${GSL_ROOT}/lib/libgslcblas.a)

# ------------------------------------------------------------------------------
# H5hut
# ------------------------------------------------------------------------------

set(H5HUT_VERSION 2.0.0rc6)
set(H5HUT_TAR "H5hut-${H5HUT_VERSION}.tar.gz")
set(H5HUT_URL "http://amas.web.psi.ch/Downloads/H5hut/${H5HUT_TAR}")
set(H5HUT_INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/h5hut)
set(H5HUT_SRC_DIR ${CMAKE_BINARY_DIR}/_deps/h5hut-src)

ExternalProject_Add(h5hut_external
  URL ${H5HUT_URL}
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  PREFIX ${H5HUT_INSTALL_DIR}
  CONFIGURE_COMMAND
    <SOURCE_DIR>/autogen.sh &&
    env CC=mpicc CXX=mpicxx
    <SOURCE_DIR>/configure
      --prefix=<INSTALL_DIR>
      --enable-shared
      --enable-parallel
      --with-pic
      --with-hdf5=${HDF5_ROOT}
      CPPFLAGS=-I${HDF5_INCLUDE_DIR}
      LDFLAGS=-L${HDF5_ROOT}/lib
  BUILD_COMMAND make -j${NPROC}
  INSTALL_COMMAND make install
  BUILD_IN_SOURCE TRUE
)

ExternalProject_Get_Property(h5hut_external INSTALL_DIR)
set(H5HUT_ROOT ${INSTALL_DIR})
set(H5HUT_INCLUDE_DIR ${H5HUT_ROOT}/include)
set(H5HUT_LIBRARIES ${H5HUT_ROOT}/lib/libh5hut.so)  # adjust for .dylib on macOS

# ------------------------------------------------------------------------------
# IPPL library
# ------------------------------------------------------------------------------

if ("${IPPL_PLATFORMS}" STREQUAL "CUDA")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wno-deprecated-gpu-targets")
endif()

# Disable compile time assert (used by IPPL)
add_definitions (-DNOCTAssert)

message(STATUS "Found IPPL_DIR: ${IPPL_DIR}")

if(NOT IPPL_VERSION)
    set(IPPL_VERSION "3.2.0")
    message(STATUS "Defaulting to IPPL-${IPPL_VERSION}")
endif()

# Allow user to specify branch/tag, default to master
set(IPPL_GIT_TAG "master" CACHE STRING "Branch or tag for IPPL (default: master)")
message(STATUS "Fetching IPPL branch/tag: ${IPPL_GIT_TAG}")

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose build type" FORCE)
endif()
message(STATUS "Build type is: ${CMAKE_BUILD_TYPE}")

FetchContent_Declare(
    ippl
    GIT_REPOSITORY https://github.com/IPPL-framework/ippl.git
    GIT_TAG ${IPPL_GIT_TAG}
    GIT_SHALLOW TRUE
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

FetchContent_MakeAvailable(ippl)
message(STATUS "IPPL include path: ${ippl_SOURCE_DIR}/src")

# ------------------------------------------------------------------------------
# Boost library
# ------------------------------------------------------------------------------
include(ExternalProject)

set(BOOST_VERSION "1.82.0")
string(REPLACE "." "_" BOOST_VERSION_UNDERSCORE "${BOOST_VERSION}")
set(BOOST_TAR "boost_${BOOST_VERSION_UNDERSCORE}.tar.gz")
set(BOOST_URL "https://netcologne.dl.sourceforge.net/project/boost/boost/${BOOST_VERSION}/${BOOST_TAR}")

set(BOOST_INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/boost)
set(BOOST_SRC_DIR ${CMAKE_BINARY_DIR}/_deps/boost/src)
file(MAKE_DIRECTORY ${BOOST_SRC_DIR})

ExternalProject_Add(boost_external
    PREFIX ${BOOST_INSTALL_DIR}
    URL ${BOOST_URL}
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    DOWNLOAD_DIR ${BOOST_SRC_DIR}
    SOURCE_DIR ${BOOST_SRC_DIR}/boost_${BOOST_VERSION_UNDERSCORE}

    # Configure (bootstrap)
    CONFIGURE_COMMAND
        <SOURCE_DIR>/bootstrap.sh
        --prefix=<INSTALL_DIR>
        --with-toolset=gcc
        --without-libraries=python,mpi,wave,url,context,coroutine,fiber
    BUILD_COMMAND
        <SOURCE_DIR>/b2
        toolset=gcc
        cxxflags=-fPIC
        cflags=-fPIC
        cxxstd=17
        link=shared,static
        threading=multi
        variant=release
        --prefix=<INSTALL_DIR>
        install -j1
    INSTALL_COMMAND <SOURCE_DIR>/b2 install
    BUILD_IN_SOURCE TRUE
)

# Expose properties for your targets
ExternalProject_Get_Property(boost_external INSTALL_DIR)
set(BOOST_ROOT ${INSTALL_DIR})
set(BOOST_INCLUDE_DIR ${BOOST_ROOT}/include)
set(BOOST_LIBRARIES
    ${BOOST_ROOT}/lib/libboost_system.a
    ${BOOST_ROOT}/lib/libboost_filesystem.a
    ${BOOST_ROOT}/lib/libboost_program_options.a
)
message(STATUS "Boost include dir: ${BOOST_INCLUDE_DIR}")

# ------------------------------------------------------------------------------
# Your OPALX target
# ------------------------------------------------------------------------------

add_subdirectory(src)

add_dependencies(opalx fftw_external hdf5_external h5hut_external boost_external gsl_external ippl)

# Optionally, if IPPL has a library target, you can link it like this:
# target_link_libraries(opalx PRIVATE ippl)
message(STATUS "IPPL include path: ${ippl_SOURCE_DIR}/src")

target_include_directories(opalx PRIVATE
  ${FFTW_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR}
  ${H5HUT_INCLUDE_DIR}
  ${GSL_INCLUDE_DIR}
  "${ippl_SOURCE_DIR}/src"
  ${BOOST_INCLUDE_DIR}
  ${BOOST_SRC_DIR}/boost_${BOOST_VERSION_UNDERSCORE}  # raw source headers
)

target_link_libraries(opalx PRIVATE
  ${FFTW_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${H5HUT_LIBRARIES}
    ${GSL_LIBRARIES}
    ippl
    ${BOOST_LIBRARIES}
  #Boost::filesystem
  #Boost::program_options
  #h5hut
  #gtest
  #gtest_main
)



