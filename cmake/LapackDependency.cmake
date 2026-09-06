# Pinned reference BLAS/LAPACK and C interface for small host-side optics maps.
# No system/vendor selection: configure/build is reproducible through FetchContent.
include(FetchContent)
include(CheckLanguage)
check_language(Fortran)
if(NOT CMAKE_Fortran_COMPILER)
  message(FATAL_ERROR
    "Fetched LAPACK requires a Fortran compiler. Install gfortran (or equivalent) "
    "and configure with -DCMAKE_Fortran_COMPILER=/path/to/compiler.")
endif()
enable_language(Fortran)

function(opalx_fetch_lapack)
  # Function scope prevents generic upstream options from changing OPALX/IPPL.
  set(BUILD_TESTING OFF)
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  set(BUILD_SINGLE OFF)
  set(BUILD_DOUBLE ON)
  set(BUILD_COMPLEX OFF)
  set(BUILD_COMPLEX16 OFF)
  set(BUILD_INDEX64 OFF)
  set(BUILD_INDEX64_EXT_API OFF)
  set(BUILD_DEPRECATED OFF)
  set(LAPACKE ON)
  set(LAPACKE_BUILD_SINGLE OFF)
  set(LAPACKE_BUILD_DOUBLE ON)
  set(LAPACKE_BUILD_COMPLEX OFF)
  set(LAPACKE_BUILD_COMPLEX16 OFF)
  set(LAPACKE_WITH_TMG OFF)
  set(CBLAS OFF)
  set(USE_OPTIMIZED_BLAS OFF)
  set(USE_OPTIMIZED_LAPACK OFF)
  set(BLAS_LIBRARIES "")
  set(LAPACK_LIBRARIES "")
  set(BLAS_FOUND FALSE)
  set(LATESTLAPACK_FOUND FALSE)
  set(USE_XBLAS OFF)
  set(BUILD_HTML_DOCUMENTATION OFF)
  set(BUILD_MAN_DOCUMENTATION OFF)
  FetchContent_Declare(opalx_lapack
    URL https://codeload.github.com/Reference-LAPACK/lapack/tar.gz/refs/tags/v3.12.1
    URL_HASH SHA256=2ca6407a001a474d4d4d35f3a61550156050c48016d949f0da0529c0aa052422
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  FetchContent_MakeAvailable(opalx_lapack)
endfunction()
opalx_fetch_lapack()
add_library(OPALX::lapacke ALIAS lapacke)
