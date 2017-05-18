#
# Find GSL includes and library
#
# GSL
# It can be found at:
#     http://amas.web.psi.ch/tools/GSL/index.html
#
# GSL_INCLUDE_DIR - where to find ippl.h
# GSL_LIBRARY     - qualified libraries to link against.
# GSL_FOUND       - do not attempt to use if "no" or undefined.

FIND_PATH(GSL_INCLUDE_DIR gsl/gsl_fft.h
    HINTS $ENV{GSL_INCLUDE_PATH} $ENV{GSL_INCLUDE_DIR} $ENV{GSL_PREFIX}/include $ENV{GSL_DIR}/include $ENV{GSL}/include
    PATHS ENV CPP_INCLUDE_PATH
)


FIND_LIBRARY(GSL_LIBRARY libgsl.a gsl
    HINTS $ENV{GSL_LIBRARY_PATH} $ENV{GSL_LIBRARY_DIR} $ENV{GSL_PREFIX}/lib $ENV{GSL_DIR}/lib $ENV{GSL}/lib
    PATHS ENV LIBRARY_PATH
)

FIND_LIBRARY(GSL_LIBRARY_CBLAS libgslcblas.a gslcblas
    HINTS $ENV{GSL_LIBRARY_PATH} $ENV{GSL_LIBRARY_DIR} $ENV{GSL_PREFIX}/lib $ENV{GSL_DIR}/lib $ENV{GSL}/lib
    PATHS ENV LIBRARY_PATH
)

set( GSL_LIBRARY
    ${GSL_LIBRARY}
    ${GSL_LIBRARY_CBLAS}
)

IF(GSL_INCLUDE_DIR AND GSL_LIBRARY)
    SET( GSL_FOUND "YES" )
ENDIF(GSL_INCLUDE_DIR AND GSL_LIBRARY)

IF (GSL_FOUND)
   IF (NOT GSL_FIND_QUIETLY)
      MESSAGE(STATUS "Found GSL libraries: ${GSL_LIBRARY}")
      MESSAGE(STATUS "Found GSL include dir: ${GSL_INCLUDE_DIR}")
   ENDIF (NOT GSL_FIND_QUIETLY)
ELSE (GSL_FOUND)
   IF (GSL_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find GSL!")
   ENDIF (GSL_FIND_REQUIRED)
ENDIF (GSL_FOUND)
