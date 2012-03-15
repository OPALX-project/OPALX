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
  PATHS $ENV{GSL_PREFIX}/include
  NO_DEFAULT_PATH
)
FIND_PATH(GSL_INCLUDE_DIR gsl/gsl_fft.h
  /usr/include
  /usr/local/include
)

FIND_LIBRARY(GSL_LIBRARY gsl
  PATHS $ENV{GSL_PREFIX}/lib
  NO_DEFAULT_PATH
)
FIND_LIBRARY(GSL_LIBRARY gsl
  /usr/lib
  /usr/local/lib
)
FIND_LIBRARY(GSL_LIBRARY_CBLAS gslcblas 
  PATHS $ENV{GSL_PREFIX}/lib
  NO_DEFAULT_PATH
)
FIND_LIBRARY(GSL_LIBRARY_CBLAS gslcblas 
  /usr/lib
  /usr/local/lib
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
      MESSAGE(STATUS "Found GSL: ${GSL_LIBRARY}")
   ENDIF (NOT GSL_FIND_QUIETLY)
ELSE (GSL_FOUND)
   IF (GSL_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find GSL!")
   ENDIF (GSL_FIND_REQUIRED)
ENDIF (GSL_FOUND)
