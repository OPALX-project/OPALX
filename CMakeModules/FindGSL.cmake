#
# Find GSL includes and library
#
# FindGSL() shipped with CMake is somehow broken in newer versions (>= 3.9,
# maybe also 3.7 and 3.8). It works only if GSL_ROOT_DIR is set. Hints do
# not work.
#
# :FIXME: For the time being we use our own version.
#
# GSL_INCLUDE_DIR - where to find ippl.h
# GSL_LIBRARY     - GSL library to link against.
# GSL_CBLAS_LIBRARY - GSL CBlas library to link against
# GSL_FOUND       - do not attempt to use if "no" or undefined.

FIND_PATH (GSL_INCLUDE_DIR gsl/gsl_fft.h
    HINTS $ENV{GSL_INCLUDE_PATH} $ENV{GSL_INCLUDE_DIR} $ENV{GSL_PREFIX}/include $ENV{GSL_DIR}/include $ENV{GSL}/include
    PATHS ENV CPP_INCLUDE_PATH
)


FIND_LIBRARY (GSL_LIBRARY gsl
    HINTS $ENV{GSL_LIBRARY_PATH} $ENV{GSL_LIBRARY_DIR} $ENV{GSL_PREFIX}/lib $ENV{GSL_DIR}/lib $ENV{GSL}/lib
    PATHS ENV LIBRARY_PATH
)

FIND_LIBRARY (GSL_CBLAS_LIBRARY gslcblas
    HINTS $ENV{GSL_LIBRARY_PATH} $ENV{GSL_LIBRARY_DIR} $ENV{GSL_PREFIX}/lib $ENV{GSL_DIR}/lib $ENV{GSL}/lib
    PATHS ENV LIBRARY_PATH
)

IF (GSL_INCLUDE_DIR AND GSL_LIBRARY)
    SET( GSL_FOUND "YES" )
ENDIF()

IF (GSL_FOUND)
   IF (NOT GSL_FIND_QUIETLY)
      MESSAGE(STATUS "Found GSL libraries: ${GSL_LIBRARY}")
      MESSAGE(STATUS "Found GSL include dir: ${GSL_INCLUDE_DIR}")
   ENDIF (NOT GSL_FIND_QUIETLY)
ELSE ()
   IF (GSL_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find GSL!")
   ENDIF ()
ENDIF ()
