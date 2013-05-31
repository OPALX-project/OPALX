#
# Find H5hut includes and library
#
# H5Hut
# It can be found at:
#     http://amas.web.psi.ch/tools/H5hut/index.html
#
# H5Hut_INCLUDE_DIR - where to find H5hut.h
# H5Hut_LIBRARY     - qualified libraries to link against.
# H5Hut_FOUND       - do not attempt to use if "no" or undefined.

FIND_PATH(H5Hut_INCLUDE_DIR H5hut.h
    HINTS $ENV{H5hut}/src $ENV{H5hut}/include
    PATHS ENV C_INCLUDE_PATH
)

FIND_LIBRARY(H5Hut_LIBRARY H5hut
    HINTS $ENV{H5hut}/lib
    PATHS ENV LIBRARY_PATH
)

FIND_LIBRARY(H5Hut_LIBRARY_C H5hutC
    HINTS $ENV{H5hut}/lib
    PATHS ENV LIBRARY_PATH
)


IF(H5Hut_INCLUDE_DIR AND H5Hut_LIBRARY)
    SET( H5Hut_FOUND "YES" )
ENDIF(H5Hut_INCLUDE_DIR AND H5Hut_LIBRARY)

IF (H5Hut_FOUND)
   IF (NOT H5Hut_FIND_QUIETLY)
      MESSAGE(STATUS "Found H5Hut libraries: ${H5Hut_LIBRARY}; ${H5Hut_LIBRARY_C}")
      MESSAGE(STATUS "Found H5Hut include dir: ${H5Hut_INCLUDE_DIR}")
   ENDIF (NOT H5Hut_FIND_QUIETLY)
ELSE (H5Hut_FOUND)
   IF (H5Hut_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find H5Hut!")
   ENDIF (H5Hut_FIND_REQUIRED)
ENDIF (H5Hut_FOUND)
