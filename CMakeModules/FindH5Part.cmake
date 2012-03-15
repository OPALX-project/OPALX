#
# Find H5Part includes and library
#
# H5Part 
# It can be found at:
#     http://amas.web.psi.ch/tools/H5Part/index.html
#
# H5Part_INCLUDE_DIR - where to find ippl.h
# H5Part_LIBRARY     - qualified libraries to link against.
# H5Part_FOUND       - do not attempt to use if "no" or undefined.

FIND_PATH(H5Part_INCLUDE_DIR H5Part.h
  /usr/local/include
  /usr/include
  $ENV{H5Part}/src
  $ENV{H5Part}/include
)

FIND_LIBRARY(H5Part_LIBRARY pH5Part
  /usr/local/lib
  /usr/lib
  $ENV{H5Part}/src
  $ENV{H5Part}/lib
)

IF(H5Part_INCLUDE_DIR AND H5Part_LIBRARY)
    SET( H5Part_FOUND "YES" )
ENDIF(H5Part_INCLUDE_DIR AND H5Part_LIBRARY)

IF (H5Part_FOUND)
   IF (NOT H5Part_FIND_QUIETLY)
      MESSAGE(STATUS "Found H5Part: ${H5Part_LIBRARY}")
   ENDIF (NOT H5Part_FIND_QUIETLY)
ELSE (H5Part_FOUND)
   IF (H5Part_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find H5Part!")
   ENDIF (H5Part_FIND_REQUIRED)
ENDIF (H5Part_FOUND)
