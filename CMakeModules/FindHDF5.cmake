#
# Find HDF5 includes and library
#
# HDF5 
# It can be found at:
#     http://amas.web.psi.ch/tools/HDF5/index.html
#
# HDF5_INCLUDE_DIR - where to find ippl.h
# HDF5_LIBRARY     - qualified libraries to link against.
# HDF5_FOUND       - do not attempt to use if "no" or undefined.

FIND_PATH(HDF5_INCLUDE_DIR hdf5.h
  /usr/local/include
  /usr/include
  $ENV{HDF5_INCLUDE_PATH}
)

FIND_LIBRARY(HDF5_LIBRARY hdf5
  /usr/local/lib
  /usr/lib
  $ENV{HDF5_LIBRARY_PATH}
)

IF(HDF5_INCLUDE_DIR AND HDF5_LIBRARY)
    SET( HDF5_FOUND "YES" )
ENDIF(HDF5_INCLUDE_DIR AND HDF5_LIBRARY)

IF (HDF5_FOUND)
   IF (NOT HDF5_FIND_QUIETLY)
      MESSAGE(STATUS "Found HDF5: ${HDF5_LIBRARY}")
   ENDIF (NOT HDF5_FIND_QUIETLY)
ELSE (HDF5_FOUND)
   IF (HDF5_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find HDF5!")
   ENDIF (HDF5_FIND_REQUIRED)
ENDIF (HDF5_FOUND)
