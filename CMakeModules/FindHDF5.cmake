#
# Find HDF5 includes and library
#
# HDF5 
# It can be found at:
#     http://amas.web.psi.ch/tools/HDF5/index.html
#
# HDF5_INCLUDE_DIR - where to find hdf5.h
# HDF5_LIBRARY     - qualified libraries to link against.
# HDF5_FOUND       - do not attempt to use if "no" or undefined.

SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

FIND_PATH(HDF5_INCLUDE_DIR hdf5.h
  PATHS $ENV{HDF5_INCLUDE_PATH}
  NO_DEFAULT_PATH
)

FIND_PATH(HDF5_INCLUDE_DIR hdf5.h
  /usr/include
  /usr/local/include
)

FIND_LIBRARY(HDF5_LIBRARY libhdf5.a
  PATHS $ENV{HDF5_LIBRARY_PATH}
  NO_DEFAULT_PATH
)

FIND_LIBRARY(HDF5_LIBRARY libhdf5.a
  /usr/lib 
  /usr/local/lib
)

#SET(HDF5_LIBRARY
#    $ENV{HDF5_LIBRARY_PATH}/libhdf5.a
#)

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
