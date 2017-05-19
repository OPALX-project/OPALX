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
    HINTS $ENV{H5HUT_INCLUDE_PATH} $ENV{H5HUT_INCLUDE_DIR} $ENV{H5HUT}/include $ENV{H5HUT_PREFIX}/include $ENV{H5HUT_DIR}/include $ENV{H5hut}/include
    PATHS ENV C_INCLUDE_PATH
)

FIND_LIBRARY(H5Hut_LIBRARY libH5hut.a H5hut
    HINTS $ENV{H5HUT_LIBRARY_PATH} $ENV{H5HUT_LIBRARY_DIR} $ENV{H5HUT}/lib $ENV{H5HUT_PREFIX}/lib $ENV{H5HUT_DIR}/lib $ENV{H5hut}/lib
    PATHS ENV LIBRARY_PATH
)

IF(H5Hut_INCLUDE_DIR AND H5Hut_LIBRARY)
    SET( H5Hut_FOUND "YES" )
ENDIF(H5Hut_INCLUDE_DIR AND H5Hut_LIBRARY)

IF (H5Hut_FOUND)
    IF (NOT H5Hut_FIND_QUIETLY)
        MESSAGE(STATUS "Found H5Hut libraries: ${H5Hut_LIBRARY}")
        MESSAGE(STATUS "Found H5Hut include dir: ${H5Hut_INCLUDE_DIR}")
    ENDIF (NOT H5Hut_FIND_QUIETLY)
ELSE (H5Hut_FOUND)
    IF (H5Hut_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find H5Hut!")
    ENDIF (H5Hut_FIND_REQUIRED)
ENDIF (H5Hut_FOUND)

include (CheckIncludeFile)
CHECK_INCLUDE_FILE (H5_file_attribs.h HAVE_API2_FUNCTIONS "-DPARALLEL_IO")

IF (HAVE_API2_FUNCTIONS)
    MESSAGE (STATUS "H5hut version is >= 2")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_H5HUT2")
ELSE (HAVE_API2_FUNCTIONS)
    MESSAGE (STATUS "H5hut version is 1.x")
ENDIF (HAVE_API2_FUNCTIONS)

# Local Variables:
# mode:cmake
# cmake-tab-width: 4
# indent-tabs-mode:nil
# End:
