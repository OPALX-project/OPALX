#
# Find IPPL package
# https://
#
# IPPL_INCLUDE_DIR
# IPPL_LIBRARY_DIR
# IPPL_FOUND

message (STATUS "Hint IPPL include dir: $ENV{IPPL_INCLUDE_PATH}")

find_path (IPPL_INCLUDE_DIR Ippl.h
    HINTS $ENV{IPPL_INCLUDE_PATH} 
    PATHS ENV C_INCLUDE_PATH 
)

find_path (IPPL_LIBRARY_DIR libippl.a
    HINTS $ENV{IPPL_LIBRARY_PATH}
    PATHS ENV LIBRARY_PATH
)

if (IPPL_INCLUDE_DIR AND IPPL_LIBRARY_DIR)
    set (IPPL_FOUND "YES")
endif ()

if (IPPL_FOUND)
    if (NOT IPPL_FIND_QUIETLY)
        message (STATUS "Found IPPL include dir: ${IPPL_INCLUDE_DIR}")
        message (STATUS "Found IPPL library dir: ${IPPL_LIBRARY_DIR}")
    endif ()
else (IPPL_FOUND)
    if (IPPL_FIND_REQUIRED)
        if (NOT IPPL_INCLUDE_DIR)
            message (WARNING 
                "IPPL include directory was not found! "
                "Make sure that IPPL is compiled and that "
                "the directory ippl/include/ippl has been automatically created. "
                "Also make sure that at least one of the following "
                "environment variables is set: "
                "IPPL_INCLUDE_DIR, IPPL_INCLUDE_PATH, IPPL_PREFIX, or IPPL.")
        endif ()
        if (NOT IPPL_LIBRARY_DIR)
            message (WARNING 
                "IPPL library was not found! "
                "Make sure that IPPL is compiled and that "
                "the directory ippl/lib has been automatically created. "
                "Also make sure that at least one of the following "
                "environment variables is set: "
                "IPPL_LIBRARY_DIR, IPPL_LIBRARY_PATH, IPPL_PREFIX, or IPPL.")
        endif ()
        message (STATUS "IPPL can be downloaded and compiled from https://xxxx")
        message (FATAL_ERROR "Could not find IPPL!")
    endif (IPPL_FIND_REQUIRED)
endif (IPPL_FOUND)


include (CheckIncludeFile)
SET (CMAKE_REQUIRED_INCLUDES ${IPPL_INCLUDE_DIR})
CHECK_INCLUDE_FILE (Ippl.h HAVE_API2_FUNCTIONS "-I${IPPL_INCLUDE_DIR} -DPARALLEL_IO")

IF (HAVE_API2_FUNCTIONS)
    MESSAGE (STATUS "IPPL version is OK")
ELSE (HAVE_API2_FUNCTIONS)
    MESSAGE (ERROR "IPPL >= 2 required")
ENDIF (HAVE_API2_FUNCTIONS)