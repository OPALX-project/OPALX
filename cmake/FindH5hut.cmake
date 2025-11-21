# --- FindH5hut.cmake -----------------------------------------------------
# Looks for H5hut installed by Spack or system package managers.
# Defines:
#   H5hut_FOUND
#   H5HUT_INCLUDE_DIRS
#   H5HUT_LIBRARIES
#   H5hut::H5hut

find_path(H5HUT_INCLUDE_DIRS
    NAMES H5hut.h
    PATH_SUFFIXES include
)

find_library(H5HUT_LIBRARIES
    NAMES H5hut libH5hut h5hut libh5hut
    PATH_SUFFIXES lib lib64
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    H5hut
    REQUIRED_VARS H5HUT_INCLUDE_DIRS H5HUT_LIBRARIES
)

if(H5hut_FOUND)
    add_library(H5hut::H5hut UNKNOWN IMPORTED)
    set_target_properties(H5hut::H5hut PROPERTIES
        IMPORTED_LOCATION "${H5HUT_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${H5HUT_INCLUDE_DIRS}"
    )
endif()
