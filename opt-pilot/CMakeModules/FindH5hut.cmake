#
# Find H5hut includes and library
#
# H5hut
# It can be found at:
#     http://amas.web.psi.ch/tools/H5hut/index.html
#
# H5HUT_INCLUDE_DIR - where to find H5hut.h
# H5HUT_LIBRARY     - qualified libraries to link against.
# H5HUT_FOUND       - do not attempt to use if "no" or undefined.

find_path (H5HUT_INCLUDE_DIR H5hut.h
    PATHS
	$ENV{H5HUT_DIR}/include
	$ENV{H5HUT_HOME}/include
	$ENV{H5HUT_PREFIX}/include
	ENV H5HUT_INCLUDE_DIR C_INCLUDE_PATH
    NO_DEFAULT_PATH
)

find_library (H5HUT_LIBRARY H5hut
    PATHS
	$ENV{H5HUT_DIR}/lib
	$ENV{H5HUT_HOME}/lib
	$ENV{H5HUT_PREFIX}/lib
	ENV H5HUT_LIBRARY_DIR LIBRARY_PATH
    NO_DEFAULT_PATH
)

find_package_handle_standard_args (H5hut REQUIRED_VARS H5HUT_INCLUDE_DIR)
find_package_handle_standard_args (H5hut REQUIRED_VARS H5HUT_LIBRARY)
