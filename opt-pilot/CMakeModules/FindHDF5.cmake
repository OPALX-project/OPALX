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

find_path (HDF5_INCLUDE_DIR hdf5.h
    PATHS
	$ENV{HDF5_DIR}/include
	$ENV{HDF5_HOME}/include
	$ENV{HDF5_PREFIX}/include
	ENV HDF5_INCLUDE_DIR C_INCLUDE_PATH
    NO_DEFAULT_PATH
)

find_library (HDF5_LIBRARY hdf5
    PATHS
	$ENV{HDF5_DIR}/lib
	$ENV{HDF5_HOME}/lib
	$ENV{HDF5_PREFIX}/lib
	ENV HDF5_LIBRARY_DIR LIBRARY_PATH
    NO_DEFAULT_PATH
)

find_package_handle_standard_args (hdf5 REQUIRED_VARS HDF5_INCLUDE_DIR)
find_package_handle_standard_args (hdf5 REQUIRED_VARS HDF5_LIBRARY)

