#
# Find GSL includes and library
#
# GSL
# It can be found at:
#     http://amas.web.psi.ch/tools/GSL/index.html
#
# GSL_INCLUDE_DIR - where to find ippl.h
# GSL_LIBRARY     - qualified libraries to link against.
# GSL_FOUND       - do not attempt to use if "no" or undefined.

find_path (GSL_INCLUDE_DIR gsl/gsl_fft.h
    PATHS
	$ENV{GSL_DIR}/include
	$ENV{GSL_HOME}/include
	$ENV{GSL_PREFIX}/include
	ENV GSL_INLCUDE_DIR
    NO_DEFAULT_PATH
    )

find_library (GSL_LIBRARY gsl
    PATHS
	$ENV{GSL_DIR}/lib
 	$ENV{GSL_HOME}/lib
 	$ENV{GSL_PREFIX}/lib
	ENV  GSL_LIBRARY_DIR LIBRARY_PATH
    NO_DEFAULT_PATH
)

find_library (GSL_LIBRARY_CBLAS gslcblas
    PATHS
	$ENV{GSL_DIR}/lib
	$ENV{GSL_HOME}/lib
	$ENV{GSL_PREFIX}/lib
	ENV  GSL_LIBRARY_DIR LIBRARY_PATH
    NO_DEFAULT_PATH
)

find_package_handle_standard_args (GSL REQUIRED_VARS GSL_INCLUDE_DIR)
find_package_handle_standard_args (GSL REQUIRED_VARS GSL_LIBRARY)
find_package_handle_standard_args (GSL REQUIRED_VARS GSL_LIBRARY_CBLAS)
set (GSL_LIBRARY ${GSL_LIBRARY} ${GSL_LIBRARY_CBLAS})