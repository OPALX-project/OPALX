#
# Find CLASSIC includes and library

# CLASSIC_INCLUDE_DIR - where to find opal.h
# CLASSIC_LIBRARY     - qualified libraries to link against.
# CLASSIC_FOUND       - do not attempt to use if "no" or undefined.

find_path (CLASSIC_INCLUDE_DIR opal.h
    HINTS $ENV{OPAL_DIR}/include $ENV{OPAL_HOME}/include ENV OPAL_INCLUDE_DIR
    PATHS ENV C_INCLUDE_PATH CPLUS_INCLUDE_PATH
    )

find_library (CLASSIC_LIBRARY CLASSIC
    HINTS $ENV{OPAL_DIR}/lib $ENV{OPAL_HOME}/lib ENV OPAL_LIBRARY_DIR
    PATHS ENV LIBRARY_PATH
)

find_package_handle_standard_args (CLASSIC REQUIRED_VARS CLASSIC_INCLUDE_DIR)
find_package_handle_standard_args (CLASSIC REQUIRED_VARS CLASSIC_LIBRARY)