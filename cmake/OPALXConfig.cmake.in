@PACKAGE_INIT@

include(CMakeFindDependencyMacro)

find_dependency(MPI REQUIRED)

# -----------------------------------------------------------------------------
# To support FetchContent downloaded kokkos/heffte:
# Detect if this file is being sourced from build tree or install tree
# -----------------------------------------------------------------------------
get_filename_component(_opalxconfig_dir "${CMAKE_CURRENT_LIST_DIR}" REALPATH)
get_filename_component(_build_dir       "@CMAKE_BINARY_DIR@"       REALPATH)
get_filename_component(_install_dir     "@CMAKE_INSTALL_PREFIX@"   REALPATH)
if(_opalxconfig_dir MATCHES "^${_install_dir}")
    set(_INSTALL_TREE 1)
elseif(_opalxconfig_dir MATCHES "^${_build_dir}")
    set(_BUILD_TREE 1)
endif()

# -----------------------------------------------------------------------------
# Make sure OPALX targets are setup
# -----------------------------------------------------------------------------
include("${CMAKE_CURRENT_LIST_DIR}/OPALXTargets.cmake")
