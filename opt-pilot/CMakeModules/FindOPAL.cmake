#
# Find OPAL includes and library
#
# OPAL (Object Oriented Parallel Accelerator Library) is a tool for
# charged-particle optic calculations in large accelerator structures and
# beam lines including 3D space charge. OPAL is built from first principles
# as a parallel application, OPAL admits simulations of any scale: on the
# laptop and up to the largest High Performance Computing (HPC) clusters
# available today. Simulations, in particular HPC simulations, form the
# third pillar of science, complementing theory and experiment.
# OPAL includes various beam line element descriptions and methods for single
# particle optics, namely maps up to arbitrary order, symplectic integration
# schemes and lastly time integration. OPAL is based on IPPL (Independent
# Parallel Particle Layer) which adds parallel capabilities. Main functions
# inherited from IPPL are: structured rectangular grids, fields and parallel
# FFT and particles with the respective interpolation operators. Other
# features are, expression templates and massive parallelism (up to 8000
# processors) which makes is possible to tackle the largest problems in the
# field.
# It can be found at:
#     http://amas.web.psi.ch/docs/index.html
#

# OPAL_INCLUDE_DIR - where to find opal.h
# OPAL_LIBRARY     - qualified libraries to link against.
# OPAL_FOUND       - do not attempt to use if "no" or undefined.

find_path (OPAL_INCLUDE_DIR opal.h
    PATHS
	$ENV{OPAL_DIR}/include
	$ENV{OPAL_HOME}/include
	$ENV{OPAL_PREFIX}/include
	ENV OPAL_INCLUDE_DIR C_INCLUDE_PATH CPLUS_INCLUDE_PATH
    NO_DEFAULT_PATH
    )

find_library (OPAL_LIBRARY OPAL
    PATHS
	$ENV{OPAL_DIR}/lib
	$ENV{OPAL_HOME}/lib
	$ENV{OPAL_PREFIX}/lib
	ENV OPAL_LIBRARY_DIR LIBRARY_PATH
    NO_DEFAULT_PATH
)

find_package_handle_standard_args (OPAL REQUIRED_VARS OPAL_INCLUDE_DIR)
find_package_handle_standard_args (OPAL REQUIRED_VARS OPAL_LIBRARY)

if (OPAL_FOUND)
    set (OPAL_INCLUDE_DIR
        ${OPAL_INCLUDE_DIR}
        ${OPAL_INCLUDE_DIR}/Classic
        )
endif ()
