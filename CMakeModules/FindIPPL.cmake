#
# Find IPPL includes and library
#
# IPPL (Independent Parallel Particle Layer) is an object-oriented framework for
# particle based applications in computational science requiring high-performance
# parallel computers. It is a library of C++ classes designed to represent common
# abstractions in these applications. IPPL is based on data-parallel programming
# at the highest abstraction layer. Codes developed on serial workstations port to
# all supported architectures, including parallel processors.
# It can be found at:
#     http://amas.web.psi.ch/tools/IPPL/index.html
#
# IPPL_INCLUDE_DIR - where to find ippl.h
# IPPL_LIBRARY     - qualified libraries to link against.
# IPPL_FOUND       - do not attempt to use if "no" or undefined.

FIND_PATH(IPPL_ROOT Ippl.h
  $ENV{IPPL_ROOT}/src
)

IF(IPPL_ROOT)
    SET( IPPL_FOUND "YES" )
ENDIF(IPPL_ROOT)

IF (IPPL_FOUND)
   IF (NOT IPPL_FIND_QUIETLY)
      MESSAGE(STATUS "Found IPPL root dir: ${IPPL_ROOT}")
ENDIF (NOT IPPL_FIND_QUIETLY)
ELSE (IPPL_FOUND)
   IF (IPPL_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find IPPL!")
   ENDIF (IPPL_FIND_REQUIRED)
ENDIF (IPPL_FOUND)
