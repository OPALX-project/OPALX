cmake_minimum_required (VERSION 2.6)

project (CLASSIC)
set (CLASSIC_VERSION_MAJOR 5)
set (CLASSIC_VERSION_MINOR 1.4)

set (PACKAGE \"classic\")
set (PACKAGE_BUGREPORT \"andreas.adelmann@psi.ch\")
set (PACKAGE_NAME \"classic\")
set (PACKAGE_TARNAME \"classic\")
set (PACKAGE_VERSION \"5.1.4\")

configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)


#FIXME: get from Makefile.def (IPPL)
set (CMAKE_CXX_FLAGS
    "-DUSE_H5FED -DPARALLEL_IO ${IPPL_CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS}"
)

set (CLASSIC_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
macro (add_cl_sources)
    file (RELATIVE_PATH _relPath "${CLASSIC_SRC_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND CLASSIC_SRCS "${_relPath}/${_src}")
        else()
            list (APPEND CLASSIC_SRCS "${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SRCS to parent directory
        set (CLASSIC_SRCS ${CLASSIC_SRCS} PARENT_SCOPE)
    endif()
endmacro()

add_subdirectory(AbsBeamline)
add_subdirectory(Algebra)
add_subdirectory(Algorithms)
add_subdirectory(BeamlineCore)
add_subdirectory(BeamlineGeometry)
add_subdirectory(Beamlines)
add_subdirectory(Channels)
add_subdirectory(ComponentWrappers)
add_subdirectory(Construction)
add_subdirectory(Fields/SectorMagneticFieldMap)
add_subdirectory(Fields)
add_subdirectory(Filters)
add_subdirectory(FixedAlgebra)
add_subdirectory(MemoryManagement)
add_subdirectory(Parser)
add_subdirectory(Physics)
add_subdirectory(Solvers)
add_subdirectory(Structure)
add_subdirectory(Utilities)


include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${IPPL_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${GSL_INCLUDE_DIR}
    ${IPPL_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/src  # libclassic should not have any dependencies on libopal
)

add_library( CLASSIC ${CLASSIC_SRCS} )
