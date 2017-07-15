cmake_minimum_required (VERSION 2.8.10)
set (DKS_VERSION_OPAL \"${DKS_VERSION}\")

set (CMAKE_CXX_FLAGS
    "${IPPL_CMAKE_CXX_FLAGS} -DUSEH5FEDV2 -DPARALLEL_IO ${CMAKE_CXX_FLAGS}"
    )

IF(ENABLE_SAAMG_SOLVER)
    set (CMAKE_CXX_FLAGS
        "-DHAVE_SAAMG_SOLVER ${CMAKE_CXX_FLAGS}"
    )
ENDIF(ENABLE_SAAMG_SOLVER)

IF(ENABLE_AMR_SOLVER)
    # FIXME Find a better solution!
    # We get otherwise a "undefined reference" to "getLevelBld()" and "probinit_"
    set (CBOXLIB -lcboxlib)

    IF(NOT ${HOSTNAME_BASE} MATCHES "edison" AND NOT ${HOSTNAME_BASE} MATCHES "cori" AND NOT ${HOSTNAME_BASE} MATCHES "daint")
        set (CBOXLIB "-lmpi -lmpi_mpifh ${CBOXLIB}")
    ENDIF(NOT ${HOSTNAME_BASE} MATCHES "edison" AND NOT ${HOSTNAME_BASE} MATCHES "cori" AND NOT ${HOSTNAME_BASE} MATCHES "daint")

ENDIF(ENABLE_AMR_SOLVER)

SET (PACKAGE_COMPILE_OPTIONS \"${CMAKE_CXX_FLAGS}\")
configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/OPALconfig.h)
configure_file(revision.h.in ${CMAKE_CURRENT_SOURCE_DIR}/OPALrevision.h)

SET (OPAL_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET (OPAL_SRCS )
macro (ADD_OPAL_SOURCES)
    file (RELATIVE_PATH _relPath "${OPAL_SRC_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND OPAL_SRCS "${_relPath}/${_src}")
        else()
            list (APPEND OPAL_SRCS "${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SRCS to parent directory
        set (OPAL_SRCS ${OPAL_SRCS} PARENT_SCOPE)
    endif()
endmacro()

MACRO (ADD_OPAL_HEADERS )
    FILE (RELATIVE_PATH _relPath "${OPAL_SRC_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    FOREACH (_hdr ${ARGN})
        IF (_relPath)
            LIST (APPEND OPAL_HDRS "${_relPath}/${_hdr}")
        ELSE ()
            LIST (APPEND OPAL_HDRS "${_hdr}")
        ENDIF ()
    ENDFOREACH ()
    IF (_relPath)
        # propagate HDRS to parent directory
        SET (OPAL_HDRS ${OPAL_HDRS} PARENT_SCOPE)
    ENDIF ()
ENDMACRO ()

add_subdirectory (AbstractObjects)
add_subdirectory (Algebra)
add_subdirectory (Algorithms)
add_subdirectory (Aperture)
add_subdirectory (Attributes)
add_subdirectory (BasicActions)
add_subdirectory (Classic)
add_subdirectory (Distribution)
add_subdirectory (Editor)
add_subdirectory (Elements)
add_subdirectory (Errors)
add_subdirectory (Expressions)
add_subdirectory (Lines)
add_subdirectory (OpalConfigure)
add_subdirectory (OpalParser)
add_subdirectory (Match)
add_subdirectory (PhysicsActions)
add_subdirectory (Solvers)
add_subdirectory (Steppers)
add_subdirectory (Structure)
add_subdirectory (Tables)
add_subdirectory (Track)
add_subdirectory (Utilities)
add_subdirectory (ValueDefinitions)

IF(ENABLE_AMR_SOLVER)
    add_subdirectory (Amr)
ENDIF(ENABLE_AMR_SOLVER)

ADD_OPAL_SOURCES(opal.cpp)
ADD_OPAL_SOURCES(changes.cpp)

include_directories ( BEFORE
    ${OPAL_SRC_DIR}
    ${CLASSIC_SRC_DIR}
    ${IPPL_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${GSL_INCLUDE_DIR}
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
    ${IPPL_SOURCE_DIRS}
    ${CCSE_INCLUDE_DIRS}
    ${GTEST_INCLUDE_DIR}
    ${OPTP_INCLUDE_DIR}
)

link_directories (
    ${IPPL_LIBRARY_DIR}
    ${OPTP_LIBRARY_DIR}
    ${Trilinos_LIBRARY_DIRS}
    ${Trilinos_TPL_LIBRARY_DIRS}
    ${CCSE_LIBRARY_DIR}
)

IF (ENABLE_DKS)
   include_directories (${DKS_INCLUDE_DIR})
   link_directories (${DKS_LIBRARY_DIR})
ENDIF (ENABLE_DKS)

IF(${HOSTNAME_BASE} MATCHES "edison" OR ${HOSTNAME_BASE} MATCHES "cori" OR ${HOSTNAME_BASE} MATCHES "daint")
 SET (OPAL_LIBS
    ${GSL_LIBRARY}
    ${H5Hut_LIBRARY}
    ${HDF5_LIBRARY}
    ${IPPL_LIBRARY}
    ${OPTP_LIBRARY}
    ${OPTP_LIBS}
    ${Boost_LIBRARIES}
    m
    z
)
ELSE()
 SET (OPAL_LIBS
    ${GSL_LIBRARY}
    ${H5Hut_LIBRARY}
    ${HDF5_LIBRARY}
    ${IPPL_LIBRARY}
    ${OPTP_LIBRARY}
    ${OPTP_LIBS}
    ${Boost_LIBRARIES}
    mpi
    mpi_cxx
    m
    z
)
ENDIF(${HOSTNAME_BASE} MATCHES "edison" OR ${HOSTNAME_BASE} MATCHES "cori" OR ${HOSTNAME_BASE} MATCHES "daint")

IF (ENABLE_DKS)
   set (OPAL_LIBS
       ${OPAL_LIBS}
       dks 
       ${DKS_CUDA_LIBS})
ENDIF (ENABLE_DKS)

add_library( OPALib ${OPAL_SRCS} )
set_target_properties( OPALib PROPERTIES OUTPUT_NAME OPAL )
target_link_libraries( OPALib ${OPAL_LIBS} ${Trilinos_LIBRARIES} ${CCSE_LIBRARIES} ${OTHER_CMAKE_EXE_LINKER_FLAGS} )

add_executable( opal Main.cpp )
target_link_libraries( opal ${OPAL_LIBS} OPALib ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} -lgfortran ${CCSE_LIBRARIES} ${CBOXLIB} ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${CMAKE_DL_LIBS})

install (TARGETS ${TEST_EXE} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
install (TARGETS opal RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
install (TARGETS OPALib DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

set (HDRS
    OPALconfig.h
    opal.h
    OPALrevision.h
    changes.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/src")
