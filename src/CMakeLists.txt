cmake_minimum_required (VERSION 2.8.10)

set (PACKAGE \"opal\")
set (PACKAGE_BUGREPORT \"opal@lists.psi.ch\")
set (PACKAGE_NAME \"OPAL\")
set (PACKAGE_TARNAME \"opal\")
set (PACKAGE_VERSION \"1.3.99\")

configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)


#FIXME: get from Makefile.def (IPPL)
set (CMAKE_CXX_FLAGS
    "${IPPL_CMAKE_CXX_FLAGS} -DUSEH5FEDV2 -DPARALLEL_IO ${CMAKE_CXX_FLAGS}"
    )

IF(ENABLE_SAAMG_SOLVER)
    set (CMAKE_CXX_FLAGS
        "-DHAVE_SAAMG_SOLVER ${CMAKE_CXX_FLAGS}"
    )
ENDIF(ENABLE_SAAMG_SOLVER)

IF(ENABLE_AMR_SOLVER)
    set (CMAKE_CXX_FLAGS
        "-DHAVE_AMR_SOLVER ${CMAKE_CXX_FLAGS}"
    )
ENDIF(ENABLE_AMR_SOLVER)

#MESSAGE ("Compile OPAL with: ${CMAKE_CXX_FLAGS}")

set (OPAL_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
macro (add_sources)
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

add_subdirectory (AbstractObjects)
add_subdirectory (Algebra)
add_subdirectory (Algorithms)
add_subdirectory (Aperture)
add_subdirectory (Attributes)
add_subdirectory (BasicActions)
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
add_subdirectory (Solvers/amr)
add_subdirectory (Structure)
add_subdirectory (Tables)
add_subdirectory (Track)
add_subdirectory (Utilities)
add_subdirectory (ValueDefinitions)

include_directories (Solvers/amr
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${IPPL_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${GSL_INCLUDE_DIR}
    ${CLASSIC_SOURCE_DIR}
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
    ${IPPL_SOURCE_DIRS}
    ${CCSE_INCLUDE_DIRS}
    ${GTEST_INCLUDE_DIR}
)

link_directories (
    ${IPPL_LIBRARY_DIR}
    ${CLASSIC_SOURCE_DIR}
    ${Trilinos_LIBRARY_DIRS}
    ${Trilinos_TPL_LIBRARY_DIRS}
    ${CCSE_LIBRARY_DIR}
)

set (OPAL_LIBS
    m
    ${GSL_LIBRARY}
    ${H5Hut_LIBRARY}
    ${HDF5_LIBRARY}
    ${IPPL_LIBRARY}
    z
    mpi
    mpi_cxx
    ${Boost_LIBRARIES}
)

add_library( OPALib ${OPAL_SRCS} opal.cpp )
set_target_properties( OPALib PROPERTIES OUTPUT_NAME OPAL )

add_executable( opal Main.cpp )
target_link_libraries( opal OPALib CLASSIC ${OPAL_LIBS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${CCSE_LIBRARIES} ${OTHER_CMAKE_EXE_LINKER_FLAGS} -lgfortran ${CMAKE_DL_LIBS})

install (TARGETS ${TEST_EXE} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
install (TARGETS opal RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
install (TARGETS OPALib DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
install (FILES opal.h DESTINATION "${CMAKE_INSTALL_PREFIX}/include")
install (FILES Utilities/OpalException.h DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Utilities")

set (OPAL_LIBS ${OPAL_LIBS} PARENT_SCOPE)
set (OPAL_CXX_FLAGS ${CMAKE_CXX_FLAGS} PARENT_SCOPE)