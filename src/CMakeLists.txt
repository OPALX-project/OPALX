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
    ${IPPL_INCLUDE_DIRS}
    ${CCSE_INCLUDE_DIRS}
)

link_directories (
    ${IPPL_LIBRARY_DIR}
    ${CLASSIC_SOURCE_DIR}
    ${DOOM_SOURCE_DIR}
    ${Trilinos_LIBRARY_DIRS}
    ${Trilinos_TPL_LIBRARY_DIRS}
    ${CCSE_LIBRARY_DIR}
)
#   "/opt/parmetis/parmetis-3.1"
#   "/opt/intel-mkl/mkl-10.0/lib/em64t"

set (OPAL_LIBS
    m
    ${GSL_LIBRARY}
    ${H5Hut_LIBRARY}
    ${HDF5_LIBRARY}
    ${IPPL_LIBRARY}
    z
    mpi
    mpi_cxx
)

add_library( OPALib ${OPAL_SRCS} opal.cpp )
set_target_properties( OPALib PROPERTIES OUTPUT_NAME OPAL )

add_executable( opal Main.cpp )
target_link_libraries( opal OPALib CLASSIC ${OPAL_LIBS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${CCSE_LIBRARIES} ${OTHER_CMAKE_EXE_LINKER_FLAGS} -lgfortran ${CMAKE_DL_LIBS})

# build unit tests; unit tests are defined in unit_tests directory, we build the
# executable here so that we can use the OPAL src and lib definitions
if (BUILD_OPAL_UNIT_TESTS)
    # touch CMakeLists to force cmake to redo the glob for new test source files
    # otherwise cmake does not detect new source files
    execute_process(COMMAND touch ${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/CMakeLists.txt)
    add_subdirectory(unit_tests)
    # include and link gtest library
    # this can be built using the script provided at
    #     unit_test/gtest/bin/install_gtest.bash
    link_directories( ${GTEST_LIB_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../classic/5.0/src)
    include_directories(${GTEST_INCLUDE_DIR})
    # Check to see if cmake finds the test files
    message(STATUS "unit test src files: ${TEST_SRCS}")
    # Build the test exe. We don't do an install on the unit test exe as it is
    # assumed that this is internal to opal
    add_executable(${TEST_EXE} ${TEST_SRCS}) # the opal and classic sources are not needed again if we link agains libOPAL and libCLASSIC!
    target_link_libraries(${TEST_EXE} OPALib CLASSIC ${OPAL_LIBS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${CCSE_LIBRARIES} -lgfortran ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${GTEST_LINK_TARGETS} )
endif (BUILD_OPAL_UNIT_TESTS)

install (TARGETS ${TEST_EXE} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
install (TARGETS opal RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
install (TARGETS OPALib DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
install (FILES opal.h DESTINATION "${CMAKE_INSTALL_PREFIX}/include")
install (FILES Utilities/OpalException.h DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Utilities")