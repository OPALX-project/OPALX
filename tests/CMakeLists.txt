cmake_minimum_required(VERSION 2.6)

set (TEST_EXE opal_unit_tests)
set (SRC_DIRS classic_src opal_src)

set (CMAKE_CXX_FLAGS
    "${OPAL_CXX_FLAGS} -DWITH_UNIT_TESTS"
)

set (TEST_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
macro (add_sources)
    file (RELATIVE_PATH _relPath "${TEST_SRC_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND TEST_SRCS_LOCAL "${_relPath}/${_src}")
        else()
            list (APPEND TEST_SRCS_LOCAL "${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SRCS to parent directory
        set (TEST_SRCS_LOCAL ${TEST_SRCS_LOCAL} PARENT_SCOPE)
    endif()
endmacro()

add_subdirectory (opal_src)
add_subdirectory (classic_src)

set(TEST_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/Main.cpp ${TEST_SRCS_LOCAL})

LINK_DIRECTORIES( ${GTEST_LIB_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${IPPL_LIBRARY_DIR}
    ${CLASSIC_SOURCE_DIR}
    ${OPAL_SOURCE_DIR}/src
    ${Trilinos_LIBRARY_DIRS}
    ${Trilinos_TPL_LIBRARY_DIRS}
    ${CCSE_LIBRARY_DIR}
    ${Boost_LIBRARIES}
)

INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${IPPL_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${GSL_INCLUDE_DIR}
    ${CLASSIC_SOURCE_DIR}
    ${OPAL_SOURCE_DIR}/src
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
    ${IPPL_SOURCE_DIRS}
    ${CCSE_INCLUDE_DIRS}
)

# Check to see if cmake finds the test files
# MESSAGE(STATUS "unit test src files: ${TEST_SRCS}")

# Build the test exe. We don't do an install on the unit test exe as it is
# assumed that this is internal to opal
ADD_EXECUTABLE(${TEST_EXE} ${TEST_SRCS}) # the opal and classic sources are not needed again if we link agains libOPAL and libCLASSIC!
TARGET_LINK_LIBRARIES(${TEST_EXE} OPALib CLASSIC ${OPAL_LIBS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${CCSE_LIBRARIES} -lgfortran ${OTHER_CMAKE_EXE_LINKER_FLAGS} ${GTEST_BOTH_LIBRARIES} -lpthread)