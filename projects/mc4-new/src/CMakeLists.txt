cmake_minimum_required (VERSION 2.6)

IF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.svn")
    include(FindSubversion)
    Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} MC4)
ELSE(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.svn")
    set (MC4_WC_REVISION 0)
    EXECUTE_PROCESS(
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMAND git rev-parse HEAD
            OUTPUT_VARIABLE GIT_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        set (MC4_WC_REVISION \"${GIT_OUT}\")
ENDIF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.svn")

set (PACKAGE \"mc4\")
set (PACKAGE_BUGREPORT \"andreas.adelmann@psi.ch\")
set (PACKAGE_NAME \"MC4\")
set (PACKAGE_STRING \"MC4\ 0.9\")
set (PACKAGE_TARNAME \"mc4\")
set (PACKAGE_VERSION \"0.9\")
set (VERSION \"0.9\")

configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)


#FIXME: get from Makefile.def (IPPL)
set (CMAKE_CXX_FLAGS
    "${IPPL_CMAKE_CXX_FLAGS} -DUSEH5FEDV2 -DPARALLEL_IO ${CMAKE_CXX_FLAGS}"
    )

IF(ENABLE_ML_SOLVER)
    set (CMAKE_CXX_FLAGS
        "-DHAVE_ML_SOLVER ${CMAKE_CXX_FLAGS}"
    )
ENDIF(ENABLE_ML_SOLVER)

#MESSAGE ("Compile MC4 with: ${CMAKE_CXX_FLAGS}")

set (MC4_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
macro (add_sources)
    file (RELATIVE_PATH _relPath "${MC4_SRC_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND MC4_SRCS "${_relPath}/${_src}")
        else()
            list (APPEND MC4_SRCS "${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SRCS to parent directory
        set (MC4_SRCS ${MC4_SRCS} PARENT_SCOPE)
    endif()
endmacro()

add_subdirectory (Parser)
add_subdirectory (Utilities)
add_subdirectory (Initializer)
add_subdirectory (PwrSpec)
add_subdirectory (ChargedParticles)
add_subdirectory (FieldSolvers)
add_subdirectory (DataSink)

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${IPPL_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${GSL_INCLUDE_DIR}
    ${IPPL_INCLUDE_DIRS}
)

link_directories (
    ${IPPL_LIBRARY_DIR}
)

set (MC4_LIBS
    m
    ${GSL_LIBRARY}
    ${H5Hut_LIBRARY_C}
    ${H5Hut_LIBRARY}
    ${HDF5_LIBRARY}
    ${IPPL_LIBRARY}
    z
)

add_executable( mc4 ${MC4_SRCS} mc4.cpp )
target_link_libraries( mc4  ${MC4_LIBS} )

