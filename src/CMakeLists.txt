cmake_minimum_required (VERSION 2.6)

IF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.svn")
    include(FindSubversion)
    Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} OPAL)
ELSE(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.svn")
    set (OPAL_WC_REVISION 0)
    EXECUTE_PROCESS(
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMAND git rev-parse HEAD
            OUTPUT_VARIABLE GIT_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        set (OPAL_WC_REVISION \"${GIT_OUT}\")
ENDIF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.svn")

set (PACKAGE \"opal\")
set (PACKAGE_BUGREPORT \"opal@lists.psi.ch\")
set (PACKAGE_NAME \"OPAL\")
set (PACKAGE_STRING \"OPAL\ 1.2.0\")
set (PACKAGE_TARNAME \"opal\")
set (PACKAGE_VERSION \"1.2.0\")
set (VERSION \"1.2.0\")

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
add_subdirectory (Structure)
add_subdirectory (Tables)
add_subdirectory (Track)
add_subdirectory (Utilities)
add_subdirectory (ValueDefinitions)

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${IPPL_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${GSL_INCLUDE_DIR}
    ${CLASSIC_SOURCE_DIR}
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
    ${IPPL_INCLUDE_DIRS}
)

link_directories (
    ${IPPL_LIBRARY_DIR}
    ${CLASSIC_SOURCE_DIR}
    ${DOOM_SOURCE_DIR}
    ${Trilinos_LIBRARY_DIRS}
    ${Trilinos_TPL_LIBRARY_DIRS}
)
#   "/opt/parmetis/parmetis-3.1"
#   "/opt/intel-mkl/mkl-10.0/lib/em64t"

set (OPAL_LIBS
    m
    ${GSL_LIBRARY}
    ${H5Hut_LIBRARY_C}
    ${H5Hut_LIBRARY}
    ${HDF5_LIBRARY}
    ${IPPL_LIBRARY}
    z
)

add_executable( opal ${OPAL_SRCS} Main.cpp )
target_link_libraries( opal CLASSIC ${OPAL_LIBS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

if (BUILD_LIBOPAL)
    add_library( OPAL ${OPAL_SRCS} opal.cpp )
endif (BUILD_LIBOPAL)
