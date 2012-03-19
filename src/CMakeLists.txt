cmake_minimum_required (VERSION 2.6)

IF(EXISTS ".svn")
    include(FindSubversion)
    Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} IPPL)
ELSE(EXISTS ".svn")
    set (IPPL_WC_REVISION 0)
    EXECUTE_PROCESS(
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMAND git rev-parse HEAD 
            OUTPUT_VARIABLE GIT_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        set (IPPL_WC_REVISION \"${GIT_OUT}\")
ENDIF(EXISTS ".svn")

configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

MESSAGE("configure: make IpplInfo.h")
EXECUTE_PROCESS(
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMAND date
        OUTPUT_VARIABLE DATE_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
    )
EXECUTE_PROCESS(
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMAND uname -a
        OUTPUT_VARIABLE UNAME_OUT OUTPUT_STRIP_TRAILING_WHITESPACE
    )

set (IPPL_COMPILE_ARCH \"$ENV{IPPL_ARCH}\")
set (IPPL_COMPILE_LINE \"${CMAKE_CXX_FLAGS}\")
set (IPPL_COMPILE_DATE \"${DATE_OUT}\")
set (IPPL_COMPILE_MACHINE \"${UNAME_OUT}\")
set (IPPL_COMPILE_OPTIONS \"${CXXDEFINES}\")
set (IPPL_COMPILE_USER \"adelmann\")
#FIXME: "CVS" version of each file:
#echo `find ${IPPL_ROOT}/src -name '*.cpp' -print` `find ${IPPL_ROOT}/src -name '*.h' -print` | perl ${IPPL_ROOT}/bin/mkver2 >> IpplVersions.h
configure_file(IpplVersions.h.in ${CMAKE_CURRENT_SOURCE_DIR}/IpplVersions.h)

MESSAGE (${CMAKE_CXX_FLAGS})

set (IPPL_SRCS_FORT
	FFT/fftpack.F
)

set (IPPL_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
macro (add_sources)
    file (RELATIVE_PATH _relPath "${IPPL_SRC_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND IPPL_SRCS "${_relPath}/${_src}")
        else()
            list (APPEND IPPL_SRCS "${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SRCS to parent directory
        set (IPPL_SRCS ${IPPL_SRCS} PARENT_SCOPE)
    endif()
endmacro()

macro (add_headers)
    file (RELATIVE_PATH _relPath "${IPPL_SRC_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_hdr ${ARGN})
        if (_relPath)
            list (APPEND IPPL_HDRS "${_relPath}/${_hdr}")
        else()
            list (APPEND IPPL_HDRS "${_hdr}")
        endif()
    endforeach()
    if (_relPath)
        # propagate HDRS to parent directory
        set (IPPL_HDRS ${IPPL_HDRS} PARENT_SCOPE)
    endif()
endmacro()

set(IPPL_BASEDIR_HDRS
  Ippl.h
  IpplVersions.h
  )
add_headers(${IPPL_BASEDIR_HDRS})

add_subdirectory(AppTypes)
add_subdirectory(DataSource)
add_subdirectory(Field)
add_subdirectory(FieldLayout)
add_subdirectory(Index)
add_subdirectory(Meshes)
add_subdirectory(Message)
add_subdirectory(Profile)
add_subdirectory(Utility)
add_subdirectory(PETE)

set (MAKE_DEF
    ../test/Makefile.def
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library( ippl ${IPPL_SRCS} ${IPPL_SRCS_FORT} )

install(FILES ${MAKE_DEF} DESTINATION lib)
install(TARGETS ippl DESTINATION lib)
install(FILES ${IPPL_BASEDIR_HDRS} DESTINATION include)
