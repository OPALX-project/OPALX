add_definitions (-DPARALLEL_IO)

get_directory_property (OPAL_COMPILE_OPTION_LIST COMPILE_OPTIONS)
string (REGEX REPLACE "([^\\]|^);" "\\1 " OPAL_COMPILE_OPTIONS "${OPAL_COMPILE_OPTION_LIST}")
configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/OPALconfig.h)

set (OPAL_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set (OPAL_SRCS )
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

macro (ADD_OPAL_HEADERS )
    file (RELATIVE_PATH _relPath "${OPAL_SRC_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_hdr ${ARGN})
        if (_relPath)
            list (APPEND OPAL_HDRS "${_relPath}/${_hdr}")
        else ()
            list (APPEND OPAL_HDRS "${_hdr}")
        endif ()
    endforeach ()
    if (_relPath)
        # propagate HDRS to parent directory
        set (OPAL_HDRS ${OPAL_HDRS} PARENT_SCOPE)
    endif ()
endmacro ()

add_opal_sources (opal.cpp)
add_opal_sources (changes.cpp)
add_subdirectory (AbstractObjects)
add_subdirectory (Algorithms)
add_subdirectory (Attributes)
add_subdirectory (BasicActions)
add_subdirectory (Classic)
add_subdirectory (Distribution)
add_subdirectory (Elements)
add_subdirectory (Expressions)
add_subdirectory (Lines)
add_subdirectory (OpalConfigure)
add_subdirectory (OpalParser)
add_subdirectory (Physics)
add_subdirectory (Solvers)
add_subdirectory (Steppers)
add_subdirectory (Structure)
add_subdirectory (Tables)
add_subdirectory (Track)
add_subdirectory (Utilities)
add_subdirectory (ValueDefinitions)

include_directories (BEFORE
    ${OPAL_SRC_DIR}
    ${CLASSIC_SRC_DIR}
    ${IPPL_INCLUDE_DIR}
    ${H5Hut_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${GSL_INCLUDE_DIR}
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
    ${IPPL_SOURCE_DIRS}
    ${GTEST_INCLUDE_DIRS}
    )

link_directories (
    ${IPPL_LIBRARY_DIR}
    ${Trilinos_LIBRARY_DIRS}
    ${Trilinos_TPL_LIBRARY_DIRS}
    )

#
# we cannot use the target name 'OPAL' for 'libOPAL'. On systems with
# case-insensitive filesystems (like macOS) this conflicts with the
# target 'opal' for the OPAL binary.
#
add_library (libOPAL ${OPAL_SRCS})
set_target_properties (libOPAL PROPERTIES OUTPUT_NAME OPAL)

target_link_libraries (libOPAL
    ${MPI_CXX_LIBRARIES}
    )

add_executable (opal Main.cpp)
target_link_libraries( opal
    libOPAL
    ${IPPL_LIBRARY}
    ${GSL_LIBRARY}
    ${GSL_CBLAS_LIBRARY}
    ${H5Hut_LIBRARY}
    ${HDF5_LIBRARIES}
    ${Boost_LIBRARIES}
    m
    z
    ${CMAKE_DL_LIBS}
    )

if (ENABLE_OPAL_FEL)
    include_directories (${MITHRA_INCLUDE_DIR})
    target_link_libraries(libOPAL ${MITHRA_LIBRARY_DIR}/libmithra.a)
endif()

install (TARGETS ${TEST_EXE} RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
install (TARGETS opal RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
install (TARGETS libOPAL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

set (HDRS
    OPALconfig.h
    opal.h
    OPALrevision.h
    changes.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/src")