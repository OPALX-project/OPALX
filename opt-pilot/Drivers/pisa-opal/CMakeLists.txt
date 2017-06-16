# this prevents compile errors (due to warnings) since Boost Spirits grammar
# definitions are interpreted as parentheses errors (overloaded operators).
set_source_files_properties(${CMAKE_SOURCE_DIR}/Expression/Parser/expression.cpp PROPERTIES COMPILE_FLAGS "-Wno-error=parentheses" )

if (Boost_FILESYSTEM_FOUND)
    set (BOOST_FILESYSTEM
        boost_filesystem
	)
else (Boost_FILESYSTEM_FOUND)
    set (BOOST_FILESYSTEM
	)
endif (Boost_FILESYSTEM_FOUND)

if (NOT ${GSL_FOUND})
    message (FATAL "GSL not found")
endif ()

if (NOT ${HDF5_FOUND})
    message (FATAL "HDF5 not found")
endif ()

if (NOT ${H5HUT_FOUND})
    message (FATAL "H5hut not found")
endif ()

if (NOT ${IPPL_FOUND})
    message (FATAL "IPPL not found")
endif ()

if (NOT ${OPAL_FOUND})
    message (FATAL "OPAL not found")
endif ()

set (OPAL_LIBS
    ${OPAL_LIBRARY}
    ${IPPL_LIBRARY}
    ${GSL_LIBRARY}
    ${H5HUT_LIBRARY}
    ${HDF5_LIBRARY}
    )

IF (ENABLE_DKS)
   set (OPAL_LIBS
       ${OPAL_LIBS}
       ${DKS_LIBRARY_DIR}/libdks.a
       )
ENDIF (ENABLE_DKS)

include_directories (
    ${CMAKE_SOURCE_DIR}
    ${H5HUT_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${GSL_INCLUDE_DIR}
    ${IPPL_INCLUDE_DIR}
    ${OPAL_INCLUDE_DIR}
    )
link_directories (
    ${IPPL_LIBRARY_DIR}
    ${OPAL_LIB_DIR}
    ${Boost_LIBRARY_DIRS}
    )
add_executable( pisa-opal.exe
    pisa-opal.cpp ${_SRCS} ${PILOT_OPAL_SRCS}
    )

target_link_libraries( pisa-opal.exe
    ${OPAL_LIBS}
    ${MPI_CXX_LIBRARIES}
    ${Trilinos_LIBRARIES}
    -lgfortran
    ${PILOT_LIBS}
    ${Boost_LIBRARIES}
    z
    pthread
    m
    )

if (NOT  ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_link_libraries( pisa-opal.exe
	"-Wl,--allow-multiple-definition"
	)
endif ()
set_target_properties( pisa-opal.exe PROPERTIES
    COMPILE_FLAGS "${IPPL_CMAKE_CXX_FLAGS} -Wall -std=c++0x")

install (TARGETS pisa-opal.exe
    DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
    )
