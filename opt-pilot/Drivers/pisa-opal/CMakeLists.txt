if (Boost_FILESYSTEM_FOUND)
    set (BOOST_FILESYSTEM
        boost_filesystem
	)
else (Boost_FILESYSTEM_FOUND)
    set (BOOST_FILESYSTEM
	)
endif (Boost_FILESYSTEM_FOUND)

# OPAL
message(STATUS "OPAL_LIBS ${OPAL_LIBS}")

set (OPAL_LIBS
    ${OPAL_LIBS}
    ${OPTP_LIBRARY}
    ${IPPL_LIBRARY}
    ${GSL_LIBRARY}
    ${H5HUT_LIBRARY}
    ${HDF5_LIBRARY}
    )

include_directories (
    ${CMAKE_SOURCE_DIR}
    ${H5HUT_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${GSL_INCLUDE_DIR}
    ${IPPL_INCLUDE_DIR}
    )
link_directories (
    ${OPTP_LIBRARY_DIR}
    ${IPPL_LIBRARY_DIR}
    ${Boost_LIBRARY_DIRS}
    )
add_executable( pisa-opal.exe
    pisa-opal.cpp ${_SRCS} ${PILOT_OPAL_SRCS}
    )

target_link_libraries( pisa-opal.exe
    ${OPAL_LIBS}
    OPALib
    ${MPI_CXX_LIBRARIES}
    ${Trilinos_LIBRARIES}
    -lgfortran
    ${OPTP_LIBS}
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
