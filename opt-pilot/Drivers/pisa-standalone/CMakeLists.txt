include_directories (
    ${CMAKE_SOURCE_DIR}
    )
link_directories (
    ${OPTP_LIBRARY_DIR}
    ${Boost_LIBRARY_DIRS}
    )

add_executable( pisa-standalone.exe
    pisa-standalone.cpp
    ${PILOT_FEVAL_SRCS}
    )
target_link_libraries ( pisa-standalone.exe
    ${OPTP_LIBS}
    ${OPTP_LIBRARY}
    ${MPI_CXX_LIBRARIES}
    )

message(STATUS "OPTP_LIBRARY ${OPTP_LIBRARY}")
message(STATUS "OPTP_LIBRARY_DIR ${OPTP_LIBRARY_DIR}")

if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_link_libraries( pisa-standalone.exe
	"-Wl,--allow-multiple-definition")
endif ()

set_target_properties( pisa-standalone.exe PROPERTIES
    COMPILE_FLAGS "-Wall -std=c++11"
    )

install (TARGETS pisa-standalone.exe
    DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
    )