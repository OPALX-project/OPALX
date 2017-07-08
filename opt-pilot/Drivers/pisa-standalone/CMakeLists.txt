
# this prevents compile errors (due to warnings) since Boost Spirits grammar
# definitions are interpreted as parentheses errors (overloaded operators).
set_source_files_properties(${CMAKE_SOURCE_DIR}/Expression/Parser/expression.cpp PROPERTIES COMPILE_FLAGS "-Wno-error=parentheses" )

include_directories (
    ${CMAKE_SOURCE_DIR}
    )
link_directories (
    ${Boost_LIBRARY_DIRS}
    )
add_executable( pisa-standalone.exe
    pisa-standalone.cpp
    ${_SRCS} ${PILOT_FEVAL_SRCS}
    )
target_link_libraries ( pisa-standalone.exe
    ${PILOT_LIBS}
    ${BOOST_FILESYSTEM}
    ${MPI_CXX_LIBRARIES}
    -lgfortran
    )

if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_link_libraries( pisa-standalone.exe
	"-Wl,--allow-multiple-definition")
endif ()

set_target_properties( pisa-standalone.exe PROPERTIES
    COMPILE_FLAGS "-Wall -std=c++0x"
    )

install (TARGETS pisa-standalone.exe
    DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
    )
