CMAKE_MINIMUM_REQUIRED (VERSION 2.8.10)
PROJECT (StandAloneSDDSReader)
SET (StandAloneSDDSReader_VERSION_MAJOR 0)
SET (StandAloneSDDSReader_VERSION_MINOR 1)

configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

set (CMAKE_CXX_FLAGS
    "${IPPL_CMAKE_CXX_FLAGS} -DPARALLEL_IO ${CMAKE_CXX_FLAGS}"
    )

SET (PARSERPATH "${CMAKE_SOURCE_DIR}/src/Structure")

INCLUDE_DIRECTORIES (
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/Classic
)

SET (SRCS
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/array.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/associate.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/ast.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/column.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/data.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/description.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/file.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/include.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/parameter.cpp
  ${CMAKE_SOURCE_DIR}/src/Structure/SDDSParser/version.cpp
  ${CMAKE_SOURCE_DIR}/src/Utilities/OpalException.cpp
  ${CMAKE_SOURCE_DIR}/src/Classic/Utilities/ClassicException.cpp
)

SET (LIBS
    ${Boost_LIBRARIES}
    ${MPI_LIBRARIES}
    z
    dl
)

ADD_EXECUTABLE( SDDSReader main.cpp ${SRCS})
TARGET_LINK_LIBRARIES( SDDSReader ${LIBS} )

INSTALL(TARGETS SDDSReader RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")