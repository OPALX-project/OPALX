set (_SRCS
    CoordinateSystemTrafo.cpp
    Quaternion.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    CoordinateSystemTrafo.h
    Quaternion.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")
