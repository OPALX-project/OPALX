set (_SRCS
    CoordinateSystemTrafo.cpp
    OpalParticle.cpp
    Quaternion.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    CoordinateSystemTrafo.h
    ListElem.h
    OpalParticle.h
    Quaternion.h
    Vektor.h
    Matrix.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")
