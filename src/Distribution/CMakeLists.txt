set (_SRCS
  Distribution.cpp
  LaserProfile.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_opal_sources(${_SRCS})

set (HDRS
    ClosedOrbitFinder.h
    Distribution.h
    Harmonics.h
    LaserProfile.h
    MapGenerator.h
    matrix_vector_operation.h
    rdm.h
    SigmaGenerator.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Distribution")