set (_SRCS
  Distribution.cpp
  LaserProfile.cpp
  MagneticField.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_opal_sources(${_SRCS})

set (HDRS
    ClosedOrbitFinder.h
    Distribution.h
    error.h
    Harmonics.h
    LaserProfile.h
    MagneticField.h
    MapGenerator.h
    matrix_vector_operation.h
    physical_error.h
    rdm.h
    SigmaGenerator.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Distribution")