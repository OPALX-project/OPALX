set (_SRCS
    LossDataSink.cpp
    MeshGenerator.cpp
    PeakFinder.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    LossDataSink.h
    MeshGenerator.h
    PeakFinder.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Structure")