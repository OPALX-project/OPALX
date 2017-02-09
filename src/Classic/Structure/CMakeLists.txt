set (_SRCS
    LossDataSink.cpp
    MeshGenerator.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    LossDataSink.h
    MeshGenerator.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Structure")