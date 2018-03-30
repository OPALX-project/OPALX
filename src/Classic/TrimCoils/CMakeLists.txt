set (_SRCS
    OpalTrimCoil.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    OpalTrimCoil.h
    TrimCoil.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/TrimCoils")
install (FILES ${_SRCS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/TrimCoils")
