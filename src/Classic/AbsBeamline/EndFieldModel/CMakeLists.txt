set (_SRCS
  Tanh.cpp 
  Enge.cpp
  Sin.cpp
  EndFieldModel.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
  Tanh.h
  Enge.h
  Sin.h
  EndFieldModel.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/AbsBeamline/EndFieldModel/")

