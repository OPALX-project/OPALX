set (_SRCS
  CorrectorWrapper.cpp
  MultipoleWrapper.cpp
  SBendWrapper.cpp
  RBendWrapper.cpp
  CyclotronWrapper.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    CorrectorWrapper.h
    CyclotronWrapper.h
    MultipoleWrapper.h
    RBendWrapper.h
    SBendWrapper.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ComponentWrappers")