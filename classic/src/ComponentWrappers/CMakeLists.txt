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

add_cl_sources(${_SRCS})