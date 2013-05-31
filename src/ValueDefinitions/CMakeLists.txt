set (_SRCS
  BoolConstant.cpp
  RealConstant.cpp
  RealVariable.cpp
  RealVector.cpp
  StringConstant.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
