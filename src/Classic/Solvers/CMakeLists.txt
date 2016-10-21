set (_SRCS
  CollimatorPhysics.cpp
  CSRWakeFunction.cpp
  CSRIGFWakeFunction.cpp
  GreenWakeFunction.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_sources(${_SRCS})
