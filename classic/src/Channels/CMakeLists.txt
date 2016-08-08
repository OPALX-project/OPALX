set (_SRCS
  Channel.cpp
  ConstChannel.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_cl_sources(${_SRCS})