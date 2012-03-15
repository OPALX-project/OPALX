set (_SRCS
  Distribution.cpp
  ranlib.cc
  halton1d_sequence.cpp
  LaserProfile.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
