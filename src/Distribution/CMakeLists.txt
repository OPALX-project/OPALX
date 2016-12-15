set (_SRCS
  Distribution.cpp
  halton1d_sequence.cpp
  LaserProfile.cpp
  MagneticField.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
