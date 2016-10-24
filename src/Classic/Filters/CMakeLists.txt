set (_SRCS
  SavitzkyGolay.cpp
  Stencil.cpp
  FixedFFTLowPass.cpp
  RelativeFFTLowPass.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_sources(${_SRCS})