set (_SRCS
  )

set (_HDRS
  f2c.h
  FFTBase.cpp
  FFTBase.h
  FFT.cpp
  FFT.h
  fftpack_FFT.h
  SCSL_FFT.h
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
add_headers(${_HDRS})

install(FILES ${_HDRS} DESTINATION include/FFT)
