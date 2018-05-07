set (_SRCS
    SampleCmd.cpp
    SVar.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_opal_sources(${_SRCS})

set (HDRS
    Individual.h
    SampleCmd.h
    SVar.h
    Sampler.h
    Sampler.tcc
    SamplingOperator.h
    Uniform.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Sample")