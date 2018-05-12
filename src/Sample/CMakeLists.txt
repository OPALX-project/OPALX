set (_SRCS
    OpalSample.cpp
    SampleCmd.cpp
    Sampler.cpp
    SequenceSampleCmd.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_opal_sources(${_SRCS})

set (HDRS
    FromFile.h
    OpalSample.h
    SampleCmd.h
    SequenceSampleCmd.h
    SampleIndividual.h
    SamplePilot.h
    Sampler.h
    SampleSequence.h
    SampleWorker.h
    SamplingOperator.h
    Uniform.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Sample")