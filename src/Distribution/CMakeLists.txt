# Input side: the DISTRIBUTION command.
set (_INPUT_SRCS
    Distribution.cpp
    )

# Run side: the values the samplers read, and the samplers.
set (_RUN_SRCS
    DistributionSettings.cpp
    EmittedFromFile.cpp
    EmittingSamplers.cpp
    FlatTop.cpp
    FromFile.cpp
    Gaussian.cpp
    MultiVariateGaussian.cpp
    OpalFlatTop.cpp
    SamplingBase.cpp
    Uniform.cpp
    )

set (HDRS
    Distribution.h
    DistributionSettings.h
    EmittedFromFile.h
    EmittingSamplers.h
    FlatTop.h
    FromFile.h
    Gaussian.h
    MultiVariateGaussian.h
    OpalFlatTop.h
    SamplingBase.hpp
    Uniform.h
    )

opalx_include_source_directory()

target_sources(opalx PRIVATE ${_INPUT_SRCS})
target_sources(opalx_run PRIVATE ${_RUN_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Distribution")
