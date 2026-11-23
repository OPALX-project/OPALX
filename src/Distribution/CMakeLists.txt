set (_SRCS
    Distribution.cpp
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

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Distribution")
