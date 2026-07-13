set (_SRCS
    Distribution.cpp
    EmittedFromFile.cpp
    FlatTop.cpp
    FromFile.cpp
    Gaussian.cpp
    LaserProfile.cpp
    MultiVariateGaussian.cpp
    OpalFlatTop.cpp
    SamplingBase.cpp
    )

set (HDRS
    Distribution.h
    EmittedFromFile.h
    FlatTop.h
    FromFile.h
    Gaussian.h
    LaserProfile.h
    MultiVariateGaussian.h
    OpalFlatTop.h
    SamplingBase.hpp
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Distribution")
