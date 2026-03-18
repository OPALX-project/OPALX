set (_SRCS
    Distribution.cpp
    FlatTop.cpp
    FromFile.cpp
    Gaussian.cpp
    LaserProfile.cpp
    MultiVariateGaussian.cpp
    SamplingBase.cpp
    )

set (HDRS
    Distribution.h
    FlatTop.h
    FromFile.h
    Gaussian.h
    LaserProfile.h
    MultiVariateGaussian.h
    SamplingBase.hpp
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Distribution")