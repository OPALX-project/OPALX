set (_SRCS
    Distribution.cpp
    FlatTop.cpp
    Gaussian.cpp
    LaserProfile.cpp
    MultiVariateGaussian.cpp
    )

set (HDRS
    Distribution.h
    FlatTop.h
    Gaussian.h
    LaserProfile.h
    MapGenerator.h
    MultiVariateGaussian.h
    matrix_vector_operation.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Distribution")