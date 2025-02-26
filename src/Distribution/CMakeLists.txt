set (_SRCS
    Distribution.cpp
    FlatTop.cpp
    Gaussian.cpp
    LaserProfile.cpp
    MultiVariateGaussian.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    Distribution.h
    FlatTop.h
    Gaussian.h
    LaserProfile.h
    MapGenerator.h
    MultiVariateGaussian.h
    matrix_vector_operation.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Distribution")
