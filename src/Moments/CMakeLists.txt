set (_SRCS
    BeamMoments.cpp
    MomentComputation.cpp
)

set (HDRS
    BeamMoments.h
    MomentRequest.h
)

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Moments")
