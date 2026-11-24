set (_SRCS
    BeamMoments.cpp
    MomentComputation.cpp
)

set (HDRS
    BeamMoments.h
    MomentRequest.h
)

opalx_include_source_directory()

target_sources(opalx_run PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Moments")
