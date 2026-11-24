set (_SRCS
    LinearCompton.cpp
    LinearBreitWheeler.cpp
    ParticleProperties.cpp
    )

set (HDRS
    LinearCompton.h
    LinearBreitWheeler.h
    MuonDecay.h
    ParticleProperties.h
    Physics.h
    )

opalx_include_source_directory()

target_sources(opalx_run PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Physics")