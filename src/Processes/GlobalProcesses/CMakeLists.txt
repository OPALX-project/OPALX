set(_SRCS
    GlobalProcess.cpp
    Decay.cpp
    MuonDecay.cpp
    PionDecay.cpp
)

set(HDRS
    GlobalProcess.h
    Decay.h
    MuonDecay.h
    PionDecay.h
)

opalx_include_source_directory()
target_sources(opalx_run PRIVATE ${_SRCS})

install(FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Processes/GlobalProcesses")
