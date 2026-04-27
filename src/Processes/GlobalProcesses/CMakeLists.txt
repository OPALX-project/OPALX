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

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(opalx PRIVATE ${_SRCS})

install(FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Processes/GlobalProcesses")
