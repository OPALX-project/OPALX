set (_SRCS
    DriftRep.cpp
    MarkerRep.cpp
    MonitorRep.cpp
    MultipoleRep.cpp
    RFCavityRep.cpp
    TravelingWaveRep.cpp
    SolenoidRep.cpp
    ProbeRep.cpp
    )

set (HDRS
    DriftRep.h
    MarkerRep.h
    MonitorRep.h
    MultipoleRep.h
    RFCavityRep.h
    TravelingWaveRep.h
    SolenoidRep.h
    ProbeRep.h
)

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/BeamlineCore")