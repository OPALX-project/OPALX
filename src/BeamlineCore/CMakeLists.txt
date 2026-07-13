set (_SRCS
    ConstantEFieldCavityRep.cpp
    DriftRep.cpp
    LaserRep.cpp
    MarkerRep.cpp
    MonitorRep.cpp
    MultipoleRep.cpp
    RBendRep.cpp
    RFCavityRep.cpp
    SBendRep.cpp
    TravelingWaveRep.cpp
    SolenoidRep.cpp
    ProbeRep.cpp
    )

set (HDRS
    ConstantEFieldCavityRep.h
    DriftRep.h
    LaserRep.h
    MarkerRep.h
    MonitorRep.h
    MultipoleRep.h
    RBendRep.h
    RFCavityRep.h
    SBendRep.h
    TravelingWaveRep.h
    SolenoidRep.h
    ProbeRep.h
)

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/BeamlineCore")
