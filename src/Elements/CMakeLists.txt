set (_SRCS
    OpalCavity.cpp
    OpalTravelingWave.cpp
    OpalDrift.cpp
    OpalElement.cpp
    OpalMarker.cpp
    OpalMultipole.cpp
    OpalMultipoleT.cpp
    OpalSolenoid.cpp
    OpalBeamline.cpp
    OpalRingDefinition.cpp      
    OpalOffset/OpalLocalCylindricalOffset.cpp
    OpalOffset/OpalLocalCartesianOffset.cpp
    OpalOffset/OpalGlobalCylindricalOffset.cpp
    OpalOffset/OpalGlobalCartesianOffset.cpp
    OpalVerticalFFAMagnet.cpp
    OpalProbe.cpp
)

set (HDRS
    OpalBeamline.h
    OpalCavity.h
    OpalDrift.h
    OpalElement.h
    OpalMarker.h
    OpalMultipole.h
    OpalMultipoleT.h
    OpalSource.h
    OpalSolenoid.h
    OpalTravelingWave.h
    OpalRingDefinition.h
    OpalOffset/OpalGlobalCartesianOffset.h
    OpalOffset/OpalGlobalCylindricalOffset.h
    OpalOffset/OpalLocalCartesianOffset.h
    OpalOffset/OpalLocalCylindricalOffset.h
    OpalVerticalFFAMagnet.h
    OpalProbe.h
)

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Elements")