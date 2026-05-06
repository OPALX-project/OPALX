set (_SRCS
    OpalCavity.cpp
    OpalConstantEFieldCavity.cpp
    OpalTravelingWave.cpp
    OpalBend.cpp
    OpalDrift.cpp
    OpalElement.cpp
    OpalLaser.cpp
    OpalMarker.cpp
    OpalMultipole.cpp
    OpalMultipoleT.cpp
    OpalQuadrupole.cpp
    OpalRBend.cpp
    OpalSBend.cpp
    OpalSolenoid.cpp
    OpalBeamline.cpp
    OpalRingDefinition.cpp      
    OpalVerticalFFAMagnet.cpp
    OpalProbe.cpp
    OpalPolynomialTimeDependence.cpp
    OpalSplineTimeDependence.cpp
    OpalSinusoidalTimeDependence.cpp
    OpalVariableRFCavity.cpp
)

set (HDRS
    OpalBeamline.h
    OpalBend.h
    OpalCavity.h
    OpalConstantEFieldCavity.h
    OpalDrift.h
    OpalElement.h
    OpalMarker.h
    OpalMultipole.h
    OpalMultipoleT.h
    OpalQuadrupole.h
    OpalRBend.h
    OpalSBend.h
    OpalSource.h
    OpalSolenoid.h
    OpalTravelingWave.h
    OpalRingDefinition.h
    OpalLaser.h
    OpalVerticalFFAMagnet.h
    OpalProbe.h
    OpalPolynomialTimeDependence.h
    OpalSplineTimeDependence.h
    OpalSinusoidalTimeDependence.h
    OpalVariableRFCavity.h
)

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Elements")
