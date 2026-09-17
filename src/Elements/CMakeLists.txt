set (_SRCS
    OpalAsymmetricEnge.cpp
    OpalCavity.cpp
    OpalCollimator.cpp
    OpalConstantEFieldCavity.cpp
    OpalConstantFocusing.cpp
    OpalTravelingWave.cpp
    OpalBend.cpp
    OpalDrift.cpp
    OpalElement.cpp
    OpalEnge.cpp
    OpalLaser.cpp
    OpalMonitor.cpp
    OpalMarker.cpp
    OpalMultipole.cpp
    OpalMultipoleT.cpp
    OpalQuadrupole.cpp
    OpalRBend.cpp
    OpalSBend.cpp
    OpalScalingFFAMagnet.cpp
    OpalSolenoid.cpp
    OpalBeamline.cpp
    PlacementResolver.cpp
    OpalVerticalFFAMagnet.cpp
    OpalProbe.cpp
    OpalPolynomialTimeDependence.cpp
    OpalSplineTimeDependence.cpp
    OpalSinusoidalTimeDependence.cpp
    OpalTanh.cpp
    OpalVariableRFCavity.cpp
)

set (HDRS
    OpalAsymmetricEnge.h
    OpalBeamline.h
    PlacementResolver.h
    OpalBend.h
    OpalCavity.h
    OpalCollimator.h
    OpalConstantEFieldCavity.h
    OpalConstantFocusing.h
    OpalDrift.h
    OpalElement.h
    OpalEnge.h
    OpalMonitor.h
    OpalMarker.h
    OpalMultipole.h
    OpalMultipoleT.h
    OpalQuadrupole.h
    OpalRBend.h
    OpalSBend.h
    OpalScalingFFAMagnet.h
    OpalSource.h
    OpalSolenoid.h
    OpalTravelingWave.h
    OpalLaser.h
    OpalVerticalFFAMagnet.h
    OpalProbe.h
    OpalPolynomialTimeDependence.h
    OpalSplineTimeDependence.h
    OpalSinusoidalTimeDependence.h
    OpalTanh.h
    OpalVariableRFCavity.h
)

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Elements")
