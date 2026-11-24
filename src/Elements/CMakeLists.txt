# Input side: the element and time-dependence commands.
set (_INPUT_SRCS
    OpalCavity.cpp
    OpalBox.cpp
    OpalCollimator.cpp
    OpalConstantEFieldCavity.cpp
    OpalConstantFocusing.cpp
    OpalTravelingWave.cpp
    OpalBend.cpp
    OpalDrift.cpp
    OpalFieldmapElement.cpp
    OpalElement.cpp
    OpalLaser.cpp
    OpalMonitor.cpp
    OpalMarker.cpp
    OpalMultipole.cpp
    OpalMultipoleT.cpp
    OpalQuadrupole.cpp
    OpalRBend.cpp
    OpalSBend.cpp
    OpalSolenoid.cpp
    OpalPolynomialTimeDependence.cpp
    OpalSplineTimeDependence.cpp
    OpalSinusoidalTimeDependence.cpp
    OpalVariableRFCavity.cpp
)

# Run side: the runtime beamline and the placement of its elements.
set (_RUN_SRCS
    OpalBeamline.cpp
    PlacementResolver.cpp
)

set (HDRS
    OpalBeamline.h
    PlacementResolver.h
    OpalBend.h
    OpalCavity.h
    OpalBox.h
    OpalCollimator.h
    OpalConstantEFieldCavity.h
    OpalConstantFocusing.h
    OpalDrift.h
    OpalFieldmapElement.h
    OpalElement.h
    OpalMonitor.h
    OpalMarker.h
    OpalMultipole.h
    OpalMultipoleT.h
    OpalQuadrupole.h
    OpalRBend.h
    OpalSBend.h
    OpalSolenoid.h
    OpalTravelingWave.h
    OpalLaser.h
    OpalPolynomialTimeDependence.h
    OpalSplineTimeDependence.h
    OpalSinusoidalTimeDependence.h
    OpalVariableRFCavity.h
)

opalx_include_source_directory()

target_sources(opalx PRIVATE ${_INPUT_SRCS})
target_sources(opalx_run PRIVATE ${_RUN_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Elements")
