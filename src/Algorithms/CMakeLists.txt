set (_SRCS
    AbstractTimeDependence.cpp
    AbstractTracker.cpp
    CavityAutophaser.cpp
    DefaultVisitor.cpp
    DistributionMoments.cpp
    Flagger.cpp
    GeneralElementMapping.cpp
    ExternalFieldRayTracker.cpp
    DeviceExternalField.cpp
    PassiveRingProbe.cpp
    OneTurnMap.cpp
    ClosedOrbitSolver.cpp
    LinearMapEigenAnalysis.cpp
    SpectralTunes.cpp
    ReferenceTracer.cpp
    ParallelTracker.cpp
    TrackReferenceStep.cpp
    PartData.cpp
    PolynomialTimeDependence.cpp
    SplineTimeDependence.cpp
    SinusoidalTimeDependence.cpp
    StepSizeConfig.cpp
    Tracker.cpp
    CoordinateSystemTrafo.cpp
    Quaternion.cpp
)

set (HDRS
    AbstractTimeDependence.h
    AbstractTracker.h
    BeamStop.h
    CavityAutophaser.h
    CompensatedSum.h
    RungeKuttaTableau.h
    DefaultVisitor.h
    DistributionMoments.h
    Flagger.h
    GeneralElementMapping.h
    ExternalFieldRayTracker.h
    DeviceExternalField.h
    PassiveProbe.h
    PassiveRingProbe.h
    OneTurnMap.h
    ClosedOrbitSolver.h
    LinearMapEigenAnalysis.h
    SpectralTunes.h
    ParallelTracker.h
    ReferenceTracer.h
    ClosedOrbitInitialState.h
    TrackReferenceStep.h
    PartData.h
    PolynomialTimeDependence.h
    SplineTimeDependence.h
    SinusoidalTimeDependence.h
    StepSizeConfig.h
    Tracker.h
    Matrix.h
    CoordinateSystemTrafo.h
    Quaternion.hpp
)

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")
