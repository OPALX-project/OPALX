set (_SRCS
    AbstractTimeDependence.cpp
    AbstractTracker.cpp
    CavityAutophaser.cpp
    DefaultVisitor.cpp
    Flagger.cpp
    IndexMap.cpp
    ExternalFieldRayTracker.cpp
    DeviceExternalField.cpp
    PassiveRingProbe.cpp
    CyclotronProbeRecorder.cpp
    OneTurnMap.cpp
    ClosedOrbitSolver.cpp
    MatchedCovariance.cpp
    MatchedDistributionCovariance.cpp
    SpaceChargeMap.cpp
    LinearMapEigenAnalysis.cpp
    TpsClosedOrbitMap.cpp
    TpsLinearMapAlgebra.cpp
    TpsLinearTransferMap.cpp
    LinearTransferMapBuilder.cpp
    SpectralTunes.cpp
    OrbitThreader.cpp
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
    CavityAutophaser.h
    CompensatedSum.h
    RungeKuttaTableau.h
    DefaultVisitor.h
    Flagger.h
    IndexMap.h
    ExternalFieldRayTracker.h
    DeviceExternalField.h
    PassiveProbe.h
    CyclotronProbe.h
    CyclotronProbeRecorder.h
    PassiveRingProbe.h
    OneTurnMap.h
    ClosedOrbitSolver.h
    LinearMapEigenAnalysis.h
    TpsClosedOrbitMap.h
    TpsLinearMapAlgebra.h
    TpsLinearTransferMap.h
    LinearTransferMapBuilder.h
    SpectralTunes.h
    OrbitThreader.h
    OrbitThreaderDiagnostics.h
    ParallelTracker.h
    ClosedOrbitInitialState.h
    ClosedOrbitOptics.h
    MatchedCovariance.h
    MatchedDistributionCovariance.h
    SpaceChargeMap.h
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
