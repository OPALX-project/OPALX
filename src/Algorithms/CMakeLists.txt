set (_SRCS
    AbstractTimeDependence.cpp
    BeamStop.cpp
    BunchReference.cpp
    CavityAutophaser.cpp
    GeneralElementMapping.cpp
    ExternalFieldRayTracker.cpp
    ReferenceTracer.cpp
    ParallelTracker.cpp
    TrackReferenceStep.cpp
    PartData.cpp
    PolynomialTimeDependence.cpp
    SplineTimeDependence.cpp
    SinusoidalTimeDependence.cpp
    StepSizeConfig.cpp
    CoordinateSystemTrafo.cpp
    Quaternion.cpp
)

set (HDRS
    AbstractTimeDependence.h
    BeamStop.h
    BunchReference.h
    CavityAutophaser.h
    CompensatedSum.h
    RungeKuttaTableau.h
    GeneralElementMapping.h
    ExternalFieldRayTracker.h
    ParallelTracker.h
    ReferenceTracer.h
    TrackReferenceStep.h
    PartData.h
    PolynomialTimeDependence.h
    SplineTimeDependence.h
    SinusoidalTimeDependence.h
    StepSizeConfig.h
    Matrix.h
    CoordinateSystemTrafo.h
    Quaternion.hpp
)

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")
