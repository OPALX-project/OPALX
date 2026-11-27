set (_SRCS
    AbstractTimeDependence.cpp
    BeamStop.cpp
    BunchReference.cpp
    CavityAutophaser.cpp
    GeneralElementMapping.cpp
    ExternalFieldRayTracker.cpp
    ReferenceTracer.cpp
    ParallelTracker.cpp
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
    GeneralElementMapping.h
    ExternalFieldRayTracker.h
    ParallelTracker.h
    ReferenceTracer.h
    PartData.h
    PolynomialTimeDependence.h
    SplineTimeDependence.h
    SinusoidalTimeDependence.h
    StepSizeConfig.h
    Matrix.h
    CoordinateSystemTrafo.h
    Quaternion.hpp
)

opalx_include_source_directory()

target_sources(opalx_run PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")
