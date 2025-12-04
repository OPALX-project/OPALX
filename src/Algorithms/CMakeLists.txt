set (_SRCS
    AbstractTimeDependence.cpp
    AbstractTracker.cpp
    CavityAutophaser.cpp
    DefaultVisitor.cpp
    DistributionMoments.cpp
    Flagger.cpp
    IndexMap.cpp
    OrbitThreader.cpp
    OpalParticle.cpp
    ParallelTracker.cpp
    PartBins.cpp
    PartData.cpp
    PolynomialTimeDependence.cpp
    SplineTimeDependence.cpp
    StepSizeConfig.cpp
    Tracker.cpp
    CoordinateSystemTrafo.cpp
    Quaternion.cpp
)

set (HDRS
    AbstractTimeDependence.h
    AbstractTracker.h
    BoostMatrix.h
    CavityAutophaser.h
    DefaultVisitor.h
    DistributionMoments.h
    Flagger.h
    IndexMap.h
    OrbitThreader.h
    OpalParticle.h
    ParallelTracker.h
    PartBins.h
    PartBunch.hpp
    PartData.h
    PolynomialTimeDependence.h
    SplineTimeDependence.h
    StepSizeConfig.h
    Tracker.h
    Matrix.h
    CoordinateSystemTrafo.h
    Quaternion.hpp
)

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")