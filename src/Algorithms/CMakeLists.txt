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
    PartBunch.cpp
#    PartBunch.hpp
    PartBins.cpp
    PartData.cpp
    PolynomialTimeDependence.cpp
    SplineTimeDependence.cpp
    StepSizeConfig.cpp
    Tracker.cpp
)   

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

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
    PartBunch.h
 #   PartBunch.h
    PartData.h
    PolynomialTimeDependence.h
    SplineTimeDependence.h
    StepSizeConfig.h
    Tracker.h
    Matrix.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")
