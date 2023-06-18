set (_SRCS
    AbstractTimeDependence.cpp
    AbstractTracker.cpp
    CoordinateSystemTrafo.cpp
    DefaultVisitor.cpp
    DistributionMoments.cpp
    Flagger.cpp
    PartBunch.cpp
    PartBunchBase.hpp
    PartBins.cpp
    PartData.cpp
    OpalParticle.cpp
    PolynomialTimeDependence.cpp
    SplineTimeDependence.cpp
    Tracker.cpp
    Quaternion.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    AbstractTimeDependence.h
    AbstractTracker.h
    CoordinateSystemTrafo.h
    DefaultVisitor.h
    DistributionMoments.h
    Flagger.h
    ListElem.h
    PartBins.h
    PartBunch.h
    PartBunchBase.h
    PartData.h
    OpalParticle.h
    PBunchDefs.h
    PolynomialTimeDependence.h
    Quaternion.h
    SplineTimeDependence.h
    Tracker.h
    Vektor.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")
