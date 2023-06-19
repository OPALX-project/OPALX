set (_SRCS
    AbstractTimeDependence.cpp
    AbstractTracker.cpp
    CavityAutophaser.cpp
    DefaultVisitor.cpp
    Flagger.cpp
    IndexMap.cpp
    OrbitThreader.cpp
    ParallelTTracker.cpp
    PartBunch.cpp
    PartBunchBase.hpp
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
    CavityAutophaser.h
    DefaultVisitor.h
    Flagger.h
    IndexMap.h
    OrbitThreader.h
    ParallelTTracker.h
    PartBins.h
    PartBunch.h
    PartBunchBase.h
    PartData.h
    PolynomialTimeDependence.h
    SplineTimeDependence.h
    StepSizeConfig.h
    Tracker.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")
