set (_SRCS
  AbstractMapper.cpp
  AbstractTimeDependence.cpp
  AbstractTracker.cpp
  CoordinateSystemTrafo.cpp
  DefaultVisitor.cpp
  Flagger.cpp
  IdealMapper.cpp
  LinearMapper.cpp
  MapIntegrator.cpp
  Mapper.cpp
  OrbitTracker.cpp
  PartBunch.cpp
  PartBins.cpp
  PartBinsCyc.cpp
  PartData.cpp
  OpalParticle.cpp
  PolynomialTimeDependence.cpp
  Surveyor.cpp
  ThinMapper.cpp
  ThinTracker.cpp
  TrackIntegrator.cpp
  Tracker.cpp
  Quaternion.cpp
  )

if ( ENABLE_AMR )
    list(APPEND _SRCS AmrPartBunch.cpp)
endif ( ENABLE_AMR )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    AbstractMapper.h
    AbstractTimeDependence.h
    AbstractTracker.h
    CoordinateSystemTrafo.h
    DefaultVisitor.h
    Flagger.h
    IdealMapper.h
    LinearMapper.h
    ListElem.h
    MapIntegrator.h
    Mapper.h
    OrbitTracker.h
    PartBinsCyc.h
    PartBins.h
    PartBunch.h
    PartBunchBase.h
    PartData.h
    OpalParticle.h
    PBunchDefs.h
    PolynomialTimeDependence.h
    Quaternion.h
    rbendmap.h
    Surveyor.h
    ThinMapper.h
    ThinTracker.h
    Tracker.h
    TrackIntegrator.h
    Vektor.h
)

if ( ENABLE_AMR )
    list(APPEND HDRS AmrPartBunch.h)
endif ( ENABLE_AMR )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Algorithms")