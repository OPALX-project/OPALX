set (_SRCS
  AbstractMapper.cpp
  AbstractTimeDependence.cpp
  AbstractTracker.cpp
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
  Particle.cpp
  PolynomialTimeDependence.cpp
  Surveyor.cpp
  ThinMapper.cpp
  ThinTracker.cpp
  TrackIntegrator.cpp
  Tracker.cpp
  Vektor.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_sources(${_SRCS})
