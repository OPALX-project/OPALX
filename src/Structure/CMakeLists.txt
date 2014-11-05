set (_SRCS
  Beam.cpp
  OpalWake.cpp
  SurfacePhysics.cpp
  SecondaryEmissionPhysics.cpp
  PriEmissionPhysics.cpp
  BoundaryGeometry.cpp
  FieldSolver.cpp
  DataSink.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
