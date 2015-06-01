set (_SRCS
  Beam.cpp
  OpalWake.cpp
  SurfacePhysics.cpp
  SecondaryEmissionPhysics.cpp
  PriEmissionPhysics.cpp
  BoundaryGeometry.cpp
  FieldSolver.cpp
  DataSink.cpp
  H5PartWrapper.cpp
  H5PartWrapperForPC.cpp
  H5PartWrapperForPS.cpp
  H5PartWrapperForPT.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})