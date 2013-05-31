set (_SRCS
  Euclid3D.cpp
  Geometry.cpp
  Matrix3D.cpp
  NullGeometry.cpp
  OffsetGeometry.cpp
  PlanarArcGeometry.cpp
  RBendGeometry.cpp
  Rotation3D.cpp
  SRotatedGeometry.cpp
  StraightGeometry.cpp
  Vector3D.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_cl_sources(${_SRCS})