set (_SRCS
  Euclid3D.cpp
  Geometry.cpp
  Matrix3D.cpp
  Euclid3DGeometry.cpp
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

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    Euclid3DGeometry.h
    Euclid3D.h
    Geometry.h
    Matrix3D.h
    NullGeometry.h
    OffsetGeometry.h
    PlanarArcGeometry.h
    RBendGeometry.h
    Rotation3D.h
    SRotatedGeometry.h
    StraightGeometry.h
    Vector3D.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/BeamlineGeometry")