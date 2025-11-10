set (_SRCS
    Euclid3D.cpp
    Geometry.cpp
    Matrix3D.cpp
    Euclid3DGeometry.cpp
    NullGeometry.cpp
    PlanarArcGeometry.cpp
    VarRadiusGeometry.cpp
    RBendGeometry.cpp
    Rotation3D.cpp
    StraightGeometry.cpp
    Vector3D.cpp
    )

set (HDRS
    Euclid3DGeometry.h
    Euclid3D.h
    Geometry.h
    Matrix3D.h
    NullGeometry.h
    PlanarArcGeometry.h
    VarRadiusGeometry.h
    RBendGeometry.h
    Rotation3D.h
    StraightGeometry.h
    Vector3D.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/BeamlineGeometry")