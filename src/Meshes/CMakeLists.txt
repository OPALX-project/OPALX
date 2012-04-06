SET (_SRCS
  Centering.cpp
  )

SET (_HDRS
  CartesianCentering.cpp
  CartesianCentering.h
  Cartesian.cpp
  Cartesian.h
  CartesianStencilSetup.h
  Centering.h
  Mesh.cpp
  Mesh.h
  UniformCartesian.cpp
  UniformCartesian.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_SOURCES (${_SRCS})
ADD_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/Meshes)
