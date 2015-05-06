SET (_SRCS
  Centering.cpp
  )

SET (_HDRS
  CartesianCentering.hpp
  CartesianCentering.h
  Cartesian.hpp
  Cartesian.h
  CartesianStencilSetup.h
  Centering.h
  Mesh.hpp
  Mesh.h
  UniformCartesian.hpp
  UniformCartesian.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_SOURCES (${_SRCS})
ADD_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/Meshes)
