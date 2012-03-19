set (_SRCS
  Centering.cpp
  )

set (_HDRS
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

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
add_headers(${_HDRS})

install(FILES ${_HDRS} DESTINATION include/Meshes)
