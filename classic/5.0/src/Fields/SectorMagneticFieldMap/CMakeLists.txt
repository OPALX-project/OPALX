set (_SRCS
  Interpolator3dGridTo1d.cpp
  Interpolator3dGridTo3d.cpp
  Mesh.cpp
  SectorField.cpp
  ThreeDGrid.cpp
  TriLinearInterpolator.cpp
  VectorMap.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_cl_sources(${_SRCS})
