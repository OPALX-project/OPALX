set (_SRCS
  Interpolator3dGridTo1d.cpp
  Interpolator3dGridTo3d.cpp
  Mesh.cpp
  MMatrix.cpp
  MVector.cpp
  PolynomialPatch.cpp
  PPSolveFactory.cpp
  SolveFactory.cpp
  SquarePolynomialVector.cpp
  ThreeDGrid.cpp
  TriLinearInterpolator.cpp
  VectorMap.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
