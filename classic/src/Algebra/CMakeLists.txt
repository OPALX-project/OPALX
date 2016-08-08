set (_SRCS
  ComplexEigen.cpp
  DoubleEigen.cpp
  DynamicFixedPoint.cpp
  NormalForm.cpp
  StaticFixedPoint.cpp
  TpsComplex.cpp
  TpsData.cpp
  TpsDouble.cpp
  TpsMonomial.cpp
  TpsTps.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_cl_sources(${_SRCS})