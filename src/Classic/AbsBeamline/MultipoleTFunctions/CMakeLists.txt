set (_SRCS
  CoordinateTransform.cpp
  DifferentialOperator.cpp
  DifferentialOperatorTwo.cpp
  Polynomial.cpp
  PolynomialSum.cpp
  RecursionRelation.cpp
  RecursionRelationTwo.cpp
  tanhDeriv.cpp
  TwoPolynomial.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
  CoordinateTransform.h
  DifferentialOperator.h
  DifferentialOperatorTwo.h
  Polynomial.h
  PolynomialSum.h
  RecursionRelation.h
  RecursionRelationTwo.h
  tanhDeriv.h
  TwoPolynomial.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/AbsBeamline/MultipoleTFunctions/")
#Turn on gcov
#set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)
#set(CMAKE_CXX_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
