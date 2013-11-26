set (_SRCS
  ArithmeticError.cpp
  AttributeError.cpp
  ClassicException.cpp
  ComplexErrorFun.cpp
  ConvergenceError.cpp
  DivideError.cpp
  DomainError.cpp
  EigenvalueError.cpp
  FormatError.cpp
  Gauss.cpp
  InverseGauss.cpp
  LogicalError.cpp
  NormalFormError.cpp
  OverflowError.cpp
  ParseError.cpp
  Random.cpp
  CLRangeError.cpp
  SingularMatrixError.cpp
  SizeError.cpp
  SwitcherError.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_cl_sources(${_SRCS})
