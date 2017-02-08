set (_SRCS
  ArithmeticError.cpp
  AttributeError.cpp
  ClassicException.cpp
  ClassicField.cpp
  CLRangeError.cpp
  ComplexErrorFun.cpp
  ConvergenceError.cpp
  DivideError.cpp
  DomainError.cpp
  EigenvalueError.cpp
  FormatError.cpp
  Gauss.cpp
  GeneralClassicException.cpp
  InverseGauss.cpp
  LogicalError.cpp
  NormalFormError.cpp
  Options.cpp
  OverflowError.cpp
  ParseError.cpp
  ClassicRandom.cpp
  RingSection.cpp
  SingularMatrixError.cpp
  SizeError.cpp
  SwitcherError.cpp
  Util.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    ArithmeticError.h
    AttributeError.h
    ClassicException.h
    ClassicField.h
    ClassicRandom.h
    CLRangeError.h
    ComplexErrorFun.h
    ConvergenceError.h
    DivideError.h
    DomainError.h
    EigenvalueError.h
    FormatError.h
    Gauss.h
    GeneralClassicException.h
    InverseGauss.h
    LogicalError.h
    NormalFormError.h
    Options.h
    OptionTypes.h
    OverflowError.h
    ParseError.h
    RingSection.h
    SingularMatrixError.h
    SizeError.h
    SwitcherError.h
    Util.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Utilities")