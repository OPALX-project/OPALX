set (_SRCS
    BeamBeamWindowAnimation.cpp
    ArithmeticError.cpp
    CLRangeError.cpp
    ExpressionRandom.cpp
    ComplexErrorFun.cpp
    ConvergenceError.cpp
    DivideError.cpp
    DomainError.cpp
    FormatError.cpp
    GeneralOpalException.cpp
    LogicalError.cpp
    Options.cpp
    OverflowError.cpp
    ParseError.cpp
    PortableGraymapReader.cpp
    Util.cpp
    SDDSParser.cpp
    OpalException.cpp
    OpalFilter.cpp
    RegularExpression.cpp
    Timer.cpp
    AbstractSpline.cpp
    LinearSpline.cpp
    CubicSpline.cpp
    )

set (HDRS
    BeamBeamWindowAnimation.h
    ArithmeticError.h
    CLRangeError.h
    ElementList.h
    ExpressionRandom.h
    ComplexErrorFun.h
    ConvergenceError.h
    DivideError.h
    DomainError.h
    FormatError.h
    GeneralOpalException.h
    LogicalError.h
    Options.h
    OverflowError.h
    ParseError.h
    PortableGraymapReader.h
    Util.h
    OpalException.h
    OpalFilter.h
    RegularExpression.h
    Timer.h
    AbstractSpline.h
    LinearSpline.h
    CubicSpline.h
    GSLSpline.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Utilities")
