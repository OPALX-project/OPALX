# Input side: the OPTION globals and the parser's error.
set (_INPUT_SRCS
    Options.cpp
    ParseError.cpp
    )

# Run side.
set (_RUN_SRCS
    ArithmeticError.cpp
    CLRangeError.cpp
    ExpressionRandom.cpp
    DomainError.cpp
    FormatError.cpp
    GeneralOpalException.cpp
    LogicalError.cpp
    OverflowError.cpp
    Util.cpp
    OpalException.cpp
    Timer.cpp
    AbstractSpline.cpp
    LinearSpline.cpp
    CubicSpline.cpp
    )

set (HDRS
    ArithmeticError.h
    CLRangeError.h
    ElementList.h
    ExpressionRandom.h
    DomainError.h
    FormatError.h
    GeneralOpalException.h
    LogicalError.h
    Options.h
    OverflowError.h
    ParseError.h
    Util.h
    OpalException.h
    Timer.h
    AbstractSpline.h
    LinearSpline.h
    CubicSpline.h
    GSLSpline.h
    )

opalx_include_source_directory()

target_sources(opalx PRIVATE ${_INPUT_SRCS})
target_sources(opalx_run PRIVATE ${_RUN_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Utilities")
