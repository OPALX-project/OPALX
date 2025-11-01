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

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/AbsBeamline/MultipoleTFunctions/")