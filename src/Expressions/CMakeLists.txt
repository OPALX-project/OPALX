set (_SRCS
    ATable.cpp
    Expressions.cpp
    SFunction.cpp
    SHash.cpp
    )

set (HDRS
    AAutomatic.h
    ABinary.h
    ADeferred.h
    AList.h
    ARefExpr.h
    ASUnary.h
    ATable.h
    AUnary.h
    AValue.h
    Expressions.h
    Indexer.h
    SAutomatic.h
    SBinary.h
    SConstant.h
    SDeferred.h
    SFunction.h
    SHash.h
    SNull.h
    SRefAttr.h
    SRefExpr.h
    SUnary.h
    SValue.h
    TFind.h
    TFunction0.h
    TFunction1.h
    TFunction2.h
    )

opalx_include_source_directory()

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Expressions")