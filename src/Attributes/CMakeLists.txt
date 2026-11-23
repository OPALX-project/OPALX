set (_SRCS
    Attributes.cpp
    Bool.cpp
    BoolArray.cpp
    PredefinedString.cpp
    Range.cpp
    Real.cpp
    RealArray.cpp
    Reference.cpp
    OpalString.cpp
    StringArray.cpp
    TokenList.cpp
    TokenListArray.cpp
    UpperCaseString.cpp
    UpperCaseStringArray.cpp
    )

set (HDRS
    Attributes.h
    BoolArray.h
    Bool.h
    PredefinedString.h
    Range.h
    RealArray.h
    Real.h
    Reference.h
    OpalString.h
    StringArray.h
    TokenListArray.h
    TokenList.h
    UpperCaseString.h
    UpperCaseStringArray.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Attributes")
