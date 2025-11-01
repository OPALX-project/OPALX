set (_SRCS
    Attributes.cpp
    Bool.cpp
    BoolArray.cpp
    Place.cpp
    PredefinedString.cpp
    Range.cpp
    Real.cpp
    RealArray.cpp
    Reference.cpp
    String.cpp
    StringArray.cpp
    TableRow.cpp
    TokenList.cpp
    TokenListArray.cpp
    UpperCaseString.cpp
    UpperCaseStringArray.cpp
    )

set (HDRS
    Attributes.h
    BoolArray.h
    Bool.h
    Place.h
    PredefinedString.h
    Range.h
    RealArray.h
    Real.h
    Reference.h
    String.h
    StringArray.h
    TableRow.h
    TokenListArray.h
    TokenList.h
    UpperCaseString.h
    UpperCaseStringArray.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Attributes")