set (_SRCS
    Attributes.cpp
    Bool.cpp
    PredefinedString.cpp
    Range.cpp
    Real.cpp
    RealArray.cpp
    OpalString.cpp
    StringArray.cpp
    UpperCaseString.cpp
    UpperCaseStringArray.cpp
    )

set (HDRS
    Attributes.h
    Bool.h
    PredefinedString.h
    Range.h
    RealArray.h
    Real.h
    OpalString.h
    StringArray.h
    UpperCaseString.h
    UpperCaseStringArray.h
    )

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Attributes")
