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

opalx_include_source_directory()

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Attributes")
