set (_SRCS
    BoolConstant.cpp
    RealConstant.cpp
    RealVariable.cpp
    RealVector.cpp
    StringConstant.cpp
    )

set (HDRS
    BoolConstant.h
    RealConstant.h
    RealVariable.h
    RealVector.h
    StringConstant.h
    )

opalx_include_source_directory()

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ValueDefinitions")