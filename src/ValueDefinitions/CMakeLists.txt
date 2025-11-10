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

target_include_directories(opalx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(opalx PRIVATE ${_SRCS})

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ValueDefinitions")