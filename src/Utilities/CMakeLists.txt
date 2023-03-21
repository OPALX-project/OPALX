set (_SRCS
    SDDSParser.cpp
    SDDSParser/array.cpp
    SDDSParser/associate.cpp
    SDDSParser/ast.cpp
    SDDSParser/column.cpp
    SDDSParser/data.cpp
    SDDSParser/description.cpp
    SDDSParser/file.cpp
    SDDSParser/include.cpp
    SDDSParser/parameter.cpp
    SDDSParser/version.cpp
#
    EarlyLeaveException.cpp
    OpalException.cpp
    OpalFilter.cpp
    RegularExpression.cpp
    Timer.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    EarlyLeaveException.h
    OpalException.h
    OpalFilter.h
    RegularExpression.h
    Timer.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Utilities")
