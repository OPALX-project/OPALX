set (_SRCS
    MacroCmd.cpp
    AbsFileStream.cpp
    CompoundStatement.cpp
    FileStream.cpp
    SimpleStatement.cpp
    Statement.cpp
    StringStream.cpp
    Token.cpp
    TokenStream.cpp
    OpalParser.cpp
    IfStatement.cpp	
    WhileStatement.cpp
    Macro.cpp
    MacroStream.cpp         
)

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    MacroCmd.h    
    AbsFileStream.h
    CompoundStatement.h
    FileStream.h
    Parser.h
    SimpleStatement.h
    Statement.h
    StringStream.h
    Token.h
    TokenStream.h
    OpalParser.h
    IfStatement.h
    WhileStatement.h
    Macro.h
    MacroStream.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Parser")