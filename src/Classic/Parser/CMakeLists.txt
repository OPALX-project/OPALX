set (_SRCS
  AbsFileStream.cpp
  ClassicParser.cpp
  FileStream.cpp
  Parser.cpp
  SimpleStatement.cpp
  Statement.cpp
  StringStream.cpp
  TerminalStream.cpp
  Token.cpp
  TokenStream.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

ADD_OPAL_SOURCES(${_SRCS})

set (HDRS
    AbsFileStream.h
    ClassicParser.h
    FileStream.h
    Parser.h
    SimpleStatement.h
    Statement.h
    StringStream.h
    TerminalStream.h
    Token.h
    TokenStream.h
)

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Parser")