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

add_sources(${_SRCS})