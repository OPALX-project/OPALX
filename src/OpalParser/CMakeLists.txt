set (_SRCS
  CompoundStatement.cpp
  IfStatement.cpp
  Macro.cpp
  MacroCmd.cpp
  MacroStream.cpp
  OpalParser.cpp
  WhileStatement.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
