set (_SRCS
  AColumn.cpp
  ARow.cpp
  ATable.cpp
  Expressions.cpp
  SCell.cpp
  SFunction.cpp
  SHash.cpp
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
