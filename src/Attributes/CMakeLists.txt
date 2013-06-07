set (_SRCS
  Attributes.cpp
  Bool.cpp
  BoolArray.cpp
  Place.cpp
  Range.cpp
  Real.cpp
  RealArray.cpp
  Reference.cpp
  opalstr.cpp
  StringArray.cpp
  TableRow.cpp
  TokenList.cpp
  TokenListArray.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
