set (_SRCS
  FlatWriter.cpp
  Line.cpp
  LineTemplate.cpp
  Replacer.cpp
  Sequence.cpp
  SequenceMember.cpp
  SequenceParser.cpp
  SequenceTemplate.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
