set (_SRCS
  Action.cpp
  Attribute.cpp
  AttributeBase.cpp
  AttributeHandler.cpp
  BeamSequence.cpp
  Definition.cpp
  Directory.cpp
  Editor.cpp
  Element.cpp
  Invalidator.cpp
  OpalData.cpp
  Object.cpp
  ObjectFunction.cpp
  PlaceRep.cpp
  RangeRep.cpp
  Table.cpp
  TableRowRep.cpp
  ValueDefinition.cpp
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
