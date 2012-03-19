set (_SRCS
  )

set (_HDRS
  Assign.cpp
  AssignDefs.h
  AssignGeneralBF.cpp
  AssignGeneralIBF.cpp
  Assign.h
  AssignTags.h
  BareField.cpp
  BareField.h
  BareFieldIterator.h
  BCond.cpp
  BCond.h
  BrickExpression.cpp
  BrickExpression.h
  BrickIterator.cpp
  BrickIterator.h
  CompressedBrickIterator.cpp
  CompressedBrickIterator.h
  Field.cpp
  Field.h
  FieldLoc.h
  FieldSpec.h
  GuardCellSizes.cpp
  GuardCellSizes.h
  IndexedBareField.cpp
  IndexedBareField.h
  IndexedField.h
  LField.cpp
  LField.h
  ReductionLoc.cpp
  ReductionLoc.h
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
add_headers(${_HDRS})

install(FILES ${_HDRS} DESTINATION include/Field)
