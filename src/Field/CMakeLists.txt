SET (_SRCS
  )

SET (_HDRS
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

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_SOURCES (${_SRCS})
ADD_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/Field)
