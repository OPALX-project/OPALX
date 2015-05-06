SET (_SRCS
  )

SET (_HDRS
  Assign.hpp
  AssignDefs.h
  AssignGeneralBF.hpp
  AssignGeneralIBF.hpp
  Assign.h
  AssignTags.h
  BareField.hpp
  BareField.h
  BareFieldIterator.h
  BCond.hpp
  BCond.h
  BrickExpression.hpp
  BrickExpression.h
  BrickIterator.hpp
  BrickIterator.h
  CompressedBrickIterator.hpp
  CompressedBrickIterator.h
  Field.hpp
  Field.h
  FieldLoc.h
  FieldSpec.h
  GuardCellSizes.hpp
  GuardCellSizes.h
  IndexedBareField.hpp
  IndexedBareField.h
  IndexedField.h
  LField.hpp
  LField.h
  ReductionLoc.hpp
  ReductionLoc.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_SOURCES (${_SRCS})
ADD_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/Field)
