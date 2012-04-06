SET (_SRCS
  )

SET (_HDRS
  SubBareField.cpp
  SubBareField.h
  SubFieldAssign.cpp
  SubFieldAssignDefs.h
  SubFieldAssign.h
  SubField.cpp
  SubField.h
  SubFieldIter.h
  SubFieldTraits.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_SOURCES (${_SRCS})
ADD_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/SubField)
