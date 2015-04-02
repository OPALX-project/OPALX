SET (_SRCS
  )

SET (_HDRS
  SubBareField.hpp
  SubBareField.h
  SubFieldAssign.hpp
  SubFieldAssignDefs.h
  SubFieldAssign.h
  SubField.hpp
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
