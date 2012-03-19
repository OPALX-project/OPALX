set (_SRCS
  )

set (_HDRS
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

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
add_headers(${_HDRS})

install(FILES ${_HDRS} DESTINATION include/SubField)
