SET (_SRCS
  Index.cpp
  )

SET (_HDRS
  IndexedSIndex.h
  Index.h
  IndexInlines.h
  LSIndex.h
  NDIndex.h
  NDIndexInlines.h
  SIndexAssign.cpp
  SIndexAssign.h
  SIndex.cpp
  SIndex.h
  SOffset.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_SOURCES (${_SRCS})
ADD_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/Index)
