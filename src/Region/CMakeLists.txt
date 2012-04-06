SET (_SRCS
  )

SET (_HDRS
  NDRegion.h
  PRegion.h
  RegionLayout.cpp
  RegionLayout.h
  Rnode.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_SOURCES (${_SRCS})
ADD_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/Region)
