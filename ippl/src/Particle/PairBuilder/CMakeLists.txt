SET (_SRCS
  )

SET (_HDRS
  BasicPairBuilder.h
  HashPairBuilder.h
  PairConditions.h
  SortingPairBuilder.h
  )

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_SOURCES (${_SRCS})
ADD_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/Particle/PairBuilder)
