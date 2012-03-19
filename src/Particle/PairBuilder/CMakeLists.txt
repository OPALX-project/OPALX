set (_SRCS
  )

set (_HDRS
  BasicPairBuilder.h
  HashPairBuilder.h
  PairConditions.h
  SortingPairBuilder.h
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
add_headers(${_HDRS})

install(FILES ${_HDRS} DESTINATION include/Particle/PairBuilder)
