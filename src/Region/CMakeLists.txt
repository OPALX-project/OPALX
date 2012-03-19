set (_SRCS
  )

set (_HDRS
  NDRegion.h
  PRegion.h
  RegionLayout.cpp
  RegionLayout.h
  Rnode.h
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
add_headers(${_HDRS})

install(FILES ${_HDRS} DESTINATION include/Region)
