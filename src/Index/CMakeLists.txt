set (_SRCS
  Index.cpp
  )

set (_HDRS
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

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
add_headers(${_HDRS})

install(FILES ${_HDRS} DESTINATION include/Index)
