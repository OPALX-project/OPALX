set (_SRCS
  )

set (_HDRS
  IpplExpressions.h
  IpplTypeComputations.h
  PETE.h
  TypeComputations.h
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
add_headers(${_HDRS})

install(FILES ${_HDRS} DESTINATION include/PETE)
