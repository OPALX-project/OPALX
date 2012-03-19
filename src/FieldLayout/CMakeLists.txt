set (_SRCS
  FieldLayoutUser.cpp
  MultiBalancer.cpp
  VRB.cpp
  )

set (_HDRS
  BinaryBalancer.cpp
  BinaryBalancer.h
  CenteredFieldLayout.cpp
  CenteredFieldLayout.h
  ConejoBalancer.cpp
  ConejoBalancer.h
  FieldLayout.cpp
  FieldLayout.h
  FieldLayoutUser.h
  MultiBalancer.h
  Vnode.h
  VnodeMultiBalancer.cpp
  VnodeMultiBalancer.h
  VRB.h
  )

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

add_sources(${_SRCS})
add_headers(${_HDRS})

install(FILES ${_HDRS} DESTINATION include/FieldLayout)
