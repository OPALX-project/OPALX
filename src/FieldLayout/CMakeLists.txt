SET (_SRCS
  FieldLayoutUser.cpp
  MultiBalancer.cpp
  VRB.cpp
  )

SET (_HDRS
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

INCLUDE_DIRECTORIES (
  ${CMAKE_CURRENT_SOURCE_DIR}
)

ADD_SOURCES (${_SRCS})
ADD_HEADERS (${_HDRS})

INSTALL (FILES ${_HDRS} DESTINATION include/FieldLayout)
